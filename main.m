clear all
close all
clc

%Airfoil polar
data = readmatrix('ARAD8pct_polar.txt');
polar.alpha = data(:,1);
polar.Cl = data(:,2);
polar.Cd = data(:,3);
polar.Cm = data(:,4);

figure
subplot(1,2,1);
plot(polar.alpha,polar.Cl);
xlabel('\alpha');
ylabel('C_{l}');

subplot(1,2,2);
plot(polar.alpha,polar.Cd);
xlabel('\alpha');
ylabel('C_{d}');


%Propeller geometry
prop.blade_root = 0.25; %non-dimensional blade starting location 
spacing = 'constant';
switch spacing
    case 'constant'
        %constant spacing
        prop.r_R = linspace(prop.blade_root,1.00,100);%constant spacing
        prop.r_R = prop.r_R(:);%make column vector
        prop.dr = prop.r_R(2:end)-prop.r_R(1:end-1); %non-dimensional length of each section
    case 'cosine'
        theta = linspace(0,pi,100);
        prop.r_R = prop.blade_root+0.5*(1-prop.blade_root).*(1-cos(theta));
        prop.r_R = prop.r_R(:);%make column vector
        prop.dr = prop.r_R(2:end)-prop.r_R(1:end-1);
    otherwise 
        disp('Error: Invalid spacing specified');
        return
end
prop.c_R = chord_distribution(prop.r_R);
prop.twist = -50*prop.r_R+35; %local twist
prop.collective_blade_twist = 46;% [degrees]
prop.pitch = pitch_distribution(prop.r_R, prop.collective_blade_twist);
prop.R = 0.70;%propeller radius [m]
prop.Nblades = 6;%number of blades

%Operating conditions
oper.U_inf = 60;% [m/s]
oper.n = 1200;% rotational speed [rpm]
oper.omega = (oper.n*2*pi)/(60); %rotational speed [rad/s]
oper.TSR = (oper.omega*prop.R)/oper.U_inf;

[air.Temp, air.speed_of_sound, air.pressure, air.density] = atmosisa(2000);

SectionResults = zeros(length(prop.r_R)-1,6);
for i=1:1
    chord = chord_distribution(0.5*(prop.r_R(i)+prop.r_R(i+1)));%non-dimensional
    pitch = pitch_distribution(chord,prop.collective_blade_twist);% [deg]    
    [AnnulusRadius, a, aprime, f_axial, f_tangential, gamma,Prandtltip,Prandtlroot,alpha] = SolveSection(i, polar, prop, air, oper, chord, pitch);
end

% dT = SectionResults(:,4).*prop.Nblades.*prop.dr.*prop.R;
% CT = sum(dT./(0.5*air.density*(oper.U_inf^2)*);length(prop.r_R)-



%% Functions
function a = AxialInductionFactor(CT)
    CT1 = 1.816;
    CT2 = 2*sqrt(CT1)-CT1;
    %Glauert correction for heavily loaded rotors
    if (CT<CT2)
        a = 0.5-0.5*sqrt(1-CT);        
    else
        top = CT-CT1;
        btm = 4*sqrt(CT1)-4;
        a = 1+(top/btm);
    end
end

function [ftip, froot, ftotal] = PrandtlTipRootCorrection(TSR, r_R, a, B, root)
    A1 = sqrt(1+((TSR*r_R)^2)/((1+a)^2));
    F1 = exp(-0.5*B*((1-r_R)/r_R)*A1);   
    F2 = exp(-0.5*B*((r_R-root)/r_R)*A1);
    ftip = (2/pi)*acos(F1);
    if (isnan(ftip))
       ftip = 0; 
    end
    froot = (2/pi)*acos(F2);
    if (isnan(froot))
       froot = 0;
    end
    
    ftotal = ftip*froot;
    if (ftotal < 0.0001)
        ftotal = 0.0001;
    end
end

function [f_axial, f_tangential, gamma,alpha] = SectionLoads(U_axial, U_tangential, twist, chord, polar, air)
    InflowAngle = atan2(U_axial, U_tangential); % [deg]
    alpha = 1*(twist - rad2deg(InflowAngle)); %need to double check formula
    U_resultant = sqrt((U_axial^2)+(U_tangential^2));
    cl = interp1(polar.alpha, polar.Cl, alpha);
    cd = interp1(polar.alpha, polar.Cd, alpha);
    lift = 0.5*air.density*(U_resultant^2)*chord*cl;
    drag = 0.5*air.density*(U_resultant^2)*chord*cd;
    f_axial = lift*cos(InflowAngle)-drag*sin(InflowAngle);
    f_tangential = lift*sin(InflowAngle)+drag*cos(InflowAngle);
    gamma = lift/(air.density*U_resultant);      
end
function aprime = AzimInductionFactor(Fazim, B, rho, r, oper, a)
    top = Fazim * B;
    btm = 2*rho*(2*pi*r)*oper.U_inf*(1+a)*oper.omega*r;
    aprime = top/btm;
end
function [AnnulusRadius, a, aprime, f_axial, f_tangential, gamma, Prandtltip,Prandtlroot, alpha]=SolveSection(index, polar, prop, air, oper, chord, twist)
    r_R1=prop.r_R(index); %non-dimensional
    r_R2=prop.r_R(index+1);
    AnnulusArea = pi*(prop.R^2)*((r_R2^2)-(r_R1^2)); %[m^2]
    AnnulusRadius = 0.5*(r_R1+r_R2)*prop.R;%taken as average, [m]
    
    a = 0;%axial induction factor
    aprime = 0.0;%tangential induction factor
    
    N = 1;%number of iterations
    epsilon = 0.0001;
    
    for i=1:N
        U_axial = oper.U_inf*(1+a);
        U_tangential = oper.omega*AnnulusRadius*(1-aprime);       
        [f_axial, f_tangential, gamma,alpha]=SectionLoads(U_axial, U_tangential, twist, chord*prop.R, polar, air);
        %Thrust coefficient at streamtube        
        load3D_axial = f_axial*prop.Nblades*prop.dr(index)*prop.R;
        CT_streamtube= load3D_axial./(0.5*air.density*(oper.U_inf^2)*AnnulusArea);
        a_new=AxialInductionFactor(CT_streamtube);
        
        [Prandtltip, Prandtlroot, Prandtltotal]=PrandtlTipRootCorrection(oper.TSR, AnnulusRadius/prop.R, a_new, prop.Nblades, prop.blade_root);
        a_new = a_new/Prandtltotal;
        a = 0.75*a+0.25*a_new;%update new value of axial induction via weighted average
        
        aprime = AzimInductionFactor(f_tangential, prop.Nblades, air.density, AnnulusRadius, oper, a);
        aprime = aprime/Prandtltotal;
        
        if (abs(a-a_new)<epsilon)
            break;
        end       
    end  
end
function c_R = chord_distribution(r_R)
    c_R = 0.18 - 0.06.*r_R;
end

function localpitch = pitch_distribution(r_R,collective_blade_twist)
    twist = -50.*r_R + 35; %local twist [deg]
    localpitch = twist + collective_blade_twist;
end

