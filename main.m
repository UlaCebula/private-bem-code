clear all
close all
clc


spacing = 'constant';
[polar, prop, oper, air] = param(spacing);

figure
subplot(1,2,1);
plot(polar.alpha,polar.Cl);
xlabel('\alpha');
ylabel('C_{l}');

subplot(1,2,2);
plot(polar.alpha,polar.Cd);
xlabel('\alpha');
ylabel('C_{d}');
%%
SectionResults = zeros(length(prop.r_R)-1,6);
for i=1:1
    prop.chord = chord_distribution(0.5*(prop.r_R(i)+prop.r_R(i+1)));%non-dimensional
    prop.pitch = pitch_distribution(prop.chord,prop.collective_blade_twist);% [deg]    
    [alpha, a, aprime, f_axial, f_tangential, gamma,Prandtl] = SolveSection(i, polar, prop, air, oper);
end

% dT = SectionResults(:,4).*prop.Nblades.*prop.dr.*prop.R;
% CT = sum(dT./(0.5*air.density*(oper.U_inf^2)*);length(prop.r_R)-



%% Functions
function a = AxialInductionFactor(CT)
    CT1 = 1.816;
    CT2 = 2*sqrt(CT1)-CT1;
    %Glauert correction for heavily loaded rotors
    if (CT<CT2)
        a = -0.5+0.5*sqrt(1+CT);        
    else
        top = CT1-CT;
        btm = 4*(sqrt(CT1)-1);
        a = -1+(top/btm);
    end
end
function aprime = AzimInductionFactor(Fazim, prop, air, AnnulusRadius, oper, a)    
    top = Fazim * prop.Nblades;
    btm = 2*air.density*(2*pi*AnnulusRadius^2)*oper.U_inf*(1+a)*oper.omega;    
    aprime = top/btm;
end
function [ftip, froot, ftotal] = PrandtlTipRootCorrection(prop, AnnulusRadius, InflowAngle)
    r_R_Annulus = AnnulusRadius/prop.R; %non-dimensional
    F1 = (-prop.Nblades/2)*((1-r_R_Annulus)/r_R_Annulus)*(1/abs(sin(InflowAngle)));
    F2 = (-prop.Nblades/2)*((r_R_Annulus-prop.blade_root)/r_R_Annulus)*(1/abs(sin(InflowAngle)));
    
    ftip = (2/pi)*acos(exp(F1));    
    if (isnan(ftip))
       ftip = 0; 
    end
    
    froot = (2/pi)*acos(exp(F2));
    if (isnan(froot))
       froot = 0;
    end
    
    ftotal = ftip*froot;
    if (ftotal < 0.0001)
        ftotal = 0.0001;
    end
end

function [f_axial, f_tangential, gamma,InflowAngle, alpha] = SectionLoads(U_axial, U_tangential, prop, polar, air)
    InflowAngle = atan2(U_axial, U_tangential); % [rad]
    alpha = prop.pitch - rad2deg(InflowAngle); 
    U_resultant = norm([U_axial, U_tangential]);
    cl = interp1(polar.alpha, polar.Cl, alpha);
    cd = interp1(polar.alpha, polar.Cd, alpha);
    lift = 0.5*air.density*(U_resultant^2)*(prop.chord*prop.R)*cl;
    drag = 0.5*air.density*(U_resultant^2)*(prop.chord*prop.R)*cd;
    f_axial = lift*cos(InflowAngle)-drag*sin(InflowAngle);
    f_tangential = lift*sin(InflowAngle)+drag*cos(InflowAngle);
    gamma = lift/(air.density*U_resultant);      
end

function [alpha, a, aprime, f_axial, f_tangential, gamma, Prandtl]=SolveSection(index, polar, prop, air, oper)
    r_R1=prop.r_R(index); %non-dimensional
    r_R2=prop.r_R(index+1);
    AnnulusArea = pi*(((r_R2*prop.R)^2)-((r_R1*prop.R)^2)); %[m^2]
    AnnulusRadius = 0.5*(r_R1+r_R2)*prop.R;%average radius [m]
    
    a = 0.3;%axial induction factor
    aprime = 0;%tangential induction factor
    
    N = 1;%number of iterations
    epsilon = 0.0001;
    
    for i=1:N
        U_axial = oper.U_inf*(1+a);
        U_tangential = oper.omega*AnnulusRadius*(1-aprime);       
        [f_axial, f_tangential, gamma,InflowAngle, alpha]=SectionLoads(U_axial, U_tangential,prop, polar, air);
        %Thrust coefficient at streamtube        
        load3D_axial = f_axial*prop.Nblades*(prop.dr(index)*prop.R);
        CT_streamtube= load3D_axial./(0.5*air.density*(oper.U_inf^2)*AnnulusArea);
        a_new=AxialInductionFactor(CT_streamtube);
        
        [~, ~, Prandtl]=PrandtlTipRootCorrection(prop, AnnulusRadius, InflowAngle);
        a_new = a_new/Prandtl;
        a = 0.75*a+0.25*a_new;%update new value of axial induction via weighted average
        
        aprime = AzimInductionFactor(f_tangential, prop, air, AnnulusRadius, oper, a);
        aprime = aprime/Prandtl;
%         aprime = 0.75*aprime+0.25*aprime_new;    
        if (abs(a-a_new)<epsilon)
            break;
        end       
    end  
end
