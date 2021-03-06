clear all
close all
clc

DEBUG = 0;
spacing = 'constant';
type = 'WT';
[polar, prop, oper, air, propellertype] = param(spacing, type);

% %------------
% yaw = [15];
% AzimuthAngle = deg2rad([0:1:360]);

% %-------------------
% figure
% subplot(1,2,1);
% plot(polar.alpha,polar.Cl);
% xlabel('\alpha');
% ylabel('C_{l}');
% 
% subplot(1,2,2);
% plot(polar.alpha,polar.Cd);
% xlabel('\alpha');
% ylabel('C_{d}');
%%
SectionResults = zeros(length(prop.r_R)-1,8);
for i=1:length(prop.r_R)-1
    SectionResults(i,:) = SolveSection(i, polar, prop, air, oper);
end
% 
dT = SectionResults(:,4).*prop.Nblades.*(prop.dr.*prop.R);
% CT = sum(dT./(0.5*air.density*(oper.U_inf^2)*(pi*prop.R^2)))
CT = sum(dT./(0.5*(oper.U_inf^2)*(pi*prop.R^2)))

dP = SectionResults(:,5).*prop.Nblades.*(prop.dr.*prop.R).*SectionResults(:,1).*(oper.omega*prop.R);
% CP = sum(dP./(0.5*air.density*(oper.U_inf^3)*(pi*prop.R^2)))
CP = sum(dP./(0.5*(oper.U_inf^3)*(pi*prop.R^2)))

%% Functions
function a = AxialInductionFactor(CT)
    CT1 = 1.816;
    CT2 = 2*sqrt(CT1)-CT1;
    %Glauert correction for heavily loaded rotors
    if (CT<CT2)
        a = 0.5 - 0.5*sqrt(1-CT);        
    else
        top = CT-CT1;
        btm = 4*sqrt(CT1)-4;
        a = 1+(top/btm);
    end
end
function aprime = AzimInductionFactor(Fazim, prop, oper, a, air, index)    
    top = Fazim * prop.Nblades;
    btm = 2*air.density*pi*oper.U_inf*(1-a)*oper.omega*2*prop.SectionRadius(index)^2;    
    aprime = top/btm;
end
function [ftip, froot, ftotal] = PrandtlTipRootCorrection(prop, oper, a, index)
    r_R_Section = prop.SectionRadius(index)/prop.R; %non-dimensional
    F1 = (-prop.Nblades/2)*((1-r_R_Section)/r_R_Section).*sqrt(1+(oper.TSR.*r_R_Section./(1-a)).^2);
    F2 = (-prop.Nblades/2)*((r_R_Section-prop.blade_root)/r_R_Section).*sqrt(1+(oper.TSR.*r_R_Section./(1-a)).^2);
    
%     F1 = (-prop.Nblades/2)*((1-r_R_Section)/r_R_Section)*(1/sin(InflowAngle));
%     F2 = (-prop.Nblades/2)*((r_R_Section-prop.blade_root)/r_R_Section)*(1/sin(InflowAngle));
    
    ftip = (2/pi).*acos(exp(F1));    
    ftip(isnan(ftip)) = 0;
    
    froot = (2/pi).*acos(exp(F2));
    froot(isnan(froot)) = 0;
    
    ftotal = ftip.*froot;
    ftotal(ftotal< 0.0001) = 0.0001;
    
end
function [f_axial, f_tangential, gamma, alpha, U_resultant] = SectionLoads(U_axial, U_tangential, prop, polar, index)
    InflowAngle = atan2(U_axial, U_tangential); % [rad]
    alpha = rad2deg(InflowAngle) - prop.sectionpitch(index); 
    U_resultant = norm([U_axial, U_tangential]);
    cl = interp1(polar.alpha, polar.Cl, alpha, 'linear', 'extrap');
    cd = interp1(polar.alpha, polar.Cd, alpha, 'linear', 'extrap');
%     lift = 0.5*air.density*(U_resultant^2)*(prop.sectionchord*prop.R)*cl;
%     drag = 0.5*air.density*(U_resultant^2)*(prop.sectionchord*prop.R)*cd;
    
    lift = 0.5*(U_resultant^2)*(prop.sectionchord(index)*prop.R)*cl;
    drag = 0.5*(U_resultant^2)*(prop.sectionchord(index)*prop.R)*cd;
    f_axial = lift*cos(InflowAngle)+drag*sin(InflowAngle);
    f_tangential = lift*sin(InflowAngle)-drag*cos(InflowAngle);
%     gamma = lift/(air.density*U_resultant);      
    gamma = lift/(U_resultant);      

end
function Results = SolveSection(index, polar, prop, air, oper, ~, ~)
%     r_R1=prop.r_R(index); %non-dimensional
%     r_R2=prop.r_R(index+1);
%     SectionArea = pi*(((r_R2*prop.R)^2)-((r_R1*prop.R)^2)); %[m^2]
%     SectionRadius = 0.5*(r_R1+r_R2)*prop.R;%average radius [m]

%     yaw = deg2rad(yaw);
%     AzimuthAngle = deg2rad(AzimuthAngle);
    
    a = 0.1;%axial induction factor
    aprime = 0.1;%tangential induction factor
    
    N = 100;%number of iterations
    epsilon = 0.00001;
    
    Ux = oper.U_inf;
    Uy = oper.omega*prop.SectionRadius(index);
    for i=1:N
        U_axial =Ux*(1-a);
        U_tangential = Uy.*(1+aprime);       
        [f_axial, f_tangential, gamma, alpha, InflowAngle]=SectionLoads(U_axial, U_tangential,prop, polar, index);
        %Thrust coefficient at streamtube        
        load3D_axial = f_axial*prop.Nblades*(prop.dr(index)*prop.R);
%         CT_streamtube= load3D_axial./(0.5*air.density*(oper.U_inf^2)*SectionArea);
        CT_streamtube= load3D_axial./(0.5*(oper.U_inf^2)*prop.SectionArea(index));
        
        %---------------TO EDIT---------------------
        a_new=AxialInductionFactor(CT_streamtube);
        
        [~, ~, Prandtl]=PrandtlTipRootCorrection(prop, oper, a_new, index);
        a_new = a_new./Prandtl;
        a = 0.75*a+0.25*a_new;%update new value of axial induction via weighted average
        
        aprime = AzimInductionFactor(f_tangential, prop, oper, a, air, index);
        aprime = aprime./Prandtl;
        %---------------------------------------------------------------------------
        if (abs(a-a_new)<epsilon)
            break;
        end       
    end 
    
    
    Results = [prop.SectionRadius(index)/prop.R, a, aprime, f_axial, f_tangential,gamma, alpha,  i];
end
function skewangle = SkewWakeAngle(yaw, a)
    skewangle = yaw.*(0.6.*a+1);
end
function askew = SkewWakeCorrection(yaw, a, AzimuthAngle)
    skewangle = SkewWakeAngle(yaw,a);
    tmp = (15*pi/64).*tan(skewangle/2).*r_R.*sin(AzimuthAngle);
    askew = a.*(1+tmp);
end