clear all
close all
clc

%Based on NREL methodology

DEBUG = 0;
spacing = 'constant';
type = 'WT';
[polar, prop, oper, air, propellertype] = param(spacing, type);

zeroyawresults = load('zeroyawdata.mat');
a_zeroyaw = zeroyawresults.SectionResults(:,2);
aprime_zeroyaw = zeroyawresults.SectionResults(:,3);
%%%
yaw = deg2rad(0); 
dPsi = 1;%[deg]
AzimuthAngle = deg2rad([0:dPsi:360]');
%%%

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
% SectionResults = zeros(length(prop.r_R)-1,8);
for i=1:length(prop.r_R)-1
%     for j=1:length(AzimuthAngle)%for each angle
        prop.sectionchord = chord_distribution(0.5*(prop.r_R(i)+prop.r_R(i+1)), prop.R, propellertype);%non-dimensional
        prop.sectionpitch = pitch_distribution(0.5*(prop.r_R(i)+prop.r_R(i+1)),prop.collective_blade_twist, propellertype);% [deg]        
        [SectionResults(i,:), a(i,:), aprime(i,:), alpha(i,:), inflowangle(i,:)] = SolveSection(i, polar, prop, air, oper, AzimuthAngle, dPsi, yaw);
%     end
end

% figure;
% plot(SectionResults(:,1), a,'k--');
% hold on
% plot(SectionResults(:,1),a_zeroyaw,'k');
% plot(SectionResults(:,1),aprime,'b--');
% plot(SectionResults(:,1),aprime_zeroyaw,'b');
%% Functions
function [ftip, froot, ftotal] = PrandtlTipRootCorrection(prop, SectionRadius, InflowAngle)
    r_R_Section = SectionRadius/prop.R; %non-dimensional
    
    %(prop, SectionRadius, oper, a)
%     F1 = (-prop.Nblades/2)*((1-r_R_Section)/r_R_Section)*sqrt(1+(oper.TSR*r_R_Section/(1-a))^2);
%     F2 = (-prop.Nblades/2)*((r_R_Section-prop.blade_root)/r_R_Section)*sqrt(1+(oper.TSR*r_R_Section/(1-a))^2);
    
    %(prop, SectionRadius, InflowAngle)
    F1 = (-prop.Nblades./2).*((1-r_R_Section)./r_R_Section).*(1./sin(InflowAngle));
    F2 = (-prop.Nblades./2).*((r_R_Section-prop.blade_root)./prop.blade_root).*(1./sin(InflowAngle));
    
    ftip = (2/pi).*acos(exp(F1));    
    if (isnan(ftip))
       ftip = 0; 
    end
    
    froot = (2/pi).*acos(exp(F2));
    if (isnan(froot))
       froot = 0;
    end
    
    ftotal = ftip.*froot;
%     if (ftotal < 0.0001)
%         ftotal = 0.0001;
%     end
    ftotal(ftotal<0.0001) = 0.0001;
end

function [Results, a_new, aprime_new, alpha, inflowangle] = SolveSection(index, polar, prop, ~, oper, AzimuthAngle, dPsi, yaw)
    dPsi = deg2rad(dPsi);
    r_R1=prop.r_R(index); %non-dimensional
    r_R2=prop.r_R(index+1);
%     SectionArea = pi*(((r_R2*prop.R)^2)-((r_R1*prop.R)^2)); %[m^2]
    SectionRadius = 0.5*(r_R1+r_R2)*prop.R;%average radius [m]
    SectionRotorSolidity = (prop.Nblades*prop.sectionchord*prop.R)/(2*pi*SectionRadius);

    %initialising
    a = 0.1;%axial induction factor
    aprime = 0.1;%tangential induction factor
    
    N = 100;%number of iterations for convergence
    epsilon = 0.00001;
    
    %calculate a and aprime for each section at a constant radius for each azimuth angle 
    for i=1:N
        %produce array of local AoA, Inflow angles, cl and cd
        [inflowangle, alpha, cx, cy] = SectionInflowAngleCalc(oper, a, aprime, SectionRadius, polar, prop, yaw, AzimuthAngle);
        [~, ~, Prandtl] = PrandtlTipRootCorrection(prop, SectionRadius, inflowangle);
        CTsection = CTCalc(SectionRotorSolidity, cx, a, inflowangle); 
        a_new = AxialInductionFactor(Prandtl, CTsection, SectionRotorSolidity, cx, inflowangle);
        aprime_new = AzimInductionFactor(Prandtl, SectionRotorSolidity, cy, inflowangle);       

        %Effect of skew
        wakeskewangle = WakeSkewAngleCalc(yaw, a_new);
        a_new = SkewedWakeCorrection(a_new,SectionRadius,prop,wakeskewangle,AzimuthAngle);

        if (abs(a-a_new)<epsilon)
            break;
        else
            a = 0.75.*a+0.25.*a_new;%,update new value of axial induction via weighted average
        end
    end 
    Results = [SectionRadius/prop.R, i];
end
function [phi, alpha, cx, cy] = SectionInflowAngleCalc(oper, a, aprime, SectionRadius, polar, prop, yaw, AzimuthAngle)
    U_axial = oper.U_inf.*(cos(yaw)).*(1-a);
    U_tan = (oper.U_inf.*(-sin(yaw).*cos(AzimuthAngle))+(oper.omega.*SectionRadius)).*(1+aprime);
    phi = atan2(U_axial,U_tan);
    alpha = rad2deg(phi) - prop.sectionpitch;    
    cl = interp1(polar.alpha, polar.Cl, alpha,'linear', 'extrap');
    cd = interp1(polar.alpha, polar.Cd, alpha,'linear', 'extrap');
    cx = cl.*cos(phi) + cd.*sin(phi);
    cy = cl.*sin(phi) - cd.*cos(phi);
end
function [Ct] = CTCalc(SectionRotorSolidity, cx, a, phi)
    top = SectionRotorSolidity.*cx.*(1-a).^2;
    btm = (sin(phi)).^2;
    Ct = top./btm;  
end
function a = AxialInductionFactor(F, CT, sigma, cx, phi)
    if (CT>0.96*F)
        temp1 = CT.*(50-36.*F)+12.*F.*(3.*F-4);
        top = 18.*F-20-3.*sqrt(temp1);
        btm = 36.*F-50;
        a = top./btm;
    else
        btm = sigma.*cx;
        top = 4.*F.*(sin(phi)).^2;        
        a = (1+(top./btm)).^-1;
    end
end
function aprime = AzimInductionFactor(F, sigma, cy, phi)
    top = 4.*F.*sin(phi).*cos(phi);
    btm = sigma.*cy;
    aprime = (-1+(top./btm)).^-1;
end
function wakeskewangle = WakeSkewAngleCalc(yaw, a)
    %yaw: [rad]
    wakeskewangle = yaw.*(0.6.*a+1);
end
function a_skew = SkewedWakeCorrection(a,SectionRadius,prop,chi,AzimuthAngle)
    r_R = SectionRadius./prop.R;
    temp = (15*pi/64).*r_R.*tan(chi/2).*cos(AzimuthAngle);
    a_skew = a.*(1+temp);
end