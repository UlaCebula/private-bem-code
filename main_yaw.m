clear all
close all
clc

DEBUG = 0;
spacing = 'constant';
type = 'WT';
[polar, prop, oper, air, propellertype] = param(spacing, type);

zeroyawresults = load('zeroyawdata.mat');
radialdist_zeroyaw = zeroyawresults.SectionResults(:,1);
a_zeroyaw = zeroyawresults.SectionResults(:,2);
aprime_zeroyaw = zeroyawresults.SectionResults(:,3);
%%%
yaw = deg2rad([15]); 
dPsi = 1;%[deg]
AzimuthAngle = deg2rad([0:dPsi:360]');

%%
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
% SectionResults = zeros(length(prop.r_R)-1,4);
for j=1:length(yaw)
    for i=1:length(prop.r_R)-1
        prop.sectionchord = chord_distribution(0.5*(prop.r_R(i)+prop.r_R(i+1)), prop.R, propellertype);%non-dimensional
        prop.sectionpitch = pitch_distribution(0.5*(prop.r_R(i)+prop.r_R(i+1)),prop.collective_blade_twist, propellertype);% [deg]        
        [SectionResults(i,:,j), InflowAngles(i,:,j), Cx(i,:,j),Cy(i,:,j)] = GetAveFactors.SolveSection(i, polar, prop, air, oper, AzimuthAngle, yaw(j), dPsi);
    end
end
%% calculate axial induction factor that varies with azimuth based on average tangential induction




%%
%%%%%%%%%%%%%%%%%%%%%%%%%Post Processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% plot(SectionResults(:,1,1),SectionResults(:,2,1),'k');
% hold on
% plot(radialdist_zeroyaw, a_zeroyaw,'k--');
% plot(SectionResults(:,1,2),SectionResults(:,2,2),'b');
% plot(SectionResults(:,1,3),SectionResults(:,2,3),'r');
% legend('coupled(yaw=0)','ref data for yaw=0','coupled(yaw=15deg','coupled(yaw=30deg')
%% Functions

%% Function handles for solving for variation of axial induction with azimuth (not in use)
function phi = inflowangleFunc(yaw, a, aprime, prop, oper, SectionRadius)
    K_c = K(yaw, a);%constant for constant yaw, 0 for 0 deg yaw
    chi = WakeSkewAngle(yaw, a);%constant for constant yaw, 0 for 0 deg yaw
    r_R = SectionRadius./prop.R;
    F = FlowExpansionFunction(r_R);%constant for constant radius
    Vx = oper.U_inf;
    Vy = oper.omega*SectionRadius;
    
    toptemp1 = @(AzimuthAngle) Vx.*(cos(yaw)-a.*(1+F.*K_c.*sin(AzimuthAngle)));
    toptemp2 = @(AzimuthAngle) Vy.*aprime.*cos(AzimuthAngle).*sin(chi).*(1+sin(AzimuthAngle).*sin(chi));
    top = @(AzimuthAngle) toptemp1(AzimuthAngle) + toptemp2(AzimuthAngle);
    
    
    bottemp1 = @(AzimuthAngle) Vy.*(1+aprime.*cos(chi).*(1+sin(AzimuthAngle).*sin(chi)));
    bottemp2 = @(AzimuthAngle) (Vx.*cos(AzimuthAngle)).*(a.*tan(chi/2).*(1+F.*K_c.*sin(AzimuthAngle)) - sin(yaw));
    btm = @(AzimuthAngle) bottemp1(AzimuthAngle) + bottemp2(AzimuthAngle);
    phi = @(AzimuthAngle) atan(top(AzimuthAngle)./btm(AzimuthAngle));
end
function [cx, cy] = ForceCoeffFunc(cl, cd, phi)
    cx = @(AzimuthAngle) cl.*cos(phi(AzimuthAngle)) + cd.*sin(phi(AzimuthAngle));
    cy = @(AzimuthAngle) cl.*sin(phi(AzimuthAngle)) - cd.*cos(phi(AzimuthAngle));
end
function Uresultant2 = ResultantVelocityBladeElementFunc(SectionRadius, oper, yaw, a, aprime)
    chi = WakeSkewAngle(yaw, a);
    Vx = oper.U_inf;
    Vy = oper.omega*SectionRadius;
    %1st term
    temp1_1 = @(AzimuthAngle) Vx*(cos(yaw)-a);
    temp1_2 = @(AzimuthAngle) Vy.*aprime.*sin(chi).*cos(AzimuthAngle).*(1+sin(AzimuthAngle).*sin(chi));
    temp1 = @(AzimuthAngle) temp1_1(AzimuthAngle) + temp1_2(AzimuthAngle);
    
    %2nd term
    temp2_1 = @(AzimuthAngle) Vy.*(1+aprime.*cos(chi).*(1+sin(AzimuthAngle).*sin(chi)));
    temp2_2 = @(AzimuthAngle) Vx*cos(AzimuthAngle).*(a*tan(chi/2)-sin(yaw));
    temp2 = @(AzimuthAngle) temp2_1(AzimuthAngle) + temp2_2(AzimuthAngle);
    Uresultant2 = @(AzimuthAngle) (temp1(AzimuthAngle).^2)+(temp2(AzimuthAngle).^2);
end
function CTelement = CTelementFunc(W2, sigmasection, cx, oper)
    U2 = oper.U_inf^2;
    CTelement = @(AzimuthAngle) (W2(AzimuthAngle)./U2).*cx(AzimuthAngle).*sigmasection;
end
function CQelement = CQelementFunc(W2, sigmasection, cx, cy, oper, chi)
    U2 = oper.U_inf^2;
    temp = @(AzimuthAngle) cy(AzimuthAngle).*cos(chi)-cx(AzimuthAngle).*sin(chi).*cos(AzimuthAngle);
    CQelement = @(AzimuthAngle) (W2(AzimuthAngle)./U2).*sigmasection.*(temp(AzimuthAngle));
end
%% old functions (not in use)
% function [a] = CalculateNewAxialInduction(W, SectionRotorSolidity, Cx, dPsi, Prandtl, yaw, chi, oper)
%     %solve integral discretely
%     NormalisedW = (W.^2)/(oper.U_inf^2); 
%     temp = SectionRotorSolidity.*sum(Cx.*NormalisedW.*dPsi);
%     %solve for a
%     fun = @(a) 8*pi*a*Prandtl*(cos(yaw)+(tan(chi/2)*sin(yaw))-a*Prandtl*(sec(chi/2))^2)-temp;
% %     fun = @(a) 	8*pi*a*Prandtl*sqrt(1-a*Prandtl*((2*cos(yaw))-(a*Prandtl)))-temp;
%     a = fzero(fun, 0.);    
% end
% function aprime = CalculateNewAzimInduction(W, SectionRotorSolidity, Cx, Cy, dPsi, AzimuthAngle, chi, yaw, Prandtl, a, oper, SectionRadius, prop)    
%     NormalisedW = (W.^2)./(oper.U_inf^2); 
%     r_R = SectionRadius/prop.R;
%     temp1 = sum(NormalisedW.*(cos(AzimuthAngle).*sin(chi).*Cx)+(cos(chi).*Cy).*dPsi);
%     temp = SectionRotorSolidity.*temp1;
%     
% %     fun = @(aprime) (4*aprime*Prandtl*(cos(yaw)-(a*Prandtl))*oper.TSR*r_R*pi*(1+(cos(chi)).^2))-temp;
%     fun = @(aprime) 4*aprime*Prandtl*oper.TSR*r_R*pi*(cos(yaw)-(a*Prandtl))*(1+(cos(chi))^2)-temp;
% 
%     aprime = fzero(fun, 0);  
% 
% end
