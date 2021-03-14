clear all
close all
clc

DEBUG = 0;
spacing = 'constant';
type = 'WT';
[polar, prop, oper, air, propellertype] = param(spacing, type);
%% Results from Carlos' code
zeroyawresults = load('zeroyawdata.mat');
radialdist_zeroyaw = zeroyawresults.SectionResults(:,1);
a_zeroyaw = zeroyawresults.SectionResults(:,2);
aprime_zeroyaw = zeroyawresults.SectionResults(:,3);
%%
yaw = deg2rad([15]); 
dPsi = 1;%[deg]
AzimuthAngle = deg2rad([0:dPsi:360]');

%% polar plots for airfoil
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
%Allocate
SectionResults = zeros(length(prop.r_R)-1,5);%r_R|average a|average aprime|iteration for each r_R|section rotor solidity
AveAxialInductions = zeros(length(prop.r_R)-1, length(yaw));
AveTangInductions = zeros(length(prop.r_R)-1, length(yaw));
ResultsForDifferentYaw = cell(length(yaw), 1);


for i=1:length(prop.r_R)-1
    prop.sectionchord(i,1) = chord_distribution(0.5*(prop.r_R(i)+prop.r_R(i+1)), prop.R, propellertype);%non-dimensional
    prop.sectionpitch(i,1) = pitch_distribution(0.5*(prop.r_R(i)+prop.r_R(i+1)),prop.collective_blade_twist, propellertype);% [deg] 
end
%Solve for each element/annulus
for j=1:length(yaw)
    for i=1:length(prop.r_R)-1
        [SectionResults(i,:), ~, Cx(i,:),~] = GetAveFactors.SolveSection(i, polar, prop, air, oper, AzimuthAngle, yaw(j), dPsi);
    end
    ResultsForDifferentYaw(j,1) = {SectionResults};
    SectionRadialDist = SectionResults(:,1); % non-dimensional
    SectionRotorSolidity = SectionResults(:,5); 
    AveAxialInductions(:,j) = SectionResults(:,2);
    AveTangInductions(:,j) = SectionResults(:,3);
    yawindeg = rad2deg(yaw(j));
    Finish = ['Run for ',num2str(yawindeg),' degrees finished.'];
    disp(Finish);
end
figure;
plot(SectionRadialDist,AveAxialInductions(:,1),'k');
% hold on
% plot(SectionRadialDist,AveAxialInductions(:,2),'b--');
% % plot(SectionRadialDist,AveAxialInductions(:,3),'r-.');
% legend('a(yaw=0)','a (yaw=15deg)');
% legend('a(yaw=0)','a (yaw=15deg)','a (yaw=30deg)')


figure;
plot(SectionRadialDist,AveTangInductions(:,1),'k');
hold on
% plot(SectionRadialDist,AveTangInductions(:,2),'b--');
% plot(SectionRadialDist,AveTangInductions(:,3),'r-.');
% legend('a''(yaw=0)','a'' (yaw=15deg)');

%% calculate axial induction factor that varies with azimuth based on average tangential induction
%length(SectionRadialDist)
for i=1:length(SectionRadialDist)
    [a_Azim(i,:), iterations(i,1)] = SolveElement(prop, oper, polar, SectionRadialDist(i), ...
        SectionRotorSolidity(i), AveTangInductions(i), i, AzimuthAngle, yaw);
    Finish = ['Radius index: ' num2str(i)];
    disp(Finish);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%Post Processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions
%solve for each element using average aprime
function [InflowAngles, alphas, cl, cd] = AnnularRingProperties(yaw, a, aprime, prop, oper, SectionRadius, AzimuthAngle, polar, index)
    [InflowAngles] = inflowangle(yaw, a, aprime, prop, oper, SectionRadius, AzimuthAngle);%should be an array of phi for constant radius
    alphas = rad2deg(InflowAngles) - prop.sectionpitch(index);
    cl = interp1(polar.alpha, polar.Cl, alphas,'linear', 'extrap');
    cd = interp1(polar.alpha, polar.Cd, alphas,'linear', 'extrap');
end
function [phi] = inflowangle(yaw, a, aprime, prop, oper, SectionRadius, AzimuthAngle)
    %need to double check
    K_c = K(yaw, a);%constant for constant yaw, 0 for 0 deg yaw
    chi = WakeSkewAngle(yaw, a);%constant for constant yaw, 0 for 0 deg yaw
    r_R = SectionRadius./prop.R;%constant
    F = FlowExpansionFunction(r_R);%constant for constant radius
    temp1 = a.*(1+F.*K_c.*sin(AzimuthAngle));
    top = oper.U_inf.*(cos(yaw)-temp1)+(oper.omega.*SectionRadius.*aprime.*cos(AzimuthAngle).*sin(chi).*(1+sin(AzimuthAngle).*sin(chi)));
    temp1 = oper.omega.*SectionRadius.*(1+aprime.*cos(chi).*(1+sin(AzimuthAngle).*sin(chi)));
    temp2 = a.*tan(chi/2).*(1+F.*K_c.*sin(AzimuthAngle)) - sin(yaw);
    btm = temp1 + oper.U_inf.*cos(AzimuthAngle).*temp2;
    phi = atan(top./btm);
end
function [a_final, iterations] = SolveElement...
    (prop, oper, polar, r_R, sigma, aprime, index, AzimuthAngle, yaw)
    % r_R: SectionRadius [non-dimensional], single value
    % sigma: SectionRotorSolidity, single value
    % aprime: average aprime, single value
    % AzimuthAngles: vector of azimuth angles
    % index: keep track of radial position
    % yaw: single value
    
    SectionRadius = r_R.*prop.R; %[m]
    
    N = 100;%number of iterations for convergence
    epsilon = 0.0001;%tolerance
    a_final = zeros(1, length(AzimuthAngle));
    
%     for j=1:length(AzimuthAngle)
        %initialising for each Azimiuth Angle (solving for a only)
        a = 0.1.*ones(length(AzimuthAngle),1);%axial induction factor
        for i=1:N
            [InflowAngle, alphas, cl, cd] = AnnularRingProperties(yaw, a, aprime, prop, oper, SectionRadius, AzimuthAngle, polar, index);        
%             chi = WakeSkewAngle(yaw, a);

            W2 = ResultantVelocityBladeElement(SectionRadius, oper, AzimuthAngle, yaw, a, aprime);%should be array of W here
            Cx = cl.*cos(InflowAngle)+cd.*sin(InflowAngle);
%             Cy = cl.*sin(InflowAngle)-cd.*cos(InflowAngle);
            [~, ~, Prandtl]=PrandtlTipRootCorrection(prop, SectionRadius, oper, a);
            CtElmAzim = sigma.*(W2./(oper.U_inf^2)).*Cx;
            
            % Only 1 equation required since there is only 1 unknown
            % solve for new a using fsolve
            % determine momentum eqn to be used based on value of a
            a_new = solveCoeffFunc(a, Prandtl, yaw, CtElmAzim);

            if (max(abs(a-a_new))<epsilon)
                break;
            else
                a = 0.75*a+0.25*a_new;            
            end
        end
        a_final(1,:) = a_new;
        iterations = i;
%         Finish = ['Azimnth index: ' num2str(j)];
%         disp(Finish);
%     end
end 
function [ftip, froot, ftotal] = PrandtlTipRootCorrection(prop, SectionRadius, oper, a)
    r_R_Section = SectionRadius/prop.R; %non-dimensional
    F1 = (-prop.Nblades./2).*((1-r_R_Section)./r_R_Section).*sqrt(1+(oper.TSR.*r_R_Section./(1-a)).^2);
    F2 = (-prop.Nblades./2).*((r_R_Section-prop.blade_root)./r_R_Section).*sqrt(1+(oper.TSR.*r_R_Section./(1-a)).^2);
    
%     F1 = (-prop.Nblades/2).*((1-r_R_Section)./r_R_Section).*(1./sin(InflowAngle));
%     F2 = (-prop.Nblades/2).*((r_R_Section-prop.blade_root)./r_R_Section).*(1./sin(InflowAngle));
    
    ftip = (2./pi).*acos(exp(F1));    
%     if (isnan(ftip))
%        ftip = 0; 
%     end
    ftip(isnan(ftip)) = 0;
    froot = (2./pi).*acos(exp(F2));
%     if (isnan(froot))
%        froot = 0;
%     end
    froot(isnan(froot)) = 0;
    ftotal = ftip.*froot;
    ftotal(ftotal<0.0001) = 0.0001;
    
end
function Uresultant2 = ResultantVelocityBladeElement(SectionRadius, oper, AzimuthAngle, yaw, a, aprime)
    chi = WakeSkewAngle(yaw, a);
    Vx = oper.U_inf;
    Vy = oper.omega*SectionRadius;
    %1st term
    temp1_1 = Vx.*(cos(yaw)-a);
    temp1_2 = Vy.*aprime.*sin(chi).*cos(AzimuthAngle).*(1+sin(AzimuthAngle).*sin(chi));
    temp1 =  temp1_1 + temp1_2;
    
    %2nd term
    temp2_1 =  Vy.*(1+aprime.*cos(chi).*(1+sin(AzimuthAngle).*sin(chi)));
    temp2_2 = Vx.*cos(AzimuthAngle).*(a.*tan(chi./2)-sin(yaw));
    temp2 =  temp2_1 + temp2_2;
    Uresultant2 = (temp1.^2)+(temp2.^2);
end
function CtFunc = CoeffMomentumFuncGlauert(a_new, F, yaw, CtRHS)    
    %momentum theory from Glauert
   
    tmp= 2*cos(yaw) - a_new*F;
    CtFunc= CtRHS-(4*a_new*F*sqrt(1-a_new*F*tmp));
      
       
end
function CtFunc = CoeffMomentumFuncPropBrake(a_new, F, yaw, CtRHS)    
     
    %use different Ct_momentum functions for different a   
    
    %propeller-brake
    CtFunc= CtRHS-(4*a_new*F*(a_new*F-1)*cos(yaw));
            
end
function CtFunc = CoeffMomentumFuncEmpirical(a_new, F, yaw, CtRHS)    
    %use different Ct_momentum functions for different a   
    beta = 0.4;

    [f0, fprime_0, f1]=fcoeff(F, yaw);
    k2 = double((f1-f0-fprime_0*(1-beta))/((1-beta)^2));
    k1 = double(fprime_0-2*k2*beta);
    k0 = double(f1-k1-k2);
        
    CtFunc= CtRHS-(k0+(k1*a_new)+(k2*(a_new^2)));
               
end
function VectorFunc = CoeffMomentumVectorFunc(a_new, F, yaw, a_test, CtRHS)
    beta = 0.4;
    for i=1:length(CtRHS)
        if a_test(i)<=beta
            VectorFunc(i) = CtRHS(i)-(4.*a_new(i).*F(i).*...
                sqrt(1-a_new(i).*F(i).*(2*cos(yaw)-a_new(i)*F(i))));
        elseif a_test(i)>1            
            VectorFunc(i) = CtRHS(i)-(4.*a_new(i).*F(i).*(a_new(i).*F(i)-1).*cos(yaw));
        else
            [f0, fprime_0, f1]=fcoeff(F(i), yaw);
            k2 = double((f1-f0-fprime_0*(1-beta))/((1-beta)^2));
            k1 = double(fprime_0-2*k2*beta);
            k0 = double(f1-k1-k2);

            VectorFunc(i) = CtRHS(i)-(k0+(k1.*a_new(i))+(k2.*(a_new(i).^2)));  
        end
    end
end
function  [f0, fprime_0, f1] = fcoeff(F, yaw)
    beta = 0.4;
    syms x
    g = 4*x*F*sqrt(1-x*F*(2*cos(yaw) - x*F));
    y = diff(g);
    f0 = double(vpa(subs(g,x,beta)));
    fprime_0=double(vpa(subs(y,x,beta)));
    f1 = double(2*cos(yaw));
end
function a_new = solveCoeffFunc(a_previous, F, yaw, CtRHS)
    options = optimoptions('fsolve','Display','iter');
    options.Algorithm = 'trust-region-reflective';
    options.JacobPattern = speye(length(CtRHS));
    options.PrecondBandWidth = 0;
    initial = 0.1*ones(length(CtRHS),1);
%     beta = 0.4;
    %create vector of functions for one rev: output vector of a_new
    
    a_new = fsolve(@(a_new) CoeffMomentumVectorFunc(a_new, F, yaw, a_previous, CtRHS), initial, options);

%     %-----
%     if a_previous<=beta
% %         disp('Glauert');
%         a_new = fsolve(@(a_new) CoeffMomentumFuncGlauert(a_new, F, yaw, CtRHS), initial, options);
%     elseif a_previous>1
% %         disp('propBrake');
%         a_new = fsolve(@(a_new) CoeffMomentumFuncPropBrake(a_new, F, yaw, CtRHS), initial, options);
%     else
% %         disp('Empirical');
%         a_new = fsolve(@(a_new) CoeffMomentumFuncEmpirical(a_new, F, yaw, CtRHS), initial, options);
%     end
end    
function F = FlowExpansionFunction(r_R)
    %curve fit to original function
    F = 0.5.*(r_R+0.4.*(r_R.^3)+0.4.*(r_R.^5));
end
function chi = WakeSkewAngle(yaw, a)
    %yaw: [rad]
    chi = yaw.*(0.6.*a+1);
end
function K_c = K(yaw, a)
    chi = WakeSkewAngle(yaw, a);
    K_c = 2.*tan(chi./2);
end
