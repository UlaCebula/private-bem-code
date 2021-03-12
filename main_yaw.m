clear all
close all
clc

DEBUG = 0;
spacing = 'constant';
type = 'WT';
[polar, prop, oper, air, propellertype] = param(spacing, type);

zeroyawresults = load('zeroyawdata.mat');
a_zeroyaw = zeroyawresults.SectionResults(:,2);
aprime_zeroyaw = zeroyawresults.SectionResults(:,3);
%%%
yaw = deg2rad(15); 
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
% SectionResults = zeros(length(prop.r_R)-1,8);
for i=1:length(prop.r_R)-1
    prop.sectionchord = chord_distribution(0.5*(prop.r_R(i)+prop.r_R(i+1)), prop.R, propellertype);%non-dimensional
    prop.sectionpitch = pitch_distribution(0.5*(prop.r_R(i)+prop.r_R(i+1)),prop.collective_blade_twist, propellertype);% [deg]        
    [Cx, Cy, SectionResults(i,:)] = SolveSection(i, polar, prop, air, oper, AzimuthAngle, yaw, dPsi);
end
figure;
plot(new(:,1),new(:,2),'k');
% hold on
% plot(new(:,1), a_zeroyaw,'b');

%% Functions
function [ftip, froot, ftotal] = PrandtlTipRootCorrection(prop, SectionRadius, oper, a)
    r_R_Section = SectionRadius/prop.R; %non-dimensional
    F1 = (-prop.Nblades/2)*((1-r_R_Section)/r_R_Section)*sqrt(1+(oper.TSR*r_R_Section/(1-a))^2);
    F2 = (-prop.Nblades/2)*((r_R_Section-prop.blade_root)/r_R_Section)*sqrt(1+(oper.TSR*r_R_Section/(1-a))^2);
    
%     F1 = (-prop.Nblades/2).*((1-r_R_Section)./r_R_Section).*(1./sin(InflowAngle));
%     F2 = (-prop.Nblades/2).*((r_R_Section-prop.blade_root)./r_R_Section).*(1./sin(InflowAngle));
    
    ftip = (2/pi).*acos(exp(F1));    
    if (isnan(ftip))
       ftip = 0; 
    end
    
    froot = (2/pi).*acos(exp(F2));
    if (isnan(froot))
       froot = 0;
    end
    
    ftotal = ftip.*froot;
    if (ftotal < 0.0001)
        ftotal = 0.0001;
    end
end
function [InflowAngles, alphas, cl, cd] = AnnularRingProperties(yaw, a, aprime, prop, oper, SectionRadius, AzimuthAngle, polar, index)
    [InflowAngles] = inflowangle(yaw, a, aprime, prop, oper, SectionRadius, AzimuthAngle);%should be an array of phi for constant radius
    alphas = rad2deg(InflowAngles) - prop.sectionpitch;
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
function Uresultant2 = ResultantVelocityBladeElement(SectionRadius, oper, AzimuthAngle, yaw, a, aprime)
    chi = WakeSkewAngle(yaw, a);
    Vx = oper.U_inf;
    Vy = oper.omega*SectionRadius;
    %1st term
    temp1_1 = Vx*(cos(yaw)-a);
    temp1_2 = Vy.*aprime.*sin(chi).*cos(AzimuthAngle).*(1+sin(AzimuthAngle).*sin(chi));
    temp1 =  temp1_1 + temp1_2;
    
    %2nd term
    temp2_1 =  Vy.*(1+aprime.*cos(chi).*(1+sin(AzimuthAngle).*sin(chi)));
    temp2_2 = Vx*cos(AzimuthAngle).*(a*tan(chi/2)-sin(yaw));
    temp2 =  temp2_1 + temp2_2;
    Uresultant2 = (temp1.^2)+(temp2.^2);
end
function Func = CoeffMomentumFuncGlauert(a_new, F, yaw, CtRHS, CqRHS, oper, SectionRadius, chi)    
    %a(1) = a, a(2) = aprime   
    %momentum theory from Glauert

   
    tmp= 2*cos(yaw) - a_new(1)*F;
    CtFunc= CtRHS-(4*a_new(1)*F*sqrt(1-a_new(1)*F*tmp));
       
    
    Vx = oper.U_inf;
    Vy = oper.omega*SectionRadius;
    CqFunc =  CqRHS-(2*(Vy/Vx)*a_new(2)*F*(cos(yaw)-a_new(1)*F)*(1+(cos(chi))^2));  
    
    Func = [CtFunc;CqFunc];%vector of 2 functions to be solved
    
end
function Func = CoeffMomentumFuncPropBrake(a_new, F, yaw, CtRHS, CqRHS, oper, SectionRadius, chi)    
    %a(1) = a, a(2) = aprime   
    %use different Ct_momentum functions for different a   
    
    %propeller-brake
    CtFunc= CtRHS-(4*a_new(1)*F*(a_new(1)*F-1)*cos(yaw));
        
    Vx = oper.U_inf;
    Vy = oper.omega*SectionRadius;
    CqFunc =  CqRHS-(2*(Vy/Vx)*a_new(2)*F*(cos(yaw)-a_new(1)*F)*(1+(cos(chi))^2));  
    
    Func = [CtFunc;CqFunc];%vector of 2 functions to be solved
    
end
function Func = CoeffMomentumFuncEmpirical(a_new, F, yaw, CtRHS, CqRHS, oper, SectionRadius, chi)    
    %a(1) = a, a(2) = aprime   
    %use different Ct_momentum functions for different a   
    beta = 0.4;

  
    [f0, fprime_0, f1]=fcoeff(F, yaw);
    k2 = double((f1-f0-fprime_0*(1-beta))/((1-beta)^2));
    k1 = double(fprime_0-2*k2*beta);
    k0 = double(f1-k1-k2);
    
    
    CtFunc= CtRHS-(k0+(k1*a_new(1))+(k2*(a_new(1)^2)));
    
    
    Vx = oper.U_inf;
    Vy = oper.omega*SectionRadius;
    CqFunc =  CqRHS-(2*(Vy/Vx)*a_new(2)*F*(cos(yaw)-a_new(1)*F)*(1+(cos(chi))^2));  
    
    Func = [CtFunc;CqFunc];%vector of 2 functions to be solved
    
end
function  [f0, fprime_0, f1] = fcoeff(F, yaw)
    beta = 0.4;
    syms x
    g = 4*x*F*sqrt(1-x*F*(2*cos(yaw) - x*F));
    y = diff(g);
    f0 = double(vpa(subs(g,x,beta)));
    fprime_0=double(vpa(subs(y,x,beta)));
    f1 = double(2*cos(yaw));
%     test = [f0, fprime_0, f1];
end
function a_new = solveCoeffFunc(a_previous, F, yaw, CtRHS, CqRHS, oper, SectionRadius, chi)
    options = optimoptions('fsolve','Display','off');
    initial = [0.3;0.1];
    beta = 0.4;
    if a_previous<=beta
%         disp('Glauert');
        a_new = fsolve(@(a_new) CoeffMomentumFuncGlauert(a_new, F, yaw, CtRHS, CqRHS, oper, SectionRadius, chi), initial, options);
    elseif a_previous>1
%         disp('propBrake');
        a_new = fsolve(@(a_new) CoeffMomentumFuncPropBrake(a_new, F, yaw, CtRHS, CqRHS, oper, SectionRadius, chi), initial, options);
    else
%         disp('Empirical');
        a_new = fsolve(@(a_new) CoeffMomentumFuncEmpirical(a_new, F, yaw, CtRHS, CqRHS, oper, SectionRadius, chi), initial, options);
    end

end
function [Results] = SolveSection(index, polar, prop, ~, oper, AzimuthAngle, yaw, dPsi)
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
    
    for i=1:N
        %produce array of local AoA, Inflow angles, cl and cd
        [InflowAngles, alphas, cl, cd] = AnnularRingProperties(yaw, a, aprime, prop, oper, SectionRadius, AzimuthAngle, polar, index);        
        chi = WakeSkewAngle(yaw, a);
        
        W2 = ResultantVelocityBladeElement(SectionRadius, oper, AzimuthAngle, yaw, a, aprime);%should be array of W here
        Cx = cl.*cos(InflowAngles)+cd.*sin(InflowAngles);
        Cy = cl.*sin(InflowAngles)-cd.*cos(InflowAngles);
        [~, ~, Prandtl]=PrandtlTipRootCorrection(prop, SectionRadius, oper, a);
        
        %calculate integrals on RHS here

        CtIntegral = (W2./(oper.U_inf^2)).*Cx;        
        CtElementAzimAverage = (SectionRotorSolidity/(2*pi))*trapz(AzimuthAngle, CtIntegral);
        temp = Cy.*cos(chi)-Cx.*sin(chi).*cos(AzimuthAngle);
        CqIntegral = (W2./(oper.U_inf^2)).*(temp);
        CqElementAzimAverage = (SectionRotorSolidity/(2*pi))*trapz(AzimuthAngle, CqIntegral);
        %solve for new a and aprime using fsolve
        NewValues = solveCoeffFunc(a, Prandtl, yaw, CtElementAzimAverage, CqElementAzimAverage, oper, SectionRadius, chi);
        
        a_new = NewValues(1);
        aprime_new = NewValues(2);
                 
        if (abs(a-a_new)<epsilon)
            break;
        else
            a = 0.75*a+0.25*a_new;            
        end
       
    end 
    
    Results = [SectionRadius/prop.R, a_new, aprime_new, i];
end
function F = FlowExpansionFunction(r_R)
    %curve fit to original function
    F = 0.5.*(r_R+0.4.*(r_R.^3)+0.4.*(r_R.^5));

end
function chi = WakeSkewAngle(yaw, a)
    %yaw: [rad]
    chi = yaw*(0.6*a+1);
end
function K_c = K(yaw, a)
    chi = WakeSkewAngle(yaw, a);
    K_c = 2*tan(chi/2);
end
%% Function handles for integration in terms of Azimuth Angle (not in use)
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
