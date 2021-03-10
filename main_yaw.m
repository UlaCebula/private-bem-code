clear all
close all
clc

DEBUG = 0;
spacing = 'constant';
type = 'WT';
[polar, prop, oper, air, propellertype] = param(spacing, type);

%%%
yaw = deg2rad(15); 
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
    prop.sectionchord = chord_distribution(0.5*(prop.r_R(i)+prop.r_R(i+1)), prop.R, propellertype);%non-dimensional
    prop.sectionpitch = pitch_distribution(0.5*(prop.r_R(i)+prop.r_R(i+1)),prop.collective_blade_twist, propellertype);% [deg]        
    [SectionResults(i,:), Inflowangles(i,:), alphas(i,:), temp1] = SolveSection(i, polar, prop, air, oper, AzimuthAngle, yaw, dPsi);
end


%% Functions
function [ftip, froot, ftotal] = PrandtlTipRootCorrection(prop, SectionRadius, oper, a)
    r_R_Section = SectionRadius/prop.R; %non-dimensional
    F1 = (-prop.Nblades/2)*((1-r_R_Section)/r_R_Section)*sqrt(1+(oper.TSR*r_R_Section/(1-a))^2);
    F2 = (-prop.Nblades/2)*((r_R_Section-prop.blade_root)/r_R_Section)*sqrt(1+(oper.TSR*r_R_Section/(1-a))^2);
    
%     F1 = (-prop.Nblades/2)*((1-r_R_Section)/r_R_Section)*(1/sin(InflowAngle));
%     F2 = (-prop.Nblades/2)*((r_R_Section-prop.blade_root)/r_R_Section)*(1/sin(InflowAngle));
    
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
function [temp1, InflowAngles, alphas, cl, cd] = AnnularRingProperties(yaw, a, aprime, prop, oper, SectionRadius, AzimuthAngle, polar, index)
    [temp1, InflowAngles] = inflowangle(yaw, a, aprime, prop, oper, SectionRadius, AzimuthAngle);%should be an array of phi for constant radius
    alphas = rad2deg(InflowAngles) - prop.sectionpitch;
    for i=1:length(alphas)
        if alphas(i)>max(polar.alpha)
            alphas(i) = max(polar.alpha);
        end
        if alphas(i)<min(polar.alpha)
            alphas(i) = min(polar.alpha);
        end
    end
%     alphas
    cl = interp1(polar.alpha, polar.Cl, alphas);
    cd = interp1(polar.alpha, polar.Cd, alphas);
end
function [temp1, phi] = inflowangle(yaw, a, aprime, prop, oper, SectionRadius, AzimuthAngle)
    %Equation 3.131
    K_c = K(yaw, a);%constant for constant yaw, 0 for 0 deg yaw
    chi = WakeSkewAngle(yaw, a);%constant for constant yaw, 0 for 0 deg yaw
    r_R = SectionRadius./prop.R;
    F = FlowExpansionFunction(r_R);%constant for constant radius
    temp1 = a.*(1+F.*K_c.*sin(AzimuthAngle));
    top = oper.U_inf.*(cos(yaw)-temp1)+(oper.omega.*SectionRadius.*aprime.*cos(AzimuthAngle).*sin(chi).*(1+sin(AzimuthAngle).*sin(chi)));
    temp1 = oper.omega.*SectionRadius.*(1+aprime.*cos(chi).*(1+sin(AzimuthAngle).*sin(chi)));
    temp2 = a.*tan(chi/2).*(1+F.*K_c.*sin(AzimuthAngle)) - sin(yaw);
    btm = temp1 + oper.U_inf.*cos(AzimuthAngle).*temp2;
    phi = atan(top./btm);
end
function Uresultant = ResultantVelocityBladeElement(SectionRadius, prop, oper, AzimuthAngle, yaw, a, aprime)
    % Equation 3.139 in wind energy book
    r_R = SectionRadius/prop.R;
    chi = WakeSkewAngle(yaw, a);
    temp1 = (cos(yaw)-a)+oper.TSR.*r_R.*aprime.*cos(AzimuthAngle).*sin(chi).*(1+sin(AzimuthAngle).*sin(chi));
    temp2 = (oper.TSR.*r_R.*(1+aprime.*cos(chi).*(1+sin(AzimuthAngle).*sin(chi)))+cos(AzimuthAngle).*(atan(chi/2)-sin(yaw)));
    Uresultant = sqrt((oper.U_inf.^2).*((temp1.^2)+(temp2.^2)));
end
function [a] = CalculateNewAxialInduction(W, SectionRotorSolidity, Cx, dPsi, Prandtl, yaw, chi, oper)
    %solve integral discretely
    NormalisedW = (W.^2)/(oper.U_inf^2); 
    temp = SectionRotorSolidity.*sum(Cx.*NormalisedW.*dPsi);
%     Prandtl;
%     yaw;
%     chi;
%     cos(yaw)+tan(chi/2)*sin(yaw)-0*Prandtl*(sec(chi/2))^2;
    %solve for a
%     fun = @(a) (8*pi*a*Prandtl*sqrt(1-a*Prandtl*(2*cos(yaw)-(a*Prandtl))))-temp;
    fun = @(a) 	8*pi*a*Prandtl*sqrt(1-a*Prandtl*((2*cos(yaw))-(a*Prandtl)))-temp;
    a = fzero(fun, 0.3)    
end
function aprime = CalculateNewAzimInduction(W, SectionRotorSolidity, Cx, Cy, dPsi, AzimuthAngle, chi, yaw, Prandtl, a, oper, SectionRadius, prop)    
    NormalisedW = (W.^2)/(oper.U_inf^2); 
    r_R = SectionRadius/prop.R;
    temp1 = sum(NormalisedW.*(cos(AzimuthAngle).*sin(chi).*Cx)+(cos(chi).*Cy).*dPsi);
    temp = SectionRotorSolidity.*temp1;
    
    fun = @(aprime) (4*aprime*Prandtl*(cos(yaw)-(a*Prandtl))*oper.TSR*r_R*pi*(1+(cos(chi)).^2))-temp;
    aprime = fzero(fun, 0.5);  

end
function [Results, InflowAngles, alphas, temp1] = SolveSection(index, polar, prop, ~, oper, AzimuthAngle, yaw, dPsi)
    dPsi = deg2rad(dPsi);
    r_R1=prop.r_R(index); %non-dimensional
    r_R2=prop.r_R(index+1);
    SectionArea = pi*(((r_R2*prop.R)^2)-((r_R1*prop.R)^2)); %[m^2]
    SectionRadius = 0.5*(r_R1+r_R2)*prop.R;%average radius [m]
    %initialising
    a = 0.1;%axial induction factor
    aprime = 0.1;%tangential induction factor
    
    N = 100;%number of iterations for convergence
    epsilon = 0.00001;
    
    for i=1:N
        %produce array of local AoA, Inflow angles, cl and cd
        [temp1, InflowAngles, alphas, cl, cd] = AnnularRingProperties(yaw, a, aprime, prop, oper, SectionRadius, AzimuthAngle, polar, index);
        
        chi = WakeSkewAngle(yaw, a);
        W = ResultantVelocityBladeElement(SectionRadius, prop, oper, AzimuthAngle, yaw, a, aprime);%should be array of W here
        Cx = cl.*cos(InflowAngles)+cd.*sin(InflowAngles);
        Cy = cl.*sin(InflowAngles)-cd.*cos(InflowAngles);
        SectionRotorSolidity = (prop.Nblades*prop.sectionchord*prop.R)/(2*pi*SectionRadius);
        [~, ~, Prandtl]=PrandtlTipRootCorrection(prop, SectionRadius, oper, a);
        
        %%calculate new value of a here
        a_new = CalculateNewAxialInduction(W, SectionRotorSolidity, Cx, dPsi, Prandtl, yaw, chi, oper);

        %%calculate new value of aprime here
        [~, ~, Prandtl]=PrandtlTipRootCorrection(prop, SectionRadius, oper, a_new);
        aprime= CalculateNewAzimInduction(W, SectionRotorSolidity, Cx, Cy, dPsi, AzimuthAngle, chi, yaw, Prandtl, a_new, oper, SectionRadius, prop);
        
%         a = 0.75*a+0.25*a_new;%update new value of axial induction via weighted average
        

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Thrust coefficient at streamtube        U_axial, U_tangential
%         load3D_axial = f_axial*prop.Nblades*(prop.dr(index)*prop.R);
% %         CT_streamtube= load3D_axial./(0.5*air.density*(oper.U_inf^2)*SectionArea);
%         CT_streamtube= load3D_axial./(0.5*(oper.U_inf^2)*SectionArea);
%         a_new=AxialInductionFactor(CT_streamtube);
%         
% %         [~, ~, Prandtl]=PrandtlTipRootCorrection(prop, SectionRadius, oper, a_new);
%         a_new = a_new/Prandtl;
%         a = 0.75*a+0.25*a_new;%update new value of axial induction via weighted average
%         
%         aprime = AzimInductionFactor(f_tangential, prop, SectionRadius, oper, a, air);
%         aprime = aprime/Prandtl;
        
        if (abs(a-a_new)<epsilon)
            break;
        end
        a = 0.75*a+0.25*a_new;
        
        if a>0.95
            a = 0.95;
        end
    end 
    
    Results = [SectionRadius/prop.R, a_new, aprime,  i];
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