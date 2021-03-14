clear all
close all
clc

%Based on NREL methodology

% LOAD FROM PARAM.M
DEBUG = 0;
spacing = 'constant';
type = 'WT';
[polar, prop, oper, air, propellertype] = param(spacing, type);

%% carlos' data
zeroyawresults = load('zeroyawdata.mat');
radialdist_zeroyaw = zeroyawresults.SectionResults(:,1);
a_zeroyaw = zeroyawresults.SectionResults(:,2);
aprime_zeroyaw = zeroyawresults.SectionResults(:,3);
%% yaw/Azimuth angles
yaw = deg2rad([0;15;30]); 
dPsi = 1;%[deg]
AzimuthAngle = deg2rad([0:dPsi:360]');%starts from 12 o'clock position
%% polars

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
for j=1:length(yaw)
    for i=1:length(prop.SectionRadius)
            % TODO: extract all parameters required to calculate
                % - Spanwise distribution of angle of attack and inflow angle
                % - Spanwise distribution of axial and azimuthal inductions
                % - Spanwise distribution of thrust and azimuthal loading
                % - Total thrust and torque
                % - For the cases of yawed HAWT, also plot the azimuthal variation (suggestion: polar contour plot)
            [a(i,:), aprime(i,:),] = SolveSection(i, polar, prop,  oper, AzimuthAngle, yaw(j));
            %a and aprime should be a matrix of size(number of sections, number of AzimuthAngles)
    end
end
% CALCULATION OF VARIOUS COEFFICIENTS AND REQUIRED VALUES GO BELOW HERE
%% Functions
function [ftip, froot, ftotal] = PrandtlTipRootCorrection(prop, InflowAngle)
    r_R_Section = prop.SectionRadius./prop.R; %non-dimensional
    
    %(prop, SectionRadius, oper, a)
%     F1 = (-prop.Nblades/2)*((1-r_R_Section)/r_R_Section)*sqrt(1+(oper.TSR*r_R_Section/(1-a))^2);
%     F2 = (-prop.Nblades/2)*((r_R_Section-prop.blade_root)/r_R_Section)*sqrt(1+(oper.TSR*r_R_Section/(1-a))^2);
    
    %(prop, SectionRadius, InflowAngle)
    F1 = (-prop.Nblades./2).*((1-r_R_Section)./r_R_Section).*(1./sin(InflowAngle));
    F2 = (-prop.Nblades./2).*((r_R_Section-prop.blade_root)./prop.blade_root).*(1./sin(InflowAngle));
    
    ftip = (2/pi).*acos(exp(F1));    
    ftip(isnan(ftip)) = 0;
    
    froot = (2/pi).*acos(exp(F2));
    froot(isnan(froot)) = 0;
    
    ftotal = ftip.*froot;
    ftotal(ftotal<0.0001) = 0.0001;
end
function [a_new, aprime_new] = SolveSection(index, polar, prop, oper, AzimuthAngle, yaw)
    %index: to keep track of radial position of elements
    
    %NOTE: prop.SectionRadius is dimensional [m]
    
    %initialising
    %induction factors should be vector of length(AzimuthAngle)
    a = 0.3.*ones(length(AzimuthAngle),1);%axial induction factor
    aprime = 0.1.*ones(length(AzimuthAngle),1);%tangential induction factor
    
    N = 100;%number of iterations for convergence
    epsilon = 0.0001;%decrease tolerance to speed up convergence if need be
    
    % calculate Ux and Uy, taking yaw and azimuth angle into account
    % TODO: check equation for Ux and Uy
    % Uy will be vector of length(AzimuthAngle)
    Ux = oper.U_inf.*(cos(yaw));
    Uy = oper.U_inf.*(-sin(yaw).*cos(AzimuthAngle))+(oper.omega.*prop.SectionRadius);
    
    %calculate a and aprime for each section at a constant radius for each azimuth angle 
    for i=1:N
        %produce array of local AoA, Inflow angles, cl and cd
        [inflowangle, alpha, cx, cy] = SectionInflowAngleCalc(oper, a, aprime, SectionRadius, polar, prop, yaw, AzimuthAngle);
        [~, ~, Prandtl] = PrandtlTipRootCorrection(prop, inflowangle);

        %TODO: calculate correct value of a_new here
        a_new = AxialInductionFactor(Prandtl, CTsection, SectionRotorSolidity, cx, inflowangle);

        % TODO: Correct effect of skew for a_new here, check if functions
        % are correct
        a_new = SkewedWakeCorrection(a_new,prop,AzimuthAngle, yaw);
        %TODO: calculate correct value of aprime_new here
        aprime_new = AzimInductionFactor(Prandtl, prop, cy, inflowangle);       
        
        % This break condition should work...
        if (max(abs(a-a_new))<epsilon)
            break;
        else
            a = 0.75.*a+0.25.*a_new;%,update new value of axial induction via weighted average
        end
    end
    
    
    Results = [SectionRadius/prop.R, i];% i represents max iterations done per fixed radius
end
function [phi, alpha, cx, cy] = SectionInflowAngleCalc(oper, a, aprime, SectionRadius, polar, prop, yaw, AzimuthAngle)
    U_axial = Ux.*(1-a);%
    U_tan = Uy.*(1+aprime);%
    %TODO: calculate U_resultant as well  
    
    
    phi = atan2(U_axial,U_tan);
    alpha = rad2deg(phi) - prop.sectionpitch;%vector    
    cl = interp1(polar.alpha, polar.Cl, alpha,'linear', 'extrap');
    cd = interp1(polar.alpha, polar.Cd, alpha,'linear', 'extrap');
    cx = cl.*cos(phi) + cd.*sin(phi);
    cy = cl.*sin(phi) - cd.*cos(phi);
end

function a = AxialInductionFactor(F, CT, sigma, cx, phi)
    % TODO: return correct value of a
end
function aprime = AzimInductionFactor(F, prop, cy, phi)
    % TODO: return correct value of a
end
function wakeskewangle = WakeSkewAngleCalc(yaw, a)
    % TODO: check if this is correct
    %yaw: [rad]
    wakeskewangle = yaw.*(0.6.*a+1);%scalar/vector
end
function a_skew = SkewedWakeCorrection(a,prop,AzimuthAngle,yaw)
    % TODO: check if this is correct
    chi = WakeSkewAngleCalc(yaw, a);
    r_R = prop.SectionRadius./prop.R;
    temp = (15*pi/64).*r_R.*tan(chi/2).*sin(AzimuthAngle);
    a_skew = a.*(1+temp);%should be vector
end
