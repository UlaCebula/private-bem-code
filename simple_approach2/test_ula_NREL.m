%% import polars
clc
close all
clear all

opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [2, Inf];
opts.Delimiter = " ";
opts.VariableNames = ["Alpha", "Cl", "Cd", "Cm"];
opts.VariableTypes = ["double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
DU95W180 = readtable("C:\Users\ula\Desktop\Delft\Q3\Rotor-Wake\private-bem-code-main\DU95W180.txt", opts);
clear opts
polar_alpha = DU95W180.Alpha;
polar_cl = DU95W180.Cl;
polar_cd = DU95W180.Cd;

%plot polars of the airfoil C-alfa and Cl-Cd
% figure()
% subplot(1,2,1);
% plot(polar_alpha,polar_cl);
% xlabel('\alpha');
% ylabel('C_{l}');
% xlim([-30 30]);
% grid on
% grid minor
% 
% subplot(1,2,2);
% plot(polar_alpha,polar_cd);
% xlabel('\alpha');
% ylabel('C_{d}');
% % xlim([0 0.1]);
% grid on
% grid minor

%%
% define the blade geometry
delta_r_R = 0.01;
r_R = [0.2:delta_r_R:1+delta_r_R/2];

% blade shape
pitch = 2; % degrees
chord_distribution = 3*(1-r_R)+1; % meters
twist_distribution = -14*(1-r_R)+pitch; % degrees

% define flow conditions
Uinf = 1; % unperturbed wind speed in m/s
TSR = 8; % tip speed ratio
Radius = 50;
Omega = Uinf*TSR/Radius;
NBlades = 3;
yaw=deg2rad(0);

TipLocation_R =  1;
RootLocation_R =  0.2;

% solve BEM model
results =zeros(length(r_R)-1,6);

for i=1:(length((r_R))-1)
    chord = interp1(r_R, chord_distribution, (r_R(i)+r_R(i+1))/2);
    twist = interp1(r_R, twist_distribution, (r_R(i)+r_R(i+1))/2);
    results(i,:) = solveStreamtube(Uinf, r_R(i), r_R(i+1), RootLocation_R, TipLocation_R , Omega, Radius, NBlades, chord, twist, polar_alpha, polar_cl, polar_cd, yaw);
end

%%
% plot results
areas = (r_R(2:end).^2-(r_R(1:end-1).^2))*pi*Radius^2;
dr = (r_R(2:end)-r_R(1:end-1))*Radius;
CT = sum(dr*results(:,4)*NBlades/(0.5*Uinf^2*pi*Radius^2));
CP = sum(dr*results(:,5)*results(:,3)*NBlades*Radius*Omega/(0.5*Uinf^3*pi*Radius^2));

fprintf("CT is %d\n", CT)
fprintf("CP is %d\n", CP)

figure(2)
title('Axial and tangential induction')
hold on
plot(results(:,3), results(:,1), 'm-')
plot(results(:,3), results(:,2), 'm--')
grid on
grid minor
xlabel('r/R')

% figure()
% title('Normal and tagential force, non-dimensioned')
% hold on
% plot(results(:,3), results(:,4)/(0.5*Uinf^2*Radius), 'r-')
% plot(results(:,3), results(:,5)/(0.5*Uinf^2*Radius), 'g--')
% grid on
% grid minor
% xlabel('r/R')

% figure()
% hold on
% title('Circulation distribution, non-dimensioned')
% plot(results(:,3), results(:,6)/(pi*Uinf^2/(NBlades*Omega)), 'r-')
% grid on
% grid minor
% xlabel('r/R')


%% Functions
function CT=CTfunction(a, glauert)
% This function calculates the thrust coefficient as a function of induction factor 'a'
% 'glauert' defines if the Glauert correction for heavily loaded rotors should be used; default value is false
CT = zeros(length(a));
CT = 4*a.*(1-a);
if glauert
    CT1=1.816;
    a1=1-sqrt(CT1)/2;
    CT(a>a1) = CT1-4*(sqrt(CT1)-1)*(1-a(a>a1));
end
end

function a = ainduction(CT)
% This function calculates the induction factor 'a' as a function of thrust coefficient CT
% including Glauert's correction
a = zeros(1,length(CT));
CT1=1.816;
CT2=2*sqrt(CT1)-CT1;
a(CT>=CT2) = 1 + (CT(CT>=CT2)-CT1)/(4*(sqrt(CT1)-1));
a(CT<CT2) = 0.5-0.5*sqrt(1-CT(CT<CT2));
end

function F_total=PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, TSR, NBlades, axial_induction)
%     This function calcualte steh combined tip and root Prandtl correction at agiven radial position 'r_R' (non-dimensioned by rotor radius),
%     given a root and tip radius (also non-dimensioned), a tip speed ratio TSR, the number lf blades NBlades and the axial induction factor
temp1 = -NBlades/2*(tipradius_R-r_R)./r_R*sqrt(1+((TSR*r_R).^2)./((1-axial_induction).^2));
Ftip = 2/pi*acos(exp(temp1));
Ftip(isnan(Ftip)) = 0;
temp1 = NBlades/2*(rootradius_R-r_R)./r_R*sqrt( 1+ ((TSR*r_R).^2)./((1-axial_induction).^2));
Froot = 2/pi*acos(exp(temp1));
Froot(isnan(Froot)) = 0;
F_total = Froot.*Ftip;
end

% define function to determine load in the blade element
function [inflowangle, fnorm , ftan, gamma, cx, cy] = loadBladeElement(r_R, chord, twist, polar_alpha, polar_cl, polar_cd, Uaxial, Utan_yaw)
%     calculates the load in the blade element
vmag2 = Uaxial.^2 + Utan_yaw.^2;
inflowangle = atan2(Uaxial, Utan_yaw);
alpha = twist + inflowangle*180/pi;
cl = interp1(polar_alpha, polar_cl, alpha);
cd = interp1(polar_alpha, polar_cd, alpha);
cx = cl.*cos(inflowangle)+cd.*sin(inflowangle);
cy = cl.*sin(inflowangle)-cd.*cos(inflowangle);
lift = 0.5*vmag2.*cl*chord;
drag = 0.5*vmag2.*cd*chord;
fnorm = lift.*cos(inflowangle)+drag.*sin(inflowangle);
ftan = lift.*sin(inflowangle)-drag.*cos(inflowangle);
gamma = 0.5*sqrt(vmag2).*cl*chord;
end

function result = solveStreamtube(Uinf, r1_R, r2_R, rootradius_R, tipradius_R , Omega, Radius, NBlades, chord, twist, polar_alpha, polar_cl, polar_cd, yaw)
Area = pi*((r2_R*Radius)^2-(r1_R*Radius)^2); % area streamtube
r_R = (r1_R+r2_R)/2; % centroide
a = 0.01; % axial induction
aline = 0.01; % tangential induction factor
psi = deg2rad([0:1:359]);

Niterations = 1000;
Erroriterations =0.00001; % error limit for iteration rpocess, in absolute value of induction

for i=1:Niterations
    Uaxial = Uinf*cos(yaw)*(1-a);   % axial velocity at rotor
    Utan_yaw = (1+aline).*(Omega*r_R*Radius-Uinf*sin(yaw)*cos(psi)); % tangential velocity at rotor
    [inflowangle, fnorm, ftan, gamma, cx, cy] = loadBladeElement(r_R,chord, twist, polar_alpha, polar_cl, polar_cd, Uaxial, Utan_yaw);
    
    CTsegment = 4*a.*sqrt(1-a.*(2*cos(psi)-a));
    sigma = NBlades*chord ./ (2*pi*r_R*Radius);
    F = PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, Omega*Radius/Uinf, NBlades, a);
    kappa = sigma*cx./(4*F.*sin(inflowangle).^2);
    kappap = sigma*cy./(4*F.*sin(inflowangle).*cos(inflowangle));
    
    % if kappa > 2/3
    gam1 = 2*F.*kappa-(10/9-F);
    gam2 = 2*F.*kappa-F.*(4/3-F);
    gam3 = 2*F.*kappa-(25/9-2*F);
    a_new = (gam1-sqrt(gam2))./gam3;
%     if kappa<=2/3
%         a_new = kappa/(1+kappa);
%     if kappa<1
%         a_new = kappa/(kappa-1);
%     end
    
    chi = yaw*(1+0.6*a_new);
    a_new = a_new.*(1+15*pi/64.*r_R*tan(chi/2).*sin(psi));
    aline_new = kappap./(1-kappap);
    a = 0.75*a+0.25*a_new; % for improving convergence, weigh current and previous iteration of axial induction
    aline = 0.75*aline+0.25*aline_new;

    if all((abs(a-a_new)) < Erroriterations)
        % print("iterations")
        % print(i)
        result = [mean(a), mean(aline), r_R, mean(fnorm) , mean(ftan), mean(gamma)];
        return
    end
    if i==Niterations
       disp("Not converged") 
    end
end
end