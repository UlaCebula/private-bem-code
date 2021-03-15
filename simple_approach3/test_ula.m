% plot CT as a function of induction "a", with and without Glauert correction
% define a as a range
% a = [-0.5:0.1:1];
% CTmom = CTfunction(a,0); % CT without correction
% CTglauert = CTfunction(a, 1); % CT with Glauert's correction
% a2 = ainduction(CTglauert);
%
% figure()
% hold on
% plot(a, CTmom, 'k-')
% plot(a, CTglauert, 'b--')
% plot(a, CTglauert.*(1-a), 'g--')
% xlabel('a')
% grid on
% grid minor
% ylabel(r'$C_T$ and $C_P$')

%%
% plot Prandtl tip, root and combined correction for a number of blades and induction 'a', over the non-dimensioned radius
% r_R = [0.1:0.01:1];
% a = zeros(size(r_R))+0.3;
% [Prandtl, Prandtltip, Prandtlroot] = PrandtlTipRootCorrection(r_R, 0.1, 1, 7, 3, a)
%
% figure()
% hold on
% plot(r_R, Prandtl, 'r-')
% plot(r_R, Prandtltip, 'g.')
% plot(r_R, Prandtlroot, 'b.')
% xlabel('r/R')
% grid on
% grid minor

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
plot(results(:,3), results(:,1), 'g-')
plot(results(:,3), results(:,2), 'g--')
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

function a = ainduction(CT, a_old, yaw)
% This function calculates the induction factor 'a' as a function of thrust coefficient CT
% including Glauert's correction
a = zeros(1,length(CT));
CT1=1.816;
CT2=2*sqrt(CT1)-CT1;
a(CT>=CT2) = 1 + (CT(CT>=CT2)-CT1)/(4*(sqrt(CT1)-1));
% a(CT<CT2) = 0.5-0.5*sqrt(1-CT(CT<CT2));
% a(CT<CT2) = CT(CT<CT2)./(4*sqrt(1-a_old(CT<CT2).*(2*cos(yaw) - a_old(CT<CT2))));

for i=1:length(a)
    if a(i)<CT2
        a_func = @(x) 4*x*sqrt(1-x*(2*cos(yaw) - x))-CT(i);
        a(i) = fzero(a_func, a_old(i));
    end
end
% a_func = @(x) 4*x.*sqrt(1-x.*(2*cos(yaw) - x))-CT;
% a(CT<CT2) = arrayfun(@(x0) fsolve(@a_func,x0), a(CT<CT2));

end

function [F_total, Ftip, Froot]=PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, TSR, NBlades, axial_induction)
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
function [fnorm , ftan, gamma] = loadBladeElement(r_R, chord, twist, polar_alpha, polar_cl, polar_cd, Uaxial, Utan_yaw)
%     calculates the load in the blade element
vmag2 = Uaxial.^2 + Utan_yaw.^2;
inflowangle = atan2(Uaxial, Utan_yaw);
alpha = twist + inflowangle*180/pi;
cl = interp1(polar_alpha, polar_cl, alpha);
cd = interp1(polar_alpha, polar_cd, alpha);
lift = 0.5*vmag2.*cl*chord;
drag = 0.5*vmag2.*cd*chord;
fnorm = lift.*cos(inflowangle)+drag.*sin(inflowangle);
ftan = lift.*sin(inflowangle)-drag.*cos(inflowangle);
gamma = 0.5*sqrt(vmag2).*cl*chord;
end

function result = solveStreamtube(Uinf, r1_R, r2_R, rootradius_R, tipradius_R , Omega, Radius, NBlades, chord, twist, polar_alpha, polar_cl, polar_cd, yaw)
%     solve balance of momentum between blade element load and loading in the streamtube
%     input variables:
%     Uinf - wind speed at infinity
%     r1_R,r2_R - edges of blade element, in fraction of Radius ;
%     rootradius_R, tipradius_R - location of blade root and tip, in fraction of Radius ;
%     Radius is the rotor radius
%     Omega -rotational velocity
%     NBlades - number of blades in rotor
Area = pi*((r2_R*Radius)^2-(r1_R*Radius)^2); % area streamtube
r_R = (r1_R+r2_R)/2; % centroide    
psi = deg2rad([0:1:359]);

% initiatlize variables
a = zeros(1,length(psi)); % axial induction
aline = zeros(1,length(psi)); % tangential induction factor

Niterations = 10000;
Erroriterations =0.00001; % error limit for iteration rpocess, in absolute value of induction

for i=1:Niterations
    %      ///////////////////////////////////////////////////////////////////////
    %      // this is the block "Calculate velocity and loads at blade element"
    %      ///////////////////////////////////////////////////////////////////////
    Uaxial = Uinf*cos(yaw)*(1-a);
    Utan_yaw = r_R*Radius*Omega*(1+aline)-Uinf*sin(yaw)*cos(psi);
    
    % calculate loads in blade segment in 2D (N/m)
    [fnorm, ftan, gamma] = loadBladeElement(r_R,chord, twist, polar_alpha, polar_cl, polar_cd, Uaxial, Utan_yaw);
    load3Daxial =fnorm*Radius*(r2_R-r1_R)*NBlades; % 3D force in axial direction
    
    % ///////////////////////////////////////////////////////////////////////
    % //the block "Calculate velocity and loads at blade element" is done
    % ///////////////////////////////////////////////////////////////////////
    
    % ///////////////////////////////////////////////////////////////////////
    % // this is the block "Calculate new estimate of axial and azimuthal induction"
    % ///////////////////////////////////////////////////////////////////////
    % // calculate thrust coefficient at the streamtube
    CT = load3Daxial/(0.5*Area*Uinf^2);
    
    % calculate new axial induction, accounting for Glauert's correction
    anew =  ainduction(CT, a, yaw);
    
    % correct new axial induction with Prandtl's correction
    [Prandtl, Prandtltip, Prandtlroot] = PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, Omega*Radius/Uinf, NBlades, anew);
    if (Prandtl < 0.0001)
        Prandtl = 0.0001; % avoid divide by zero
    end
    anew = anew./Prandtl; % correct estimate of axial induction
    a = 0.75*a+0.25*anew; % for improving convergence, weigh current and previous iteration of axial induction
    
    % calculate aximuthal induction
    aline = ftan*NBlades/(2*pi*Uinf*(1-a)*Omega*2*(r_R*Radius)^2);
    aline = aline./Prandtl; % correct estimate of azimuthal induction with Prandtl's correction
    % ///////////////////////////////////////////////////////////////////////////
    % // end of the block "Calculate new estimate of axial and azimuthal induction"
    % ///////////////////////////////////////////////////////////////////////
    
    %// test convergence of solution, by checking convergence of axial induction
    if (abs(a-anew) < Erroriterations)
        % print("iterations")
        % print(i)
        result = [mean(a), mean(aline), r_R, mean(fnorm) , mean(ftan), mean(gamma)];
        return
    end
    if i==Niterations
        disp("not converged")
    end
end
end