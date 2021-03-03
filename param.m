function [polar, prop, oper, air] = param(spacing)
%PARAM: outputs geometry and airfoil/flow properties 

%Airfoil polar
data = readmatrix('ARAD8pct_polar.txt');
polar.alpha = data(:,1);
polar.Cl = data(:,2);
polar.Cd = data(:,3);
polar.Cm = data(:,4);


%Propeller geometry
prop.blade_root = 0.25; %non-dimensional blade starting location 
switch spacing
    case 'constant'
        %constant spacing
        prop.r_R = linspace(prop.blade_root,1.00,100);%constant spacing
        prop.r_R = prop.r_R(:);%make column vector
        prop.dr = prop.r_R(2:end)-prop.r_R(1:end-1); %non-dimensional length of each section
    case 'cosine'
        theta = linspace(0,pi,100);
        prop.r_R = prop.blade_root+0.5*(1-prop.blade_root).*(1-cos(theta));
        prop.r_R = prop.r_R(:);%make column vector
        prop.dr = prop.r_R(2:end)-prop.r_R(1:end-1);
    otherwise 
        disp('Error: Invalid spacing specified');
        return
end
prop.c_R = chord_distribution(prop.r_R);
prop.twist = -50*prop.r_R+35; %local twist
prop.collective_blade_twist = 46;% [degrees]
prop.pitch = pitch_distribution(prop.r_R, prop.collective_blade_twist);
prop.R = 0.70;%propeller radius [m]
prop.Nblades = 6;%number of blades

%Operating conditions
oper.U_inf = 60;% [m/s]
oper.n = 1200;% rotational speed [rpm]
oper.omega = (oper.n*2*pi)/(60); %rotational speed [rad/s]
oper.TSR = (oper.omega*prop.R)/oper.U_inf;

[air.Temp, air.speed_of_sound, air.pressure, air.density] = atmosisa(2000);
end

