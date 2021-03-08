function localpitch = pitch_distribution(r_R,collective_blade_twist, propeller)
    if (propeller)
        twist = -50.*r_R + 35; %local twist [deg]
        localpitch = twist + collective_blade_twist;
    else
        twist = 14*(1-r_R);
        localpitch = twist + collective_blade_twist;
        
    end
end

