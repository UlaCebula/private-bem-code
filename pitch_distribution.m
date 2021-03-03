function localpitch = pitch_distribution(r_R,collective_blade_twist)
    twist = -50.*r_R + 35; %local twist [deg]
    localpitch = twist + collective_blade_twist;
end

