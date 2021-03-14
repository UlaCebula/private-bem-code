function [c_R] = chord_distribution(r_R, R, propeller)
%CHORD_DISTRIBUTION Summary of this function goes here
    %non-dimensional
    if (propeller)
        c_R = 0.18 - 0.06.*r_R;
    else
        c_R = (3*(1-r_R)+1)/R;
    end
end

