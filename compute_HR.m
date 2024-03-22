% Determine participant's heartrate (HR) goal based on their age

Pt_age = 65; % input age for each participant

[Pt_maxHR, Pt_goalHR, Pt_HR75, Pt_HR85] = compute_HRs(Pt_age);


% Local function:
function [maxHR, goalHR, lowerLim_HR75, upperLim_HR85] = compute_HRs(age)
    maxHR = 220 - age;
    goalHR = 0.8 * maxHR;
    lowerLim_HR75 = 0.75 * maxHR;
    upperLim_HR85 = 0.85 * maxHR;
end
