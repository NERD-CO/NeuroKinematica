%% Completed subjects/cases:
% 1: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\03_09_2023\Raw Electrophysiology MATLAB'
% 2: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\03_23_2023\Raw Electrophysiology MATLAB'
% 3: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\04_05_2023\Raw Electrophysiology MATLAB'
% 4: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\04_13_2023\Raw Electrophysiology MATLAB\LH'
% 5: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\04_13_2023\Raw Electrophysiology MATLAB\RH'
% 6: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\05_11_2023\Raw Electrophysiology MATLAB'
% 7: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\05_18_2023_a\Raw Electrophysiology MATLAB'
% 8: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\05_18_2023_b\Raw Electrophysiology MATLAB\LH'
% 9: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\05_18_2023_b\Raw Electrophysiology MATLAB\RH'
% 10: 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\05_31_2023\Raw Electrophysiology MATLAB'

%% Run function:

curPCname = getenv('COMPUTERNAME'); % for windows

switch curPCname
    case 'DESKTOP-I5CPDO7'
        Case_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\03_09_2023'; 
        IO_DataDir = 'X:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';  
    otherwise
        Case_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative\05_31_2023'; 
        IO_DataDir = 'Z:\RadcliffeE\Thesis_PD Neuro-correlated Kinematics\Data\Intraoperative';  
end

updateTrialInformation(IO_DataDir , 1 ,  Case_DataDir)
