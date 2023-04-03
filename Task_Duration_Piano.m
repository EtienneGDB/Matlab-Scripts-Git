clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\UnivarScatter-master'))

%% Ha
cd(['J:\Piano_Fatigue\Data_Exported'])
load(['cycles_Ha.mat'])

cd(['Z:\Bureau\Etienne\Extracted data\Piano_Fatigue\Data\Piano'])
load(['Group1_Ha.mat'])
part1_Ha = unique(Group1_Ha.m1_2(:,1));

TimePlaying = [];
for iSubjects = 1:50
    TimePlaying(iSubjects,:) = participants(iSubjects).t_do(length(participants(iSubjects).t_do))-participants(iSubjects).t_do(1);
end

TimePlaying_Part1_Ha = TimePlaying(part1_Ha);
mean(TimePlaying_Part1_Ha)
std(TimePlaying_Part1_Ha)

%% Li
cd(['J:\Piano_Fatigue\Data_Exported'])
load(['cycles_Li.mat'])

cd(['Z:\Bureau\Etienne\Extracted data\Piano_Fatigue\Data\Piano'])
load(['Group1_Li.mat'])
part1_Li = unique(Group1.m1_2(:,1));

TimePlaying = [];
for iSubjects = [1:16 18:50]
    TimePlaying(iSubjects,:) = participants(iSubjects).t_do(length(participants(iSubjects).t_do))--participants(iSubjects).t_do(1);
end

TimePlaying_Part1_Li = TimePlaying(part1_Li);
TimePlaying_Part1_Li(27:30) = NaN;
nanstd(TimePlaying_Part1_Li)

%% Plot
Colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
figure; UnivarScatter([TimePlaying_Part1_Ha TimePlaying_Part1_Li], 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
ylabel('Task duration (s)')
xticklabels({'Digital', 'Chord'})
axis([0.5 2.5 0 500])



