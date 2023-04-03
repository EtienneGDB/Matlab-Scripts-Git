clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\UnivarScatter-master'))

cd(['J:\Piano_Fatigue\Data_Exported'])
load(['cycles_Ha.mat'])

cd(['C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano'])
load(['Group1_Ha.mat'])
load(['Group2_Ha.mat'])

Muscles = {...
    'm1_2';'m2_3';'m3_4';'m4_5';'m5_6';'m6_7';'m8_9';'m9_10';'m10_11';'m11_12';'m12_13';'m13_14';'m15_16';'m16_17';'m17_18';...
    'm18_19';'m19_20';'m20_21';'m22_23';'m23_24';'m24_25';'m25_26';'m26_27';'m27_28';'m29_30';'m30_31';'m31_32';'m32_33';...
    'm33_34';'m34_35';'m36_37';'m37_38';'m38_39';'m39_40';'m40_41';'m41_42';'m43_44';'m44_45';'m45_46';'m46_47';'m47_48';...
    'm48_49';'Biceps';'Triceps';'DeltAnt';'DeltMed';'SupTrap'};

TimePlaying = [];
for iGroup = 1:2 ;
    Group_Selected = eval(['Group' num2str(iGroup) '_Ha']);... '_Ha']);
        
    part = unique(Group_Selected.(Muscles{1})(:,1));
    length(part)

    for iSubjects = 1:length(part)
        TimePlaying(iSubjects,iGroup) = participants(part(iSubjects)).t_do(length(participants(part(iSubjects)).t_do))-participants(part(iSubjects)).t_do(1);
    end
end
TimePlaying(TimePlaying(:,:)==0) = NaN;

Colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
UnivarScatter(TimePlaying/60, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black')
axis([0 2 0 14])

nanmean(TimePlaying(:,2))
nanstd(TimePlaying(:,2))
