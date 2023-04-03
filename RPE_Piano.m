clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))

Muscles = {...
    'm1_2';'m2_3';'m3_4';'m4_5';'m5_6';'m6_7';'m8_9';'m9_10';'m10_11';'m11_12';'m12_13';'m13_14';'m15_16';'m16_17';'m17_18';...
    'm18_19';'m19_20';'m20_21';'m22_23';'m23_24';'m24_25';'m25_26';'m26_27';'m27_28';'m29_30';'m30_31';'m31_32';'m32_33';...
    'm33_34';'m34_35';'m36_37';'m37_38';'m38_39';'m39_40';'m40_41';'m41_42';'m43_44';'m44_45';'m45_46';'m46_47';'m47_48';...
    'm48_49';'Biceps';'Triceps';'DeltAnt';'DeltMed';'SupTrap'};

cd(['\\10.89.24.15\e\Bureau\Etienne\Extracted data\Piano_Fatigue\Data\Piano'])
load(['Group1_Ha.mat'])
load(['Group2_Ha.mat'])
partGroup1 = unique(Group1_Ha.(Muscles{1})(:,1));
partGroup2 = unique(Group2_Ha.(Muscles{1})(:,1));

RPE_Group1 = Group1_Ha.(Muscles{1})(Group1_Ha.(Muscles{1})(:,3)==1,[1 3:4]);
RPE_Group2 = Group2_Ha.(Muscles{1})(Group2_Ha.(Muscles{1})(:,3)==1,[1 3:4]);

Baseline = [];
Finish = [];
Group = [];
for iGroup = 1:2
    Group_Selected = eval(['RPE_Group' num2str(iGroup)]);
    Participants = unique(Group_Selected(:,1));
    
    for iP = 1:length(Participants)
        temp = Group_Selected(Group_Selected(:,1)==Participants(iP),:);
        Baseline = [Baseline; temp(1,3)];
        Finish = [Finish; temp(end,3)];
        Group = [Group; iGroup];
    end
end

% Anova
within_factors(:,1) = Baseline;
within_factors(:,2) = Finish;
between_factors = Group;

boxplot(within_factors)

tbl = simple_mixed_anova(within_factors,between_factors,{'Time'},{'Group'})

%% RPE variability article
mean(Baseline(Group==1))
mean(Baseline(Group==2))

RPE_Termination = [];
Group = [];
for iGroup = 1:2
    Group_Selected = eval(['RPE_Group' num2str(iGroup)]);
    Participants = unique(Group_Selected(:,1));
    
    for iP = 1:length(Participants)
        temp = Group_Selected(Group_Selected(:,1)==Participants(iP),:);
        if iGroup == 1
            Term = length(temp);
        else
            Term = round(length(temp)/2);
        end
        RPE_Termination = [RPE_Termination; temp(Term,3)];
        Group = [Group; iGroup];
    end
end
mean(RPE_Termination(Group==1))
std(RPE_Termination(Group==1))

mean(RPE_Termination(Group==2))
std(RPE_Termination(Group==2))

[h,p,ci,stats] = ttest2(RPE_Termination(Group==1),RPE_Termination(Group==2))

