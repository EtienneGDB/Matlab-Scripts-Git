clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))

FoldersNames = dir('J:\Piano_Fatigue\Data_Exported\EMG_Max_MVC');
MVCFilenames = {...
    FoldersNames(3:52).name...
   } ;
    %Crée une deuxième ligne avec le nom du fichier sans .mat
    MVCFilenames{2,1}=[];
    for iMVCFilenames = 1:length(MVCFilenames)
        MVCFilenames{2,iMVCFilenames} = erase(MVCFilenames{1,iMVCFilenames},'.mat');
    end
    
FoldersNames = dir('J:\Piano_Fatigue\Data_Exported\EMG_Li');
HaFilenames = {...
    FoldersNames(3:52).name...
   } ;
    %Crée une deuxième ligne avec le nom du fichier sans .mat
    HaFilenames{2,1}=[];
    for iHaFilenames = 1:length(HaFilenames)
        HaFilenames{2,iHaFilenames} = erase(HaFilenames{1,iHaFilenames},'.mat');
    end

Muscles = {...
    'm1_2';'m2_3';'m3_4';'m4_5';'m5_6';'m6_7';'m8_9';'m9_10';'m10_11';'m11_12';'m12_13';'m13_14';'m15_16';'m16_17';'m17_18';...
    'm18_19';'m19_20';'m20_21';'m22_23';'m23_24';'m24_25';'m25_26';'m26_27';'m27_28';'m29_30';'m30_31';'m31_32';'m32_33';...
    'm33_34';'m34_35';'m36_37';'m37_38';'m38_39';'m39_40';'m40_41';'m41_42';'m43_44';'m44_45';'m45_46';'m46_47';'m47_48';...
    'm48_49';'Biceps';'Triceps';'DeltAnt';'DeltMed';'SupTrap'};

cd(['X:\Bureau\Etienne\Extracted data\Piano_Fatigue\Data\Piano'])
load(['TABLEDATA_Int_Li.mat'])
load(['Group1_Li.mat'])
load(['Group2_Li.mat'])
% load(['Group3.mat'])

RPE(:,1) = TABLEDATA_Int.(Muscles{1})(TABLEDATA_Int.(Muscles{1})(:,3)~=0,1);
RPE(:,2) = TABLEDATA_Int.(Muscles{1})(TABLEDATA_Int.(Muscles{1})(:,3)~=0,4);
RPE(:,3) = 0;
for iP = 1:50
    temp = RPE(RPE(:,1)==iP,:);
    RPE(RPE(:,1)==iP,3) = [1:length(temp)]';
end

partGroup1 = unique(Group1.(Muscles{1})(:,1));
partGroup2 = unique(Group2.(Muscles{1})(:,1));
% partGroup3 = unique(Group3.(Muscles{1})(:,1));

RPE_Group1 = RPE(ismember(RPE(:,1),partGroup1),:);
RPE_Group2 = RPE(ismember(RPE(:,1),partGroup2),:);
for iP = 1:length(partGroup2)
    temp = RPE_Group2(RPE_Group2(:,1)==partGroup2(iP),:);
    RPEend(iP,:) = temp(end,2);
end
mean(RPEend)
std(RPEend)
% RPE_Group3 = RPE(ismember(RPE(:,1),partGroup3),:);

Group_Selected = RPE_Group1;

[Moy,CI,STD] = grpstats(RPE_Group1(:,2),RPE_Group1(:,3),{'mean','meanci','std'});
plot_ci(1:length(Moy),[Moy,Moy-STD,Moy+STD],'PatchColor', 'b', 'PatchAlpha', 0.3, 'MainLineColor', 'b','LineWidth', 0.1, 'LineColor', 'w')
hold on
[Moy,CI,STD] = grpstats(RPE_Group2(:,2),RPE_Group2(:,3),{'mean','meanci','std'});
plot_ci(1:length(Moy),[Moy,Moy-STD,Moy+STD],'PatchColor', 'r', 'PatchAlpha', 0.3, 'MainLineColor', 'r','LineWidth', 0.1, 'LineColor', 'w')
% hold on
% [Moy,CI,STD] = grpstats(RPE_Group3(:,2),RPE_Group3(:,3),{'mean','meanci','std'});
% plot_ci(1:length(Moy),[Moy,Moy-STD,Moy+STD],'PatchColor', 'g', 'PatchAlpha', 0.3, 'MainLineColor', 'g','LineWidth', 0.1, 'LineColor', 'w')

axis([0 26 0 8])
title('RPE over time') ; xlabel('Time (min)') ; ylabel('RPE')
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12'})
xticks([0:2:25])
