clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\UnivarScatter-master'))

Muscles = {...
    'm1_2';'m2_3';'m3_4';'m4_5';'m5_6';'m6_7';'m8_9';'m9_10';'m10_11';'m11_12';'m12_13';'m13_14';'m15_16';'m16_17';'m17_18';...
    'm18_19';'m19_20';'m20_21';'m22_23';'m23_24';'m24_25';'m25_26';'m26_27';'m27_28';'m29_30';'m30_31';'m31_32';'m32_33';...
    'm33_34';'m34_35';'m36_37';'m37_38';'m38_39';'m39_40';'m40_41';'m41_42';'m43_44';'m44_45';'m45_46';'m46_47';'m47_48';...
    'm48_49';'Biceps';'Triceps';'DeltAnt';'DeltMed';'SupTrap'};

cd(['C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano'])
load(['Group1_Ha.mat'])
load(['Group2_Ha.mat'])

RPE_G1 = Group1_Ha.(Muscles{1})(Group1_Ha.(Muscles{1})(:,3)==1,[1 3:4]);
RPE_G2 = Group2_Ha.(Muscles{1})(Group2_Ha.(Muscles{1})(:,3)==1,[1 3:4]);

part1 = unique(Group1_Ha.m1_2(:,1));
part2 = unique(Group2_Ha.m1_2(:,1));

RPE_Group1 = [];
for iG1 = 1:length(part1)
    vec = [part1(iG1) 1 0];
    temp = RPE_G1(RPE_G1(:,1)==part1(iG1),:);
    RPE_Group1 = [RPE_Group1 ; vec ; temp];
end
RPE_Group2 = [];
for iG2 = 1:length(part2)
    vec = [part2(iG2) 1 0];
    temp = RPE_G2(RPE_G2(:,1)==part2(iG2),:);
    RPE_Group2 = [RPE_Group2 ; vec ; temp];
end

RPE_Normalized_G1 = [];
for iPG1 = 1:length(part1)
    temp = RPE_Group1(RPE_Group1(:,1)==part1(iPG1),:);
    RPE_Normalized_G1(:,iPG1) = interp1(1:length(temp),temp(:,3),linspace(1,length(temp),100)) ;
end

RPE_Normalized_G2 = [];
for iPG2 = 1:length(part2)
    temp = RPE_Group2(RPE_Group2(:,1)==part2(iPG2),:);
    RPE_Normalized_G2(:,iPG2) = interp1(1:length(temp),temp(:,3),linspace(1,length(temp),100)) ;
end
RPE_Normalized = [RPE_Normalized_G1 RPE_Normalized_G2];

%% Stats 2 facteurs
Intervals = [1 10 20 30 40 50 60 70 80 90 100];
RPE_Normalized_Intervals(1,:) = mean(RPE_Normalized(1:10,:));
vInc = 2;
for iInt = 2:2:10
    RPE_Normalized_Intervals(vInc,:) = mean(RPE_Normalized(Intervals(iInt:iInt+1),:));
    vInc = vInc+1;
end

for iRep = 1:6
    within_factors(:,iRep) = RPE_Normalized_Intervals(iRep,:);
end
between_factors = [repmat(1,26,1); repmat(2,23,1)];

t = table([part1; part2],within_factors(:,1),within_factors(:,2),within_factors(:,3),within_factors(:,4),...
    within_factors(:,5),within_factors(:,6),between_factors,'VariableNames',{'Part','MF1','MF20',...
    'MF40','MF60','MF80','MF100','Group'});
t.Group = categorical(t.Group);
Meas = table([1 2 3 4 5 6]','VariableNames',{'Measurement'});
rm = fitrm(t,'MF1-MF100~Group','WithinDesign',Meas);
ranovatbl = ranova(rm);
anovatbl = anova(rm);

    % comparaison de groupe à chaque instant
    for iT = 1:6
        t = table(within_factors(:,iT),between_factors,'VariableNames',{'RPE','Group'});
        t.Group = categorical(t.Group);
        [P,ANOVATAB,STATS] = anova1(within_factors(:,iT),between_factors)
        ranovatbl = ranova(rm);
        posthoc = multcompare(rm,'Measurement');
        PostHocResults = posthoc([1 7 13 19 25],:);
        [QTime] = mafdr(table2array(PostHocResults(:,5)),'BHFDR',true);
        computeCohen_d(table2array(t(1:STATS.n(1),1)),table2array(t(STATS.n(2):end,1)),'independent')
        pause()
    end
    
%% Post-hoc
for iG = 1%:2
    RPE_Normalized = [];
    within_factors = [];
    between_factors = [];
    if iG == 1
        for iRep = 1:6
            within_factors(:,iRep) = RPE_Normalized_Intervals(iRep,1:length(part1));
        end
        Colors = [0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410];
        UnivarScatter(within_factors, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
        Part = part1;
        between_factors = [repmat(1,length(part1),1)];
    else
        for iRep = 1:6
            within_factors(:,iRep) = RPE_Normalized_Intervals(iRep,length(part1)+1:end);
        end
        Colors = [0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980];
        UnivarScatter(within_factors, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
        Part = part2;
        between_factors = [repmat(1,length(part2),1)];
    end

    % par groupe pour l'effet de fatigue
    t = table([Part],within_factors(:,1),within_factors(:,2),within_factors(:,3),within_factors(:,4),...
        within_factors(:,5),within_factors(:,6),between_factors,'VariableNames',{'Part','MF1','MF20',...
        'MF40','MF60','MF80','MF100','Group'});
    t.Group = categorical(t.Group);
    Meas = table([1 2 3 4 5 6]','VariableNames',{'Measurement'});
    rm = fitrm(t,'MF1-MF100~Part','WithinDesign',Meas);
    ranovatbl = ranova(rm);
    posthoc = multcompare(rm,'Measurement');
    PostHocResults = posthoc([1 7 13 19 25],:);
    [QTime] = mafdr(table2array(PostHocResults(:,5)),'BHFDR',true);
    computeCohen_d(table2array(t(1:end,6)),table2array(t(1:end,7)),'paired')

    pause()
end
% Plot

RPE1 = [];
for iP = 1:length(part1)
    RPE1 = [RPE1; RPE_Normalized_G1(:,iP)];
end
RPE1(:,2) = repmat([1:100],1,length(part1));

RPE2 = [];
for iP = 1:length(part2)
    RPE2 = [RPE2; RPE_Normalized_G2(:,iP)];
end
RPE2(:,2) = repmat([1:100],1,length(part2));

[Moy,CI,STD] = grpstats(RPE1(:,1),RPE1(:,2),{'mean','meanci','std'});
plot_ci(1:length(Moy),[Moy,Moy-STD,Moy+STD],'PatchColor', 'b', 'PatchAlpha', 0.3, 'MainLineColor', 'b','LineWidth', 0.1, 'LineColor', 'w')
hold on
[Moy,CI,STD] = grpstats(RPE2(:,1),RPE2(:,2),{'mean','meanci','std'});
plot_ci(1:length(Moy),[Moy,Moy-STD,Moy+STD],'PatchColor', 'r', 'PatchAlpha', 0.3, 'MainLineColor', 'r','LineWidth', 0.1, 'LineColor', 'w')

axis([0 100 0 8])
title('mRPE over time') ; xlabel('Time (%)') ; ylabel('mRPE score')



