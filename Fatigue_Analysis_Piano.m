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
    
FoldersNames = dir('J:\Piano_Fatigue\Data_Exported\EMG_Ha');
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

cd(['C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano'])
load(['Tri_TABLEDATA_Int.mat'])
Tri_TABLEDATA_Int = Tri;
load(['Tri_TABLEDATA_Int_LessSixMin.mat'])
Tri_TABLEDATA_Int_LessSixMin = Tri;
load(['Tri_TABLEDATA_Int_MoreSixMin.mat'])
Tri_TABLEDATA_Int_MoreSixMin = Tri;

load(['TABLEDATA_Int.mat'])
load(['TABLEDATA_Int_LessSixMin.mat'])
load(['TABLEDATA_Int_MoreSixMin.mat'])

load(['Group1.mat'])
load(['Group2.mat'])
load(['Group3.mat'])

% Plot
Slope = [];
for iM = 1:length(Muscles)
    [Moy,CI,STD] = grpstats(TABLEDATA_Int_MoreSixMin.(Muscles{iM})(:,4),TABLEDATA_Int_MoreSixMin.(Muscles{iM})(:,[6]),{'mean','meanci','std'});
    plot_ci(1:length(Moy),[Moy,Moy-STD,Moy+STD],'PatchColor', 'b', 'PatchAlpha', 0.3, 'MainLineColor', 'b','LineWidth', 0.1, 'LineColor', 'w')
    title((Muscles{iM}),'Interpreter','none') ; xlabel('Time') ; ylabel('RPE')
    hold on
    plot(TABLEDATA_Int_MoreSixMin.(Muscles{iM})(TABLEDATA_Int_MoreSixMin.(Muscles{iM})(:,1)==8,6),TABLEDATA_Int_MoreSixMin.(Muscles{iM})(TABLEDATA_Int_MoreSixMin.(Muscles{iM})(:,1)==8,12))
    
    % Calculate Slopes
%     [r,m,b] = regression([1:length(Moy)],Moy');
%     Slope(iM) = m;
%     hold on
%     plot([1:length(Moy)]',[1:length(Moy)]'*m+b)
    pause
    close all
end
Slope = atan(Slope)/90;
% title('RPE over time') ; xlabel('Time (min)') ; ylabel('RPE')
% xticklabels({'','1','','2','','3','','4','','5'})...LessSixMin
% xticklabels({'','1','2','3','4','5','6','7','8','9','10','11','12'})...MoreSixMin
% xticks([0:2:25])

% save('C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\Slope_MoreSixMin_Mobility.mat','Slope')

[Moy,CI] = grpstats(TABLEDATA.(Muscles{iM})(:,11),TABLEDATA.(Muscles{iM})(:,[6]));
plot(Moy)
hold on
plot_ci(1:length(Moy),[Moy,Moy-CI,Moy+CI],'PatchColor', 'b', 'PatchAlpha', 0.3, 'MainLineColor', 'b','LineWidth', 0.1, 'LineColor', 'w')

% MAP
% Color Map Pvalues linear-mixed models
cd(['C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano\Piano_Fatigue'])
[num, txt, raw] = xlsread('TableResults_MeanMF.xlsx');
Pvalues = num(:,3);

seq=1:6:42;
t=[];
for i=1:7
        temp = Pvalues(seq(i):seq(i)+5);
        t(:,i) = temp;
end

imagesc(t)

% Color Map Slopes
cd(['C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano'])
load('Slope_MedFreq.mat');
Slope1 = Slope;
load('Slope_SpecEn.mat');
Slope2 = Slope;
load('Slope_Mobility.mat');
Slope3 = Slope;

seq = 1:6:length(Slope);
slope_MedFreq = [];
slope_SpecEn = [];
slope_Mobility = [];
for i = 1:7
    slope_MedFreq(:,i) = Slope1(seq(i):seq(i)+5)';
    slope_SpecEn(:,i) = Slope2(seq(i):seq(i)+5)';
    slope_Mobility(:,i) = Slope3(seq(i):seq(i)+5)';
end
min([Slope1 Slope2])
max([Slope1 Slope2])
imagesc(slope_MedFreq)
imagesc(slope_SpecEn)
imagesc(slope_Mobility)

% Regression
% Armonisation Tri
% Étape 1-2
MinSize = 1000;
for iM = 1:length(Muscles)
    Sz = size(Tri.(Muscles{iM}));
    if Sz(1) < MinSize
        MinSize = Sz(1);
        IdRef = iM;
    end
end
XRef = unique(Tri.(Muscles{IdRef})(:,1));
NRef = numel(XRef);
countRef = zeros(NRef,1);
for k = 1:NRef
    countRef(k) = sum(Tri.(Muscles{IdRef})(:,1)==XRef(k));
end
CTRef = [ XRef(:) countRef ];
    
for iM = 1:length(Muscles)
    X2 = unique(Tri.(Muscles{iM})(:,1));
    N2 = numel(X2);
    count2 = zeros(N2,1);
    for k = 1:N2
        count2(k) = sum(Tri.(Muscles{iM})(:,1)==X2(k));
    end
    CT2 = [ X2(:) count2 ];
    Vsync = setdiff(CT2,CTRef);
    if length(Vsync)~=0
        for iVsync = 1:length(Vsync)
            Tri.(Muscles{iM})(Tri.(Muscles{iM})(:,1)==Vsync(iVsync),:)=[];
        end
    end
    Vsync = setdiff(CT2(:,1),CTRef(:,1));
    if length(Vsync)~=0
        for iVsync = 1:length(Vsync)
            Tri.(Muscles{iM})(Tri.(Muscles{iM})(:,1)==Vsync(iVsync),:)=[];
        end
    end
end

% Matrix Muscles-RPE
MatVar = [];
for i = 1:length(struct2cell(Tri_TABLEDATA_Int_LessSixMin))
    MatVar(:,i) = Tri_TABLEDATA_Int_LessSixMin.(Muscles{i})(:,12);
end
MatVar2 = MatVar(:,[29 30 35 36 41 42]);

[b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(MatVar(:,:),Tri_TABLEDATA_Int_LessSixMin.(Muscles{1})(:,4),'penter',0.05,'premove',0.1);
Cste = ones(length(MatVar(:,inmodel)));
[b,bint,r,rint,stats] = regress(Tri_TABLEDATA_Int_LessSixMin.(Muscles{1})(:,4),[Cste(:,1) MatVar(:,inmodel)]);
stats
find(inmodel)

Cste = ones(length(MatVar2(:,1)));
[b,bint,r,rint,stats] = regress(Tri_TABLEDATA_Int_LessSixMin.(Muscles{1})(:,4),[Cste(:,1) MatVar2(:,:)]);
stats

plot(Tri_TABLEDATA_Int_LessSixMin.(Muscles{1})(:,2),Tri_TABLEDATA_Int_LessSixMin.(Muscles{1})(:,4),'.')
