clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Entropy Fouaz'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\UnivarScatter-master'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\waterloo\Waterloo_MATLAB_Library'))

addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\pls_Cao'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\PLSI_Zhang'))

% addpath E:\Bureau\Etienne\MATLAB\Functions
% addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\btk'))
% addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))
% addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\Matlab-Functions-Git'))
% addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\FFT'))
% addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\Entropy Fouaz'))

Task = 'Li'; % changer Li2 et Ha2 partout ('Matrix dimension must agree')
% FoldersNames = dir(['J:\Piano_Fatigue\Data_Exported\Features_XSENS_Cycles\' Task]);
FoldersNames = dir(['\\10.89.24.15\j\Piano_Fatigue\Data_Exported\Features_XSENS\' Task]);
FilenamesXSENS = {...
    FoldersNames(3:end).name...
   } ;

%Crée une deuxième ligne avec le nom du fichier sans .mat
FilenamesXSENS{2,1}=[];
for iFilenames = 1:length(FilenamesXSENS)
    FilenamesXSENS{2,iFilenames} = erase(FilenamesXSENS{1,iFilenames},'.mat');
end

% cd(['J:\Piano_Fatigue\Data_Exported'])
% load(['Cycle_Ha_Oceane.mat'])

FoldersNames = dir(['\\10.89.24.15\j\Piano_Fatigue\Data_Exported\TFR_MedianFreq_' Task]);
FilenamesEMG = {...
    FoldersNames(3:end).name...
   } ;

cd(['\\10.89.24.15\j\Piano_Fatigue\Matlab_matrix\Info_participants'])
load(['Info_participants_corrected'])
if Task == 'Ha'
    Info_participants(22) = [];
elseif Task == 'Li'
    Info_participants(17) = [];
end
% for iP = 1:50
%     Demography(iP,1) = Info_participants(iP).Sexe;
%     Demography(iP,2) = Info_participants(iP).Age;
%     Demography(iP,3) = Info_participants(iP).Poids;
%     Demography(iP,4) = Info_participants(iP).nb_annees;
%     Demography(iP,5) = Info_participants(iP).nb_heures;
%     Demography(iP,6) = mean(Info_participants(iP).MVC);
% 
%     Lat(iP,:) = Info_participants(iP).Lateralite;
%     Sex(iP,:) = Info_participants(iP).Sexe;
% end
% sum(Lat(:,:),'omitnan')
% C = 6;
% mean(Demography(:,C),'omitnan')
% std(Demography(:,C),'omitnan')
% 
% xxx = Demography(:,C)/2,20462262;
% mean(xxx,'omitnan')
% std(xxx,'omitnan')

% cd(['C:\Users\p1098713\Documents\2.Post-Doc\Data\Piano'])
% load(['Group1_Li.mat'])
% load(['Group2_Li.mat'])
% part1 = unique(Group1.m1_2(:,1));
% part2 = unique(Group2.m1_2(:,1));

Muscles = {...
    'm1_2';'m2_3';'m3_4';'m4_5';'m5_6';'m6_7';'m8_9';'m9_10';'m10_11';'m11_12';'m12_13';'m13_14';'m15_16';'m16_17';'m17_18';...
    'm18_19';'m19_20';'m20_21';'m22_23';'m23_24';'m24_25';'m25_26';'m26_27';'m27_28';'m29_30';'m30_31';'m31_32';'m32_33';...
    'm33_34';'m34_35';'m36_37';'m37_38';'m38_39';'m39_40';'m40_41';'m41_42';'m43_44';'m44_45';'m45_46';'m46_47';'m47_48';...
    'm48_49';'Biceps';'Triceps';'DeltAnt';'DeltMed';'SupTrap'};

XSENSs = {'Angular_Velocity'; 'Acceleration'; 'Module_Angular_Velocity'; 'Module_Acceleration'; 'Jerk'; 'Module_Jerk'};
Features = {'MedianFreq'; 'SpectralEntropy'; 'PowerLF'; 'PowerHF'; 'PowerTot'; 'PeakPower'; 'PeakPower_Freq'; 'Peak';...
    'Mean'; 'Median'; 'STD'; 'Percentile10'; 'Percentile25'; 'Percentile75'; 'Percentile90'};
Membres = {'L5';'T8';'Head';'Shoulder';'Arm';'Forearm';'Hand'};

% Label figure
% XSENSs = {'Angular Velocity'; 'Acceleration'; 'Module Angular Velocity'; 'Module Acceleration'; 'Jerk'; 'Module Jerk'};
% Features = {'Median Frequency'; 'Spectral Entropy'; 'Power Low Frequency'; 'Power High Frequency'; 'Total Power';...
%     'PeakPower'; 'PeakPower Frequency'; 'Maximum'; 'Mean'; 'Median'; 'STD'; 'Percentile10'; 'Percentile25';...
%     'Percentile75'; 'Percentile90'};
% Membres = {'L5';'T8';'Head';'Shoulder';'Arm';'Forearm';'Hand'};
% xticklabels([Features])
% xtickangle(45)
% xticklabels([Features{12}; Features{13}; Features{14}; Features{15}])
% xtickangle(45)
% xticklabels([XSENSs])
% xtickangle(45)
% xticklabels([Membres])
% xtickangle(45)

% For VarNames
XSENSsVarNames = {'AngularVelocity'; 'Acceleration'; 'AngularVelocityMag'; 'AccelerationMag'; 'Jerk'; 'JerkMag'};
FeaturesVarNames = {'MedianFreq'; 'SpectralEntropy'; 'PowerLF'; 'PowerHF'; 'PowerTot'; 'PeakPower'; 'PeakPowerFreq'; 'Maximum';...
    'Mean'; 'Median'; 'STD'; 'Percentile10'; 'Percentile25'; 'Percentile75'; 'Percentile90'};

for iFeature = 1:length(Features)
    Feature = Features{iFeature};
    for iXSENS = 1:length(XSENSs)
        XSENS = XSENSs{iXSENS};
        assignin('base',['Data_' Feature '_' XSENS],[])
%         assignin('base',['fullData_' Feature '_' XSENS],[])
    end
end
FreqMed_All=[];

%% Load 2 groups
cd(['\\10.89.24.15\e\Bureau\Etienne\Extracted data\Piano_Fatigue\Data\Piano'])
load(['Group1_Li.mat'])
load(['Group2_Li.mat'])
part1_Li = unique(Group1.m1_2(:,1));

load(['Group1_Ha.mat'])
load(['Group2_Ha.mat'])
part1_Ha = unique(Group1_Ha.m1_2(:,1));
part2_Ha = unique(Group2_Ha.m1_2(:,1));

if Task == 'Li'
    Part1 = part1_Li;
else
    Part1 = part1_Ha;
end

%% Extrait les données normalisées de variation dans les tables
for iP = Part1'
    VarNames = {};
%     if Task == 'Ha' & iP == 22
%     elseif Task == 'Li' & iP ==17
%     else
    if Task == 'Ha'
        RPE = Info_participants(iP).Hanon;
    else
        RPE = Info_participants(iP).Liszt;
    end
    idNaN = find(isnan(RPE));
    if sum(isnan(RPE)) >= 1 & idNaN(end) == length(RPE)
        idNaN = find(isnan(RPE));
        while idNaN(end) == length(RPE)
            RPE(end) = [];
            if sum(isnan(RPE)) >= 1
                idNaN = find(isnan(RPE));
            else
                idNaN = 0;
            end
        end
    elseif sum(isnan(RPE)) >= 1 & idNaN(end) ~= length(RPE)
        RPE(idNaN) = mean(RPE([idNaN-1 idNaN+1]));
    end
%     RPE = rmmissing(RPE);
%     RPE = [0; RPE];
%     normalize_RPE = interp1(1:length(RPE),RPE,linspace(1,length(RPE),100))';
    
%     % EMG FreqMed
%     cd(['J:\Piano_Fatigue\Data_Exported\TFR_MedianFreq_' Task(1:2)])
%     load([FilenamesEMG{1,iP}])
%     EMG_FreqMed = [];
%     for iM = 25:36  % Les muscles les plus affectés pour Li
%         temp = TFR_MedianFreq.MedianFreq.(Muscles{iM});
%         normalize_temp = interp1(1:length(temp),temp,linspace(1,length(temp),100))';
%         EMG_FreqMed = [EMG_FreqMed normalize_temp];        
%     end
%     freqmed = mean(EMG_FreqMed,2);
%     freqmednorm = (freqmed-mean(freqmed))/(std(freqmed));
%     FreqMed = [repmat(iP,100,1) freqmed freqmednorm];
%     FreqMed_All = [FreqMed_All; FreqMed];
    
%   % XSENS raw
    cd(['\\10.89.24.15\j\Piano_Fatigue\Data_Exported\Features_XSENS\' Task])
    load([FilenamesXSENS{1,iP}])
    for iFeature = 1:length(Features)
        Feature = Features{iFeature};
        FeatureVarNames = FeaturesVarNames{iFeature};
        for iXSENS = 1:length(XSENSs)
            VNames = {};
            XSENS = XSENSs{iXSENS};
            XSENSVarNames = XSENSsVarNames{iXSENS};
            normalize_temp_Sum_reduc = [];
            temp = Features_XSENS.(Feature).(XSENS);
            t=[0:length(temp(:,1))-1]*(1/60);
            normalize_temp = temp - temp(1,:);
            normalize_temp_Sum = normalize_temp;
            if size(temp,2) == 21
                VarInc = 1;
                for iC = 1:3:size(normalize_temp,2)
                    normalize_temp_Sum(:,iC) = normalize_temp_Sum(:,iC)/mean(maxk(normalize_temp_Sum(:,iC),10));
                    normalize_temp_Sum(:,iC+1) = normalize_temp_Sum(:,iC+1)/mean(maxk(normalize_temp_Sum(:,iC+1),10));
                    normalize_temp_Sum(:,iC+2) = normalize_temp_Sum(:,iC+2)/mean(maxk(normalize_temp_Sum(:,iC+2),10));
                    VNames{iC} = [FeatureVarNames '_' XSENSVarNames '_' Membres{VarInc} '_X'];
                    VNames{iC+1} = [FeatureVarNames '_' XSENSVarNames '_' Membres{VarInc} '_Y'];
                    VNames{iC+2} = [FeatureVarNames '_' XSENSVarNames '_' Membres{VarInc} '_Z'];
                    VarInc = VarInc +1;
                end
            else
                normalize_temp_Sum = normalize_temp;
                for iC = 1:size(normalize_temp_Sum,2)
                    normalize_temp_Sum(:,iC) = normalize_temp_Sum(:,iC)/mean(maxk(normalize_temp_Sum(:,iC),10));
                    VNames{iC} = [FeatureVarNames '_' XSENSVarNames '_' Membres{iC} '_M'];
                end
            end
            incDel = 0;
            for iRPE = 1:length(RPE)
                ValueRPE = iRPE*30;
                if ValueRPE+2 > length(normalize_temp_Sum)
                    incDel = incDel + 1;
                else
                    normalize_temp_Sum_reduc(iRPE,:) = median(normalize_temp_Sum(ValueRPE-2:ValueRPE+2,:));
                end
            end
%             fullData_normalize_temp_Sum = interp1(1:length(normalize_temp_Sum),normalize_temp_Sum,linspace(1,length(normalize_temp_Sum),100));
%             fullData_normalize_temp_Sum = [repmat(iP,length(fullData_normalize_temp_Sum),1) fullData_normalize_temp_Sum];
            normalize_temp_Sum = [repmat(iP,length(RPE)-incDel,1) RPE(1:length(RPE)-incDel) normalize_temp_Sum_reduc];

%    % XSENS cycles
%     cd(['J:\Piano_Fatigue\Data_Exported\Features_XSENS_Cycles\' Task])
%     load([FilenamesXSENS{1,iP}])
%     for iFeature = 1:length(Features)
%         Feature = Features{iFeature};
%         for iXSENS = 1:length(XSENSs)
%             VNames = {};
%             XSENS = XSENSs{iXSENS};
%             temp = Features_XSENS_Cycles.(Feature).(XSENS);
%             normalize_temp = interp1(1:length(temp),temp,linspace(1,length(temp),100));
%             normalize_temp = normalize_temp - normalize_temp(1,:);
%             normalize_temp_Sum = normalize_temp;
%             if size(temp,2) == 21
%                 VarInc = 1;
%                 for iC = 1:3:size(normalize_temp,2)
%                     normalize_temp_Sum(:,iC) = normalize_temp_Sum(:,iC)/mean(maxk(normalize_temp_Sum(:,iC),10));
%                     normalize_temp_Sum(:,iC+1) = normalize_temp_Sum(:,iC+1)/mean(maxk(normalize_temp_Sum(:,iC+1),10));
%                     normalize_temp_Sum(:,iC+2) = normalize_temp_Sum(:,iC+2)/mean(maxk(normalize_temp_Sum(:,iC+2),10));
%                     VNames{iC} = [Feature '_' XSENS '_' Membres{VarInc} '_X'];
%                     VNames{iC+1} = [Feature '_' XSENS '_' Membres{VarInc} '_Y'];
%                     VNames{iC+2} = [Feature '_' XSENS '_' Membres{VarInc} '_Z'];
%                     VarInc = VarInc +1;
%                 end
%             else
%                 normalize_temp_Sum = normalize_temp;
%                 for iC = 1:size(normalize_temp_Sum,2)
%                     normalize_temp_Sum(:,iC) = normalize_temp_Sum(:,iC)/mean(maxk(normalize_temp_Sum(:,iC),10));
%                     VNames{iC} = [Feature '_' XSENS '_' Membres{iC} '_M'];
%                 end
%             end
%             normalize_temp_Sum = [repmat(iP,100,1) normalize_RPE normalize_temp_Sum];

            VarNames = [VarNames,VNames];
%             fullData_Concatenate = [eval(['fullData_' Feature '_' XSENS]); fullData_normalize_temp_Sum];
            Concatenate = [eval(['Data_' Feature '_' XSENS]); normalize_temp_Sum];
            assignin('base',['Data_' Feature '_' XSENS],Concatenate)
%             assignin('base',['fullData_' Feature '_' XSENS],fullData_Concatenate)
        end
    end
%     end
end
% cd(['\\10.89.24.15\j\Piano_Fatigue\Data_Exported'])
% save('Data_STD_Jerk.mat',"Data_STD_Jerk")

%% Ordonne les données de variables dépendantes et indépendante dans des tables Z-values
% Tous les participants sont dans ces tables concaténées
Xnorm_allParticipants = [];
Ynorm_allParticipants = [];
Xnorm_allParticipants_Timenorm = [];
% fullXnorm_allParticipants = [];
for iP = Part1'
    X = [];
%     fullX = [];
    Xnorm = [];
    Labels = [];
%     if Task == 'Ha' & iP == 22
%     elseif Task == 'Li' & iP == 17
%     else
        for iFeature = 1:length(Features)
            Feature = Features{iFeature};
            for iXSENS = 1:length(XSENSs)
                XSENS = XSENSs{iXSENS};
                dat = eval(['Data_' Feature '_' XSENS]);
%                 fulldat = eval(['fullData_' Feature '_' XSENS]);
%                 for iMembres = 1:length(Membres)
%                     Lab{iMembres} = [Feature '_' XSENS '_' Membres{iMembres}];
%                 end
%                 Labels = [Labels Lab];

                temp = dat(dat(:,1)==iP,:);
%                 fulltemp = fulldat(fulldat(:,1)==iP,:);
                X = [X temp(:,3:end)];
%                 fullX = [fullX fulltemp(:,2:end)];
            end
        end
        if length(X) > 0
        Y = temp(:,2);
        
%         % Sub-class Y
%         Y(Y < 2) = 1;
%         Y(Y < 3 & Y > 1.5) = 2;
%         Y(Y < 5 & Y > 2.5) = 3;
%         Y(Y < 7 & Y > 4) = 4;
%         Y(Y > 6) = 5;
        
        % Normalize X and Y
        for iC = 1:size(X,2)
            Xnorm(:,iC) = (X(:,iC)-mean(X(:,iC)))/(std(X(:,iC))); % zcore
        end
        Xnorm(:,find(isnan(Xnorm(1,:))))=0;
        Ynorm = (Y-mean(Y))/(std(Y)); % zscore
%         if sum(isnan(Ynorm))>=1
%             break
%         end
%         fullXnorm_allParticipants = [fullXnorm_allParticipants; fullX];
        Xnorm_Timenorm = interp1(1:size(Xnorm,1),Xnorm,linspace(1,size(Xnorm,1),100));
        Xnorm_allParticipants_Timenorm = [Xnorm_allParticipants_Timenorm; Xnorm_Timenorm];
        
        Xnorm_allParticipants = [Xnorm_allParticipants; Xnorm];
        Ynorm_allParticipants = [Ynorm_allParticipants; Ynorm];
        end
%     end
end
% Xnorm_allParticipants_Li = Xnorm_allParticipants;
% Ynorm_allParticipants_Li = Ynorm_allParticipants;
% cd(['C:\Users\p1098713\Documents'])
% save(['Xnorm_allParticipants_' Task '.mat'],['Xnorm_allParticipants_' Task])
% save(['Ynorm_allParticipants_' Task '.mat'],['Ynorm_allParticipants_' Task])
% save(['VarNames_' Task '.mat'],['VarNames'])

%% Define training set and test set
MRMSE_IT = [];
MABSERROR_IT = [];
fullModel_VAREXPLAINED_IT = [];
fullModel_ABSERROR_IT = [];
ABSE_Original_IT = [];
ncomp_opt = [];
MABSERROR = [];
VIPSCORE_IT = [];
MRMSE = [];
for iIT = 1:10
Participants = [1:length(Part1)]';
k = 5;
c = cvpartition(length(Participants),'KFold',k);
% cd(['\\10.89.24.15\j\Piano_Fatigue\Data_Exported'])
% save(['cPartitionFolds_Revision_' Task '.mat'],['c'])

%% First AND second iteration with ncomp = 1:100 AND best ncomp
% Here VIP threshold does not vary ! VIP thresh = 1
for ncomp = 1:3
RMSE = [];
PRESS = [];
ABSERROR = [];
BETA = [];
VAREXPLAINED = [];
INDVIP = {};
VIPSCORE = [];
for ik = 1:k
    Yobs = [];
    idx = training(c,ik);
    idx = [true false true false true true false true ...
        true true true true false true true true true true true true false true true true true true];
    Part_Training_Set = Participants(idx);
    Part_Test_Set = Participants(idx==0);
    
    X_Training_Set = Xnorm_allParticipants(ismember(Data_STD_Jerk(:,1),Part1(Part_Training_Set)),:);
    X_Test_Set = Xnorm_allParticipants(ismember(Data_STD_Jerk(:,1),Part1(Part_Test_Set)),:);

    Y_Training_Set = Ynorm_allParticipants(ismember(Data_STD_Jerk(:,1),Part1(Part_Training_Set)),:);
    Y_Test_Set = Ynorm_allParticipants(ismember(Data_STD_Jerk(:,1),Part1(Part_Test_Set)),:);

    Row_Part_Training_Set = Data_STD_Jerk(ismember(Data_STD_Jerk(:,1),Part1(Part_Training_Set)),1);
    Row_Part_Test_Set = Data_STD_Jerk(ismember(Data_STD_Jerk(:,1),Part1(Part_Test_Set)),1);

    %% Apply PLSR on training set
%     ncomp = 2; % use only 10 components ? change this to desired number of components
    [x_loadings, y_loadings, x_scores, y_scores, beta, pctvar, mse, stats] = plsregress(X_Training_Set, Y_Training_Set, ncomp);
%     [x_loadings, y_loadings, x_scores, y_scores, beta, pctvar, mse, stats] = plsregress(X_Training_Set(:,Releavant_Var), Y_Training_Set, ncomp);
%     subplot(5,2,ik); plot(1:ncomp,cumsum(100*pctvar(2,:)),'-bo');
%     xlabel('Components')
%     ylabel('% Var Exp in Y')
%     title(['iteration' num2str(ik)])
    CumVar = cumsum(100*pctvar(2,:));
    VAREXPLAINED = [VAREXPLAINED; CumVar];
    
%     scatter(1:length(beta),beta,'filled')
%     xlabel('Variables')
%     ylabel('Beta coefficient')
    % hist(beta)
    BETA = [BETA beta];
    
    %         Ynorm = (Y-mean(Y))/(std(Y)); % zscore
    for iP = 1:length(Part_Test_Set)
        Yobs = [Yobs; dat(dat(:,1)==Part1(Part_Test_Set(iP)),2)];
    end
%     % Sub-class Yobs
%     Yobs(Yobs < 2) = 1;
%     Yobs(Yobs < 3 & Yobs > 1.5) = 2;
%     Yobs(Yobs < 5 & Yobs > 2.5) = 3;
%     Yobs(Yobs < 7 & Yobs > 4) = 4;
%     Yobs(Yobs > 6) = 5;
    
    yfitnorm = [ones(size(X_Test_Set,1),1) X_Test_Set]*beta;
%     yfitnorm = [ones(size(X_Test_Set,1),1) X_Test_Set(:,Releavant_Var)]*beta;
    yfit = (yfitnorm*std(Yobs))+mean(Yobs);
%     % Sub-class Yfit
%     yfit(yfit < 2) = 1;
%     yfit(yfit < 3 & yfit > 1.5) = 2;
%     yfit(yfit < 5 & yfit > 2.5) = 3;
%     yfit(yfit < 7 & yfit > 4) = 4;
%     yfit(yfit > 6) = 5;
    
    residuals = Yobs - yfit;

%     yfit = [ones(size(X_Test_Set,1),1) X_Test_Set]*beta;
%     residuals = Y_Test_Set - yfit;
%     stem(residuals)
%     xlabel('Observation');
%     ylabel('Residual');
    % RMSE
    squaredResiduals = residuals.^2;
    Mean_squaredResiduals = mean(squaredResiduals);
    rmse = sqrt(Mean_squaredResiduals);
    RMSE = [RMSE; rmse];
    
    press = sum(squaredResiduals);
    PRESS = [PRESS; press];
    
    abserror = mean(abs(Yobs - yfit));
    ABSERROR = [ABSERROR; abserror];
    
    % plot(x_scores*y_loadings')
    % hold on 
%     plot(yfit)
%     hold on 
%     plot(Y_Test_Set)
%     xlabel('Participant - Time');
%     ylabel('Normalized RPE');
%     pause()

    % VIP
    W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
    p = size(x_loadings,1);
    sumSq = sum(x_scores.^2,1).*sum(y_loadings.^2,1);
    vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
    indVIP = find(vipScore > 1);
    INDVIP{ik} = indVIP;
    VIPSCORE(:,ik) = vipScore;
    
%     scatter(1:length(vipScore),vipScore,'x')
%     hold on
%     scatter(indVIP,vipScore(indVIP),'rx')
%     plot([1 length(vipScore)],[1 1],'--k')
%     hold off
%     axis tight
%     xlabel('Predictor Variables')
%     ylabel('VIP Scores')

end
MRMSE(:,ncomp) = RMSE;
MABSERROR(:,ncomp) = ABSERROR;
fullModel_VAREXPLAINED = VAREXPLAINED;
fullModel_ABSERROR = ABSERROR;
ABSE_Original = mean(ABSERROR);
std(ABSERROR);
end
MRMSE_IT(iIT,:) = mean(MRMSE);
MABSERROR_IT(iIT,:) = mean(MABSERROR);
fullModel_VAREXPLAINED_IT(iIT,:) = mean(fullModel_VAREXPLAINED);
fullModel_ABSERROR_IT(iIT,:) = mean(MABSERROR);
ABSE_Original_IT(iIT,:) = ABSE_Original;
ncomp_opt(iIT) = find(mean(MABSERROR) == min(mean(MABSERROR)));
VIPSCORE_IT(:,iIT) = mean(VIPSCORE,2);
end
ncomp_opt(ncomp_opt == 2)

%% Plot full model
if Task == 'Ha'
    Colors = repmat([0 0.4470 0.7410],20,1);
else
    Colors = repmat([0.8500 0.3250 0.0980],100,1);
end
[p,tbl,stats] = anova1(fullModel_VAREXPLAINED_IT);
multcomp = multcompare(stats);
figure; UnivarScatter(fullModel_VAREXPLAINED_IT, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
for iNcomp = 1:size(fullModel_VAREXPLAINED_IT,2)
    iMod = find(multcomp(:,1)==iNcomp & multcomp(:,2)==iNcomp+1);
    if multcomp(iMod,6) < 0.05
        hold on;
        plot(iNcomp+1, 105, '*', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
    end
end
xlabel('Number of latent variables')
ylabel('RPE explained (%)')
axis([0 21 0 105])

[p,tbl,stats] = anova1(MABSERROR_IT);
multcomp = multcompare(stats);
multcomp = multcomp(multcomp(:,1)==2,:);
figure; UnivarScatter(MABSERROR_IT, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
for iMod = 1:length(multcomp)
    if multcomp(iMod,6) < 0.05
        hold on;
        plot(iMod+1, 2, '*', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
    end
end
xlabel('Number of latent variables')
ylabel('Absolute error')
axis([0 21 0.8 2.4])

XLat = find(mean(MABSERROR_IT) == min(mean(MABSERROR_IT)));
mean(fullModel_VAREXPLAINED_IT(:,XLat))
std(fullModel_VAREXPLAINED_IT(:,XLat))

mean(MABSERROR_IT(:,XLat))
std(MABSERROR_IT(:,XLat))

%%
sum_INDVIP = zeros(size(X_Training_Set,2),1);
for ivar = 1:size(X_Training_Set,2)
    for ik = 1:k
        y = sum(INDVIP{ik}==ivar);
        sum_INDVIP(ivar,1) = sum_INDVIP(ivar,1) + y;
    end
end
Releavant_Var = find(sum_INDVIP > 1);

Mean_VIPSCORE = mean(VIPSCORE_IT,2);

VarNames_dec = [];
for iV = 1:length(VarNames)
    temp = strsplit(VarNames{iV},'_');            
    VarNames_dec = [VarNames_dec temp'];
end

for iMembres = 1:length(Membres)
    VIPMembres{1,iMembres} = mean(VIPSCORE_IT(strcmp(VarNames_dec(3,:),Membres{iMembres}),:),2);
    VarNames_dec_Membres{1,iMembres} = VarNames_dec(:,strcmp(VarNames_dec(3,:),Membres{iMembres}));
end
for iVN = 1:length(VarNames_dec_Membres{1})
    C = {VarNames_dec_Membres{1}{1,iVN},VarNames_dec_Membres{1}{2,iVN},VarNames_dec_Membres{1}{4,iVN}};
    VarNames_dec_save{iVN,1} = C{1};
    VarNames_dec_save{iVN,2} = C{2};
    VarNames_dec_save{iVN,3} = C{3};
end
% save(['VIPMembres_Review_' Task '.mat'],['VIPMembres'])
% save(['VarNames_dec_save_Review_' Task '.mat'],['VarNames_dec_save'])

max10 = maxk(Mean_VIPSCORE,10);
for iS = 1:10
    idmax10(iS) = find(Mean_VIPSCORE == max10(iS));
end
VarNames(idmax10)

VarNames_Releavant_Var_dec = {};
VIPSCORE_Releavant_Var = VIPSCORE(Releavant_Var,:);
VarNames_Releavant_Var = VarNames(Releavant_Var);
for iV = 1:length(VarNames_Releavant_Var)
    temp = strsplit(VarNames_Releavant_Var{iV},'_');            
    VarNames_Releavant_Var_dec = [VarNames_Releavant_Var_dec temp'];
end
for iFeatures = 1:length(FeaturesVarNames)
    OccurenceFeatures(iFeatures) = sum(strcmp(VarNames_Releavant_Var_dec(1,:),FeaturesVarNames{iFeatures}));
    VIPFeatures{1,iFeatures} = mean(VIPSCORE_Releavant_Var(strcmp(VarNames_Releavant_Var_dec(1,:),FeaturesVarNames{iFeatures}),:),2);
end
Conc_VIPFeatures = padcat(VIPFeatures{1},VIPFeatures{2},VIPFeatures{3},VIPFeatures{4},VIPFeatures{5},VIPFeatures{6},VIPFeatures{7},VIPFeatures{8},VIPFeatures{9},VIPFeatures{10},VIPFeatures{11},VIPFeatures{12},VIPFeatures{13},VIPFeatures{14},VIPFeatures{15});
[p,tbl,stats] = anova1(Conc_VIPFeatures);
multcomp = multcompare(stats);
for iXSENSs= 1:length(XSENSsVarNames)
    OccurenceXSENSs(iXSENSs) = sum(strcmp(VarNames_Releavant_Var_dec(2,:),XSENSsVarNames{iXSENSs}));
    VIPXSENSs{1,iXSENSs} = mean(VIPSCORE_Releavant_Var(strcmp(VarNames_Releavant_Var_dec(2,:),XSENSsVarNames{iXSENSs}),:),2);
end
Conc_VIPXSENSs = padcat(VIPXSENSs{1},VIPXSENSs{2},VIPXSENSs{3},VIPXSENSs{4},VIPXSENSs{5},VIPXSENSs{6});
anova1(Conc_VIPXSENSs)
for iMembres = 1:length(Membres)
    OccurenceMembres(iMembres) = sum(count(VarNames_Releavant_Var_dec(3,:),Membres{iMembres}));
    VIPMembres{1,iMembres} = mean(VIPSCORE_Releavant_Var(strcmp(VarNames_Releavant_Var_dec(3,:),Membres{iMembres}),:),2);
end
Conc_VIPMembres = padcat(VIPMembres{1},VIPMembres{2},VIPMembres{3},VIPMembres{4},VIPMembres{5},VIPMembres{6},VIPMembres{7});
anova1(Conc_VIPMembres)

Colors = repmat([0 0.4470 0.7410],length(FeaturesVarNames),1);
% Colors = repmat([0.8500 0.3250 0.0980],length(FeaturesVarNames),1);
figure; UnivarScatter(Conc_VIPFeatures, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
xlabel('Features')
ylabel('VIP')
axis([0 16 0.5 4])

Colors = repmat([0 0.4470 0.7410],length(XSENSsVarNames),1);
% Colors = repmat([0.8500 0.3250 0.0980],length(XSENSsVarNames),1);
figure; UnivarScatter(Conc_VIPXSENSs, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
xlabel('Variables')
ylabel('VIP')
axis([0 16 0.5 4])

Colors = repmat([0 0.4470 0.7410],length(Membres),1);
% Colors = repmat([0.8500 0.3250 0.0980],length(Membres),1);
figure; UnivarScatter(Conc_VIPMembres, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
xlabel('Xsens positions')
ylabel('VIP')
axis([0 16 0.5 4])

%% Third iteration with VIP thresh varying of 0.1
% We want to have less error than original model and less relevant variables
VAREXPLAINED_IT = NaN(50,20);
MABSERROR_IT = NaN(50,20);
NRelevant = NaN(50,25);
VarName_Selected = {};
Relevant = {};
ABSE = ABSE_Original;
MVAREXPLAINED = NaN(50,20,20);
MABSERROR = NaN(50,20,20);
RPE_OBS = [];
RPE_PRED = [];
BETA_MODELS = NaN(1261,50,20);
R = [];
for iIT = 1:50
    Participants = [1:length(Part1)]';
    k = 5;
    c = cvpartition(length(Participants),'KFold',k);
    
    IT = ['IT' num2str(iIT)];
    VIP_Thresh = 1;
    inc = 1;
    VAREXPLAINED = NaN(k,20,20);
    ABSERROR = NaN(k,20,20);
    NRelevant(iIT,inc) = 1260;
    RPE_obs = NaN(1000,5,20);
    RPE_pred = NaN(1000,5,20);
    
    while min(NRelevant(iIT,~isnan(NRelevant(iIT,:)))) > 1 %ABSE <= ABSE_Original
        BETA = [];
        INDVIP = {};
        VIPSCORE = [];
        for ik = 1:k
            Yobs = [];
            idx = training(c,ik);
            Part_Training_Set = Participants(idx);
            Part_Test_Set = Participants(idx==0);
            
            X_Training_Set = Xnorm_allParticipants(ismember(Data_STD_Jerk(:,1),Part1(Part_Training_Set)),:);
            X_Test_Set = Xnorm_allParticipants(ismember(Data_STD_Jerk(:,1),Part1(Part_Test_Set)),:);
            
            Y_Training_Set = Ynorm_allParticipants(ismember(Data_STD_Jerk(:,1),Part1(Part_Training_Set)),:);
            Y_Test_Set = Ynorm_allParticipants(ismember(Data_STD_Jerk(:,1),Part1(Part_Test_Set)),:);
            
            Row_Part_Training_Set = Data_STD_Jerk(ismember(Data_STD_Jerk(:,1),Part1(Part_Training_Set)),1);
            Row_Part_Test_Set = Data_STD_Jerk(ismember(Data_STD_Jerk(:,1),Part1(Part_Test_Set)),1);
            
            %% Apply PLSR on training set
            ncomp = 20; % use only 10 components ? change this to desired number of components
            [x_loadings, y_loadings, x_scores, y_scores, beta, pctvar, mse, stats] = plsregress(X_Training_Set, Y_Training_Set, ncomp);
            
            BETA = [BETA beta];
            
            %         Ynorm = (Y-mean(Y))/(std(Y)); % zscore
            for iP = 1:length(Part_Test_Set)
                Yobs = [Yobs; dat(dat(:,1)==Part1(Part_Test_Set(iP)),2)];
            end
            % Sub-class Yobs
%             Yobs(Yobs < 2) = 1;
%             Yobs(Yobs < 3 & Yobs > 1.5) = 2;
%             Yobs(Yobs < 5 & Yobs > 2.5) = 3;
%             Yobs(Yobs < 7 & Yobs > 4) = 4;
%             Yobs(Yobs > 6) = 5;
            
%             Yobs(Yobs < 3) = 1;
%             Yobs(Yobs < 5 & Yobs > 2.5) = 2;
%             Yobs(Yobs > 4) = 3;
            
            yfitnorm = [ones(size(X_Test_Set,1),1) X_Test_Set]*beta;
            %     yfitnorm = [ones(size(X_Test_Set,1),1) X_Test_Set(:,Relevant)]*beta;
            yfit = (yfitnorm*std(Yobs))+mean(Yobs);
            % Sub-class Yobs
%             yfit(yfit < 2) = 1;
%             yfit(yfit < 3 & yfit > 1.5) = 2;
%             yfit(yfit < 5 & yfit > 2.5) = 3;
%             yfit(yfit < 7 & yfit > 4) = 4;
%             yfit(yfit > 6) = 5;
            
%             yfit(yfit < 3) = 1;
%             yfit(yfit < 5 & yfit > 2.5) = 2;
%             yfit(yfit > 4) = 3;
            
            residuals = Yobs - yfit;
            
            % VIP
            W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
            p = size(x_loadings,1);
            sumSq = sum(x_scores.^2,1).*sum(y_loadings.^2,1);
            vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
            indVIP = find(vipScore > VIP_Thresh(inc));
            INDVIP{ik} = indVIP;
            VIPSCORE(:,ik) = vipScore;
        end
        [SORT_VIPSCORE, IDSORT] = sort(mean(VIPSCORE,2),'descend');
        SORT_Variables = VarNames(IDSORT);
        
        sum_INDVIP = zeros(size(X_Training_Set,2),1);
        for ivar = 1:size(X_Training_Set,2)
            for ik = 1:k
                y = sum(INDVIP{ik}==ivar);
                sum_INDVIP(ivar,1) = sum_INDVIP(ivar,1) + y;
            end
        end
        Relevant{iIT,inc} = find(sum_INDVIP > 1);
        NRelevant(iIT,inc) = length(Relevant{iIT,inc});
        if NRelevant(iIT,inc) == 0
            break
        end
        
        for ncomp = 1:20
            %         ABSERROR = [];
            %         VAREXPLAINED = [];
            %         INDVIP = {};
            %         VIPSCORE = [];
            
            BETA = [];
            for ik = 1:k
                Yobs = [];
                idx = training(c,ik);
                Part_Training_Set = Participants(idx);
                Part_Test_Set = Participants(idx==0);
                
                X_Training_Set = Xnorm_allParticipants(ismember(Data_STD_Jerk(:,1),Part1(Part_Training_Set)),:);
                X_Test_Set = Xnorm_allParticipants(ismember(Data_STD_Jerk(:,1),Part1(Part_Test_Set)),:);
                
                Y_Training_Set = Ynorm_allParticipants(ismember(Data_STD_Jerk(:,1),Part1(Part_Training_Set)),:);
                Y_Test_Set = Ynorm_allParticipants(ismember(Data_STD_Jerk(:,1),Part1(Part_Test_Set)),:);
                
                Row_Part_Training_Set = Data_STD_Jerk(ismember(Data_STD_Jerk(:,1),Part1(Part_Training_Set)),1);
                Row_Part_Test_Set = Data_STD_Jerk(ismember(Data_STD_Jerk(:,1),Part1(Part_Test_Set)),1);
                
                %% Apply PLSR on training set
                if NRelevant(iIT,inc) < ncomp
                    ncomp = NRelevant(iIT,inc);
                end
                    
%                 if NRelevant(iIT,inc) < 4
%                     ncomp = NRelevant(iIT,inc);
%                 else
%                     ncomp = 4; % use only 10 components ? change this to desired number of components
%                 end
                %     [x_loadings, y_loadings, x_scores, y_scores, beta, pctvar, mse, stats] = plsregress(X_Training_Set, Y_Training_Set, ncomp);
                [x_loadings, y_loadings, x_scores, y_scores, beta, pctvar, mse, stats] = plsregress(X_Training_Set(:,Relevant{iIT,inc}), Y_Training_Set, ncomp);
                CumVar = cumsum(100*pctvar(2,:));
                VAREXPLAINED(ik,ncomp,inc) = CumVar(:,end);
                
                BETA = [BETA beta];
                
                %         Ynorm = (Y-mean(Y))/(std(Y)); % zscore
                for iP = 1:length(Part_Test_Set)
                    Yobs = [Yobs; dat(dat(:,1)==Part1(Part_Test_Set(iP)),2)];
                end
                % Sub-class Yobs
%                 Yobs(Yobs < 2) = 1;
%                 Yobs(Yobs < 3 & Yobs > 1.5) = 2;
%                 Yobs(Yobs < 5 & Yobs > 2.5) = 3;
%                 Yobs(Yobs < 7 & Yobs > 4) = 4;
%                 Yobs(Yobs > 6) = 5;
                
%                 Yobs(Yobs < 3) = 1;
%                 Yobs(Yobs < 5 & Yobs > 2.5) = 2;
%                 Yobs(Yobs > 4) = 3;
                
                %     yfitnorm = [ones(size(X_Test_Set,1),1) X_Test_Set]*beta;
                yfitnorm = [ones(size(X_Test_Set,1),1) X_Test_Set(:,Relevant{iIT,inc})]*beta;
                yfit = (yfitnorm*std(Yobs))+mean(Yobs);
                % Sub-class Yobs
%                 yfit(yfit < 2) = 1;
%                 yfit(yfit < 3 & yfit > 1.5) = 2;
%                 yfit(yfit < 5 & yfit > 2.5) = 3;
%                 yfit(yfit < 7 & yfit > 4) = 4;
%                 yfit(yfit > 6) = 5;
                
%                 yfit(yfit < 3) = 1;
%                 yfit(yfit < 5 & yfit > 2.5) = 2;
%                 yfit(yfit > 4) = 3;
                
                residuals = Yobs - yfit;
                if ncomp == 4
                    RPE_obs(1:length(Yobs),ik,inc) = Yobs;
                    RPE_pred(1:length(Yobs),ik,inc) = yfit;
                end
                                
                %             squaredResiduals = residuals.^2;
                %             Mean_squaredResiduals = mean(squaredResiduals);
                %             rmse = sqrt(Mean_squaredResiduals);
                %             RMSE = [RMSE; rmse];
                %
                %             press = sum(squaredResiduals);
                %             PRESS = [PRESS; press];
                
                abserror = mean(abs(Yobs - yfit));
                ABSERROR(ik,ncomp,inc) = abserror;
                
            end
            VarName_Selected{iIT,inc} = VarNames(Relevant{iIT,end});
            if ncomp == 4
                BETA_MODELS(1:length(BETA),iIT,inc) = mean(BETA,2);
            end
        end
        MVAREXPLAINED(iIT,:,inc) = mean(VAREXPLAINED(:,:,inc));
        MABSERROR(iIT,:,inc) = mean(ABSERROR(:,:,inc));
        
        inc = inc +1;
        VIP_Thresh(inc) = VIP_Thresh(inc-1) + 0.1;
    end

    for ipage = 1:20
        RPE_OBS_temp = RPE_obs(:,:,ipage);
        RPE_OBS_temp = RPE_OBS_temp(~isnan(RPE_OBS_temp));
        RPE_PRED_temp = RPE_pred(:,:,ipage);
        RPE_PRED_temp = RPE_PRED_temp(~isnan(RPE_PRED_temp));
        if length(RPE_OBS_temp) > 0
            R(iIT,ipage) = corr(RPE_OBS_temp,RPE_PRED_temp,'Type','Pearson');
        end
        if Task == 'Ha' & ipage == 3 | Task == 'Li' & ipage == 4
            RPE_OBS = [RPE_OBS ; RPE_OBS_temp];
            RPE_PRED = [RPE_PRED ; RPE_PRED_temp];
        end
    end
    
%     RPE_OBS_temp = RPE_obs(~isnan(RPE_obs(:,:,4)));
%     RPE_OBS = [RPE_OBS ; RPE_OBS_temp];
%     RPE_OBS = RPE_OBS(~isnan(RPE_OBS));
%     RPE_PRED_temp = RPE_pred(~isnan(RPE_pred(:,:,4)));
%     RPE_PRED = [RPE_PRED ; RPE_PRED_temp];
%     RPE_PRED = RPE_PRED(~isnan(RPE_PRED));

%     R(iIT) = corr(RPE_OBS_temp(~isnan(RPE_OBS_temp)),RPE_PRED_temp(~isnan(RPE_PRED_temp)),'Type','Pearson');
end
% cd(['J:\Piano_Fatigue\Data_Exported'])
% save(['BETA_MODELS_' Task '.mat'],['BETA_MODELS'])
% save(['RELEVANT_MODELS_' Task '.mat'],['Relevant'])
% save(['SORT_Variables_' Task '.mat'],['SORT_Variables'])
% sum(RPE_OBS_temp < 2)
% sum(RPE_OBS_temp < 3 & RPE_OBS_temp > 1.5)
% sum(RPE_OBS_temp < 5 & RPE_OBS_temp > 2.5)
% sum(RPE_OBS_temp < 7 & RPE_OBS_temp > 4)
% sum(RPE_OBS_temp > 6)

%%
if Task == 'Ha'
    Colors = repmat([0 0.4470 0.7410],20,1);
else
    Colors = repmat([0.8500 0.3250 0.0980],100,1);
end
TableResults = nanmean(NRelevant)';
TableResults(:,2) = std(NRelevant)';
TableResults = TableResults(~isnan(TableResults(:,1)),:);
for i = 1:length(TableResults)
    TableResults(i,3) = find(mean(MABSERROR(:,:,i)) == min(mean(MABSERROR(:,:,i))));
    TableResults(i,4) = mean(MVAREXPLAINED(:,TableResults(i,3),i));
    TableResults(i,5) = std(MVAREXPLAINED(:,TableResults(i,3),i));
    TableResults(i,6) = mean(MABSERROR(:,TableResults(i,3),i));
    TableResults(i,7) = std(MABSERROR(:,TableResults(i,3),i));
    
%     figure(1); UnivarScatter(MVAREXPLAINED(:,:,i), 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
%     xlabel('Number of latent variables')
%     ylabel('RPE explained (%)')
%     axis([0 13 0 100])
%     title(['Model reduced ' num2str(i) ' NRelevant = ' num2str(mean(NRelevant(:,i)))])
%     
%     figure(2); UnivarScatter(MABSERROR(:,:,i), 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
%     xlabel('Number of latent variables')
%     ylabel('Absolute error')
%     axis([0 13 0.8 2.4])
%     title(['Model reduced ' num2str(i) ' NRelevant = ' num2str(mean(NRelevant(:,i)))])
%     
%     pause
end

%% Stat
V1 = fullModel_ABSERROR_IT(:,2);
V2 = fullModel_VAREXPLAINED_IT(:,2);
for i = 1:length(TableResults)-1
    V1(:,i+1) = MABSERROR(:,TableResults(i,3),i);
    V2(:,i+1) = MVAREXPLAINED(:,TableResults(i,3),i);
end
XVLat = find(mean(V1) == min(mean(V1)));
if Task == 'Ha'
    Colors = repmat([0 0.4470 0.7410],20,1);
else
    Colors = repmat([0.8500 0.3250 0.0980],100,1);
end
[p,tbl,stats] = anova1(V2);
multcomp = multcompare(stats);
multcomp = multcomp(multcomp(:,2) == XVLat | multcomp(:,1) == XVLat,:);
figure(15); UnivarScatter(V2, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
xlabel('Model')
ylabel('RPE explained (%)')
axis([0 size(V2,2)+1 0 100])
title([Task ' Task'])
inc = 1;
for iTable = 1:size(V2,2)-1
    if inc == XVLat
        inc = inc + 1;
    end
    if multcomp(iTable,6) < 0.05
        hold on;
        plot(inc, 90, '*', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
    end
    inc = inc + 1;
end

[p,tbl,stats] = anova1(V1);
multcomp = multcompare(stats);
multcomp = multcomp(multcomp(:,2) == XVLat | multcomp(:,1) == XVLat,:);
figure(16); UnivarScatter(V1, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
xlabel('Model')
ylabel('Absolute error')
axis([0 size(V2,2)+1 0.8 2.4])
title([Task ' Task'])
inc = 1;
for iTable = 1:size(V1,2)-1
    if inc == XVLat
        inc = inc + 1;
    end
    if multcomp(iTable,6) < 0.05
        hold on;
        plot(inc, 1.6, '*', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
    end
    inc = inc + 1;
end

%% Correlation
C = confusionmat(RPE_OBS,round(RPE_PRED));
confusionchart(C)


for i = 1:15
    subplot(3,5,i);
    plot(R(:,i))
end

XVLat = 4;
RPE_OBS = RPE_obs(:,:,XVLat);
RPE_OBS = reshape(RPE_OBS,[length(RPE_OBS)*5 1]);
RPE_PRED = RPE_pred(:,:,XVLat);
RPE_PRED = reshape(RPE_PRED,[length(RPE_PRED)*5 1]);

figure; plot(RPE_OBS,RPE_PRED,'o')
tbl = table(RPE_OBS,RPE_PRED);
mdl = fitlm(tbl,'RPE_PRED ~ RPE_OBS');
plot(mdl)
legend(['R² = ' num2str(mdl.Rsquared.Adjusted)])

[R,P] = corrcoef(RPE_OBS(~isnan(RPE_OBS)),RPE_PRED(~isnan(RPE_PRED)));
% [R,P] = 
corr(RPE_OBS(~isnan(RPE_OBS)),RPE_PRED(~isnan(RPE_PRED)),'Type','Pearson')

figure; plot(RPE_OBS(~isnan(RPE_OBS)))
hold on
plot(RPE_PRED(~isnan(RPE_PRED)))

%% Data Circos
Mean_VIPSCORE_IT = mean(VIPSCORE_IT,2);
[SORT_VIPSCORE, IDSORT] = sort(Mean_VIPSCORE_IT,'descend');
SORT_Variables = VarNames(IDSORT);
a = sort(IDSORT(1:275),'ascend');
SORT_Variables = VarNames(a);

VIPSCORE_CIRCOS_REVIEW = Mean_VIPSCORE_IT(a,1);
Varnames_REVIEW = VarNames_dec(:,a);
Varnames_CIRCOS_REVIEW = {};
for i = 1:size(Varnames_REVIEW,2)
    Varnames_CIRCOS_REVIEW{i,1} = Varnames_REVIEW{1,i};
    Varnames_CIRCOS_REVIEW{i,2} = [Varnames_REVIEW{2,i} ' ' Varnames_REVIEW{4,i}];
    Varnames_CIRCOS_REVIEW{i,3} = Varnames_REVIEW{3,i};
end
cd(['J:\Piano_Fatigue\Data_Exported'])
save(['VIP_CIRCOS_Review_' Task '.mat'],['VIPSCORE_CIRCOS_REVIEW'])
save(['VarNames_CIRCOS_Review_' Task '.mat'],['Varnames_CIRCOS_REVIEW'])

%% Variation best variables
if Task == 'Ha'
    Colors = 'b';
else
    Colors = 'r';
end
TableEvolution = [];
TableEvolution_rowNames = {};
TableEvolution_colNames = {'Pelvis_Ev' 'Pelvis_p' 'Trunk_Ev' 'Trunk_p' 'Head_Ev' 'Head_p'...
    'Shoulder_Ev' 'Shoulder_p' 'Arm_Ev' 'Arm_p' 'Forearm_Ev' 'Forearm_p' 'Hand_Ev' 'Hand_p'};
inc = 1;
Feat = 'MedianFreq';
idFeat = 1;
for iSelVar = 1:10%length(a)
%     [Moy,CI,STD] = grpstats(Xnorm_allParticipants_Timenorm(:,a(iSelVar)),repmat([1:100],1,length(Part1)),{'mean','meanci','std'});

    [Moy,CI,STD] = grpstats(Xnorm_allParticipants_Timenorm(:,IDSORT(iSelVar)),repmat([1:100],1,length(Part1)),{'mean','meanci','std'});
    figure(1); subplot(2,5,iSelVar); plot_ci(1:length(Moy),[Moy,Moy-STD,Moy+STD],'PatchColor', Colors, 'PatchAlpha', 0.3, 'MainLineColor', Colors,'LineWidth', 0.1, 'LineColor', 'w')
    title([SORT_Variables{iSelVar}])
    xlabel('Time (%)')
    ylabel('Normalized amplitude')
    xticks([0:20:100])

    MeanVar = mean(Moy(end-9:end))-mean(Moy(1:10));
    [h,p,ci,stats] = ttest(Moy(1:10),Moy(end-9:end));
    
    newFeat = Varnames_CIRCOS_REVIEW{iSelVar,1};
    if ~strcmp(newFeat,Feat)
        Feat = newFeat;
        idFeat = iSelVar;
    end

    if iSelVar > idFeat & Varnames_CIRCOS_REVIEW{iSelVar,1} == Feat
        if any(strcmp({Varnames_CIRCOS_REVIEW{idFeat:iSelVar-1,2}},Varnames_CIRCOS_REVIEW{iSelVar,2}))
            idSame = find(strcmp({Varnames_CIRCOS_REVIEW{idFeat:iSelVar-1,2}},Varnames_CIRCOS_REVIEW{iSelVar,2}));
            Varnames_CIRCOS_REVIEW{iSelVar,4} = Varnames_CIRCOS_REVIEW{idFeat+idSame(1)-1,4};
            inc = Varnames_CIRCOS_REVIEW{iSelVar,4};
            TableEvolution_rowNames{inc} = [Varnames_CIRCOS_REVIEW{iSelVar,1} ' ' Varnames_CIRCOS_REVIEW{iSelVar,2}];
        else
            Varnames_CIRCOS_REVIEW{iSelVar,4} = inc;
            TableEvolution_rowNames{inc} = [Varnames_CIRCOS_REVIEW{iSelVar,1} ' ' Varnames_CIRCOS_REVIEW{iSelVar,2}];
        end
    else
        Varnames_CIRCOS_REVIEW{iSelVar,4} = inc;
        TableEvolution_rowNames{inc} = [Varnames_CIRCOS_REVIEW{iSelVar,1} ' ' Varnames_CIRCOS_REVIEW{iSelVar,2}];
    end
    
    if strcmp(Varnames_CIRCOS_REVIEW{iSelVar,3}, 'L5')
        pos = 1;
    elseif strcmp(Varnames_CIRCOS_REVIEW{iSelVar,3}, 'T8')
        pos = 3;
    elseif strcmp(Varnames_CIRCOS_REVIEW{iSelVar,3}, 'Head')
        pos = 5;
    elseif strcmp(Varnames_CIRCOS_REVIEW{iSelVar,3}, 'Shoulder')
        pos = 7;
    elseif strcmp(Varnames_CIRCOS_REVIEW{iSelVar,3}, 'Arm')
        pos = 9;
    elseif strcmp(Varnames_CIRCOS_REVIEW{iSelVar,3}, 'Forearm')
        pos = 11;
    elseif strcmp(Varnames_CIRCOS_REVIEW{iSelVar,3}, 'Hand')
        pos = 13;
    end
    TableEvolution(inc,pos) = MeanVar;
    TableEvolution(inc,pos+1) = p;
    inc = size(TableEvolution,1)+1;
end
TableEvolution = array2table(TableEvolution,'RowNames', TableEvolution_rowNames,'VariableNames', TableEvolution_colNames);
cd(['J:\Piano_Fatigue\Data_Exported'])
writetable(TableEvolution,['TableEvolution_' Task '.xlsx'],'WriteRowNames',true)  

for iL = 1:100:length(Xnorm_allParticipants_Timenorm)
    plot(Xnorm_allParticipants_Timenorm(iL:iL+99,IDSORT(iSelVar)))
    pause
end

[p,tbl,stats] = anova1([fullModel_ABSERROR_IT MABSERROR_IT]);
multcomp = multcompare(stats);
multcomp = multcomp(multcomp(:,1)==1,:);
if Task == 'Ha'
    Colors = repmat([0 0.4470 0.7410],length(ABSERROR)+1,1);
else
    Colors = repmat([0.8500 0.3250 0.0980],length(ABSERROR)+1,1);
end
figure; UnivarScatter([fullModel_ABSERROR_IT MABSERROR_IT(:,1:end)], 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
for iMod = 1:length(multcomp)
    if multcomp(iMod,6) < 0.05
        hold on;
        plot(iMod+1, 1.6, '*', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
    end
end
xlabel('VIP threshold')
ylabel('Absolute error')
xticklabels([1:0.1:2.7])
xlim([0 length(multcomp)+2])
ylim([0.6 2])

if Task == 'Ha'
    Colors = repmat([0 0.4470 0.7410],length(ABSERROR)+1,1);
else
    Colors = repmat([0.8500 0.3250 0.0980],length(ABSERROR)+1,1);
end
figure; UnivarScatter([repmat(1260,length(NRelevant),1) NRelevant], 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
set(gca, 'YScale', 'log')
xlabel('VIP threshold')
ylabel('Number of variables retained')
yticks([0:1:2:10,20:10:50,100:50:400,500:100:1300])
xticklabels([1:0.1:2.7])
xlim([0 length(multcomp)+2])
ylim([0.1 1260])

[p,tbl,stats] = anova1([fullModel_VAREXPLAINED_IT(:,end) VAREXPLAINED_IT]);
multcomp = multcompare(stats);
multcomp = multcomp(multcomp(:,1)==1,:);
if Task == 'Ha'
    Colors = repmat([0 0.4470 0.7410],length(ABSERROR)+1,1);
else
    Colors = repmat([0.8500 0.3250 0.0980],length(ABSERROR)+1,1);
end
figure; UnivarScatter([fullModel_VAREXPLAINED_IT(:,end) VAREXPLAINED_IT(:,1:end)], 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
for iMod = 1:length(multcomp)
    if multcomp(iMod,6) < 0.05
        hold on;
        plot(iMod+1, 90, '*', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
    end
end
xlabel('VIP threshold')
ylabel('RPE explained (%)')
xticklabels([1:0.1:2.7])
xlim([0 length(multcomp)+2])
ylim([0 90])

% yticks([0:1:10,20:10:100])
% axis([0 23 0 100])

%%
% Mean value of Peak power freq
PPF_initiation = [];
PPF_termination = [];
for iP = 1:49
    X = [];
    Xnorm = [];
    Labels = [];
%     if Task == 'Ha' & iP == 22
%     elseif Task == 'Li' & iP == 17
%     else
        for iFeature = 7
            Feature = Features{iFeature};
            for iXSENS = 1:length(XSENSs)
                XSENS = XSENSs{iXSENS};
                dat = eval(['Data_' Feature '_' XSENS]);
%                 for iMembres = 1:length(Membres)
%                     Lab{iMembres} = [Feature '_' XSENS '_' Membres{iMembres}];
%                 end
%                 Labels = [Labels Lab];

                temp = dat(dat(:,1)==iP,:);
                X = [X temp(:,3:end)];
            end
        end
        if length(X) > 0
        Y = temp(:,2);
        
        % Normalize X and Y
        for iC = 1:size(X,2)
            Xnorm(:,iC) = (X(:,iC)-mean(X(:,iC)))/(std(X(:,iC))); % zcore
        end
        Xnorm(:,find(isnan(Xnorm(1,:))))=0;
        Ynorm = (Y-mean(Y))/(std(Y)); % zscore
%         if sum(isnan(Ynorm))>=1
%             break
%         end
                
        Xnorm_allParticipants = [Xnorm_allParticipants; Xnorm];
        Ynorm_allParticipants = [Ynorm_allParticipants; Ynorm];
        end
%     end
end

%% Cross models
clear
close all
clc
cd(['J:\Piano_Fatigue\Data_Exported'])
load('Xnorm_AllParticipants_Ha.mat')
load('Ynorm_AllParticipants_Ha.mat')
load('Xnorm_AllParticipants_Li.mat')
load('Ynorm_AllParticipants_Li.mat')
BETA_MODELS_Ha = load('BETA_MODELS_Ha.mat');
BETA_Ha = BETA_MODELS_Ha.BETA_MODELS(:,:,3);
BETA_MODELS_Li = load('BETA_MODELS_Li.mat');
BETA_Li = BETA_MODELS_Li.BETA_MODELS(:,:,4);
load('Yobs_allHa.mat')
load('Yobs_allLi.mat')
Relevant_Ha = load('RELEVANT_MODELS_Ha.mat');
Relevant_Li = load('RELEVANT_MODELS_Li.mat'); 

for iIT = 1:50
yfitnorm_Ha = [ones(size(Xnorm_allParticipants_Ha,1),1) ...
    Xnorm_allParticipants_Ha(:,Relevant_Li.Relevant{iIT, 4})]*...
    BETA_MODELS_Li.BETA_MODELS(~isnan(BETA_MODELS_Li.BETA_MODELS(:,iIT,4)),iIT,4);
yfitnorm_Li = [ones(size(Xnorm_allParticipants_Li,1),1) ...
    Xnorm_allParticipants_Li(:,Relevant_Ha.Relevant{iIT, 3})]*...
    BETA_MODELS_Ha.BETA_MODELS(~isnan(BETA_MODELS_Ha.BETA_MODELS(:,iIT,3)),iIT,3);
yfit_Ha = (yfitnorm_Ha*std(Yobs_allHa))+mean(Yobs_allHa);
yfit_Li = (yfitnorm_Li*std(Yobs_allLi))+mean(Yobs_allLi);
residuals_Ha = Yobs_allHa - yfit_Ha;
residuals_Li = Yobs_allLi - yfit_Li;

% squaredResiduals_Ha = residuals_Ha.^2;
% squaredResiduals_Li = residuals_Li.^2;
% Mean_squaredResiduals_Ha = mean(squaredResiduals_Ha);
% Mean_squaredResiduals_Li = mean(squaredResiduals_Li);
% RMSE_Ha = sqrt(Mean_squaredResiduals_Ha);
% RMSE_Li = sqrt(Mean_squaredResiduals_Li);
    
ABSERROR_Ha(iIT) = mean(abs(Yobs_allHa - yfit_Ha));
ABSERROR_Li(iIT) = mean(abs(Yobs_allLi - yfit_Li));

% plot(yfit_Ha)
% hold on
% plot(Yobs_allHa)
% 
% plot(yfit_Li)
% hold on
% plot(Yobs_allLi)
end
%     % Leave-one-out cross validation
%     AllPart = Row_Part_Training_Set;
%     RMSE = [];
%     PRESS = [];
%     ABSERROR = [];
%     BETA = [];
%     for iP = 1:length(Part_Training_Set)
%             X_Training_Set_CrossVal = X_Training_Set;
%             Y_Training_Set_CrossVal = Y_Training_Set;
% 
%             % Discard paticipant
%             Row2Remove = find(AllPart==Part_Training_Set(iP));
%             X_Test_Set_CrossVal = X_Training_Set(Row2Remove,:);
%             Y_Test_Set_CrossVal = Y_Training_Set(Row2Remove,:);
%         
%             X_Training_Set_CrossVal(Row2Remove,:) = [];
%             Y_Training_Set_CrossVal(Row2Remove,:) = [];
% 
%             [x_loadings, y_loadings, x_scores, y_scores, beta, pctvar, mse, stats] = plsregress(X_Training_Set_CrossVal, Y_Training_Set_CrossVal, ncomp);
% %         plot(1:ncomp,cumsum(100*pctvar(2,:)),'-bo');
% %         xlabel('Components')
% %         ylabel('Percent Variance Explained in Y')
% 
%             BETA(:,iP) = beta;
% 
% %         Ynorm = (Y-mean(Y))/(std(Y)); % zscore
%             Ytest = dat(dat(:,1)==Part_Training_Set(iP),2);
% 
%             yfitnorm = [ones(size(X_Test_Set_CrossVal,1),1) X_Test_Set_CrossVal]*beta;
%             yfit = (yfitnorm*std(Ytest))+mean(Ytest);
%             residuals = Ytest - yfit;
% %         plot(yfitnorm)
% %         hold on
% %         plot(Y_Test_Set_CrossVal)
%         
%             figure(1); subplot(5,8,iP); plot(yfit)
%             hold on
%             plot(Ytest)
% %         title([num2str(round(Pearson(iP),2)) '-' num2str(round(Spearman(iP),2))])
%         
%             % RMSE
%             squaredResiduals = residuals.^2;
%             Mean_squaredResiduals = mean(squaredResiduals);
%             rmse = sqrt(Mean_squaredResiduals);
%             RMSE = [RMSE; rmse];
%         
%             press = sum(squaredResiduals);
%             PRESS = [PRESS; press];
%         
%             abserror = mean(abs(Ytest - yfit));
%             ABSERROR = [ABSERROR; abserror];
% %         relerror = mean((abs(Ytest - yfit)./Ytest)*100);
%     end


% Test the model on the Test Set
for iP = 1:49
    plot(dat(dat(:,1)==iP,2))
    hold on
    pause()
end

yfitTestSet = [ones(size(X_Test_Set,1),1) X_Test_Set]*beta;
Yfit = [];
Ytest = [];
for iTS = 1:length(Part_Test_Set)
    tempfit = yfitTestSet(Row_Part_Test_Set==Part_Test_Set(iTS));
    temptest = dat(dat(:,1)==Part_Test_Set(iTS),2);
    trans = (tempfit*std(temptest))+mean(temptest);
    Yfit = [Yfit; trans];
    Ytest = [Ytest; dat(dat(:,1)==Part_Test_Set(iTS),2)];
end
yfit = (yfitTestSet*std(Ytest))+mean(Ytest);
residuals = Ytest - Yfit;
plot(Ytest)
hold on
plot(Yfit)

% RMSE
squaredResiduals = residuals.^2;
Mean_squaredResiduals = mean(squaredResiduals);
rmse = sqrt(Mean_squaredResiduals);

abserror = mean(abs(Ytest - Yfit));

%% PCA
PCA_AllParticipants = [];
Y_AllParticipants = [];
for iP = 1:50
    X = [];
    Labels = [];
    if Task == 'Ha2' & iP == 22
    elseif Task == 'Li2' & iP == 17
    else
        for iFeature = 1:length(Features)
            Feature = Features{iFeature};
            for iXSENS = 1:length(XSENSs)
                XSENS = XSENSs{iXSENS};
                dat = eval(['Data_' Feature '_' XSENS]);
                for iMembres = 1:length(Membres)
                    Lab{iMembres} = [Feature '_' XSENS '_' Membres{iMembres}];
                end
                Labels = [Labels Lab];
                temp = dat(dat(:,1)==iP,:);
                X = [X temp(:,3:end)];
            end
        end
        Y = temp(:,2);
        
        % Normalize X
        for iC = 1:size(X,2)
            Xnorm(:,iC) = (X(:,iC)-mean(X(:,iC)))/(std(X(:,iC)));
        end
        Xnorm(:,find(isnan(Xnorm(1,:))))=0;

        % Apply PCA
        nPCA = 50;
        [coeff,score,latent] = pca(Xnorm);
        PCA = score(:,1:nPCA);
        
        PCA_AllParticipants = [PCA_AllParticipants; PCA];
        Y_AllParticipants = [Y_AllParticipants; Y];

%         explained(1:5,:)
        
%         % Number of components required to explain at least 95% variability
%         sum_explained = 0;
%         idx = 0;
%         while sum_explained < 95
%             idx = idx + 1;
%             sum_explained = sum_explained + explained(idx);
%         end
%         idx
%         pause()
        
%         % Cluster RPE
%         E = evalclusters(Y ,'kmeans','silhouette', 'KList', [1:10]);
%         [idx c sumd] = kmeans(Y,8);
%         a = idx(1);
%         inc = 1;
%         for iCluster = 1:length(idx)
%             b = idx(iCluster);
%             if b == a
%                 IDX(iCluster,:) = inc;
%             else
%                 a = idx(iCluster);
%                 inc = inc + 1;
%                 IDX(iCluster,:) = inc;
%             end
%         end
%         
%         % Scatter Plot
%         gscatter(score(:,1),score(:,2),IDX)
%         pause()
%         scatter(Y,score(:,1))
%         pause()
    end
end

%% Random Forest
XSENS_PCA_All = array2table(PCA_AllParticipants);
XSENS_norm_All = array2table(Xnorm_allParticipants);
Conc_XSENS_PCA = [bestPred XSENS_PCA_All(:,1:10)];

% Train Bagged Ensemble of Regression Trees
t = templateTree('NumVariablesToSample','all','PredictorSelection','interaction-curvature','Surrogate','on');
rng(1); % For reproducibility
Mdl2 = fitrensemble(Conc_XSENS_PCA,Ynorm_allParticipants,'Method','Bag','NumLearningCycles',200, ...
    'Learners',t);

yHat = oobPredict(Mdl2);
R2 = corr(Mdl2.Y,yHat)^2

% Predictor Importance Estimation
impOOB = oobPermutedPredictorImportance(Mdl2);
figure
bar(impOOB)
title('Unbiased Predictor Importance Estimates')
xlabel('Predictor variable')
ylabel('Importance')
h = gca;
h.XTickLabel = Mdl2.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';

[impGain,predAssociation] = predictorImportance(Mdl2);
figure
subplot(2,1,1); plot(1:numel(Mdl2.PredictorNames),[impOOB'])
title('Predictor Importance Estimation Comparison')
subplot(2,1,2); plot(1:numel(Mdl2.PredictorNames),[impGain'])
xlabel('Predictor variable')
ylabel('Importance')
h = gca;
h.XTickLabel = Mdl2.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
legend('OOB permuted','MSE improvement')
grid on

figure
imagesc(predAssociation)
title('Predictor Association Estimates')
colorbar
h = gca;
h.XTickLabel = Mdl2.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
h.YTickLabel = Mdl2.PredictorNames;

% Define a threshold based on preceding graph to keep variables
bestPredid = find(impGain>0.25e-4);

bestPred = table(zeros(size(XSENS_norm_All,1),1));
bestPredNames = {};
VNames = {};
for ibestPredid = 1:length(bestPredid)
    if VarNames{bestPredid(ibestPredid)}(end) == 'X' & ismember(bestPred.Properties.VariableNames(end), ...
            XSENS_norm_All(:,bestPredid(ibestPredid):bestPredid(ibestPredid)+2).Properties.VariableNames)==0
        bestPred = [bestPred XSENS_norm_All(:,bestPredid(ibestPredid):bestPredid(ibestPredid)+2)];
        VNames = VarNames(bestPredid(ibestPredid):bestPredid(ibestPredid)+2);
        bestPredNames = [bestPredNames, VNames];
    end
    if VarNames{bestPredid(ibestPredid)}(end) == 'Y' & ismember(bestPred.Properties.VariableNames(end), ...
            XSENS_norm_All(:,bestPredid(ibestPredid)-1:bestPredid(ibestPredid)+1).Properties.VariableNames)==0
        bestPred = [bestPred XSENS_norm_All(:,bestPredid(ibestPredid)-1:bestPredid(ibestPredid)+1)];
        VNames = VarNames(bestPredid(ibestPredid)-1:bestPredid(ibestPredid)+1);
        bestPredNames = [bestPredNames, VNames];
    end
    if VarNames{bestPredid(ibestPredid)}(end) == 'Z' & ismember(bestPred.Properties.VariableNames(end), ...
            XSENS_norm_All(:,bestPredid(ibestPredid)-2:bestPredid(ibestPredid)).Properties.VariableNames)==0
        bestPred = [bestPred XSENS_norm_All(:,bestPredid(ibestPredid)-2:bestPredid(ibestPredid))];
        VNames = VarNames(bestPredid(ibestPredid)-2:bestPredid(ibestPredid));
        bestPredNames = [bestPredNames, VNames];
    end
    if VarNames{bestPredid(ibestPredid)}(end) == 'M'
        bestPred = [bestPred XSENS_norm_All(:,bestPredid(ibestPredid))];
        VNames = VarNames(bestPredid(ibestPredid));
        bestPredNames = [bestPredNames, VNames];
    end    
end
bestPred(:,1) = [];

% Grow Random Forest Using Reduced Predictor Set
t = templateTree('PredictorSelection','interaction-curvature','Surrogate','on');
MdlReduced = fitrensemble(bestPred,Ynorm_allParticipants,'Method','Bag', ...
    'NumLearningCycles',200,'Learners',t);

yHatReduced = oobPredict(MdlReduced);
r2Reduced = corr(MdlReduced.Y,yHatReduced)^2

% Predictor Importance Estimation
impOOB = oobPermutedPredictorImportance(MdlReduced);
figure
bar(impOOB)
title('Unbiased Predictor Importance Estimates')
xlabel('Predictor variable')
ylabel('Importance')
h = gca;
h.XTickLabel = MdlReduced.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
bestPredNames'

[impGain,predAssociation] = predictorImportance(MdlReduced);
figure
subplot(2,1,1); plot(1:numel(MdlReduced.PredictorNames),[impOOB'])
title('Predictor Importance Estimation Comparison')
subplot(2,1,2); plot(1:numel(MdlReduced.PredictorNames),[impGain'])
xlabel('Predictor variable')
ylabel('Importance')
h = gca;
h.XTickLabel = MdlReduced.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
legend('OOB permuted','MSE improvement')
grid on

figure
imagesc(predAssociation)
title('Predictor Association Estimates')
colorbar
h = gca;
h.XTickLabel = MdlReduced.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
h.YTickLabel = MdlReduced.PredictorNames;

%% Load normalized TFR
X_all = [];
Y_all = [];
for iP = 1:50
    X = [];
    if Task == 'Ha2' & iP == 22
    elseif Task == 'Li2' & iP ==17
    else
    if Task == 'Ha2'
        RPE = Info_participants(iP).Hanon;
    else
        RPE = Info_participants(iP).Liszt;
    end
    RPE = rmmissing(RPE);
    RPE = [0; RPE];
    normalize_RPE = interp1(1:length(RPE),RPE,linspace(1,length(RPE),1000))';
    Y = normalize_RPE;

    cd(['J:\Piano_Fatigue\Data_Exported\TFRnormalized_XSENS\' Task])
    load([FilenamesTFR{1,iP}])
    VarTFR = fieldnames(TFR);
    for iVar = 2:length(VarTFR)-1
        NTFR = fieldnames(TFR.(VarTFR{iVar}));
        for iTFR = 1:length(NTFR)
            temp = TFR.(VarTFR{iVar}).(NTFR{iTFR})';
            X = [X temp];
        end
    end
    X_all = [X_all; X];
    Y_all = [Y_all; Y];
    end
end
save(['X_all.mat'],'X_all')
save(['Y_all.mat'],'Y_all')



%% Anova 2 facteurs
for iFeature = 1:length(Features)
    Feature = Features{iFeature};
    for iXSENS = 1:length(XSENSs)
        XSENS = XSENSs{iXSENS};
        dat = eval(['Data_' Feature '_' XSENS]);
        dat_part1 = []; dat_part2 = [];
        dat_part1_I0 = []; dat_part2_I0= [];
        dat_part1_I1 = []; dat_part2_I1= [];
        dat_part1_I2 = []; dat_part2_I2= [];
        dat_part1_I3 = []; dat_part2_I3= [];
        dat_part1_I4 = []; dat_part2_I4= [];
        dat_part1_I5 = []; dat_part2_I5 = [];
        for iP = 1:length(part1)
            temp = dat(dat(:,1)==part1(iP),:);
            dat_part1 = [dat_part1; temp(:,:)];
            dat_part1_I0 = [dat_part1_I0; mean(temp(1:10,:))];
            dat_part1_I1 = [dat_part1_I1; mean(temp(10:20,:))];
            dat_part1_I2 = [dat_part1_I2; mean(temp(30:40,:))];
            dat_part1_I3 = [dat_part1_I3; mean(temp(50:60,:))];
            dat_part1_I4 = [dat_part1_I4; mean(temp(70:80,:))];
            dat_part1_I5 = [dat_part1_I5; mean(temp(90:100,:))];
        end
        for iP = 1:length(part2)
            temp = dat(dat(:,1)==part2(iP),:);
            dat_part2 = [dat_part2; temp(:,:)];
            dat_part2_I0 = [dat_part2_I0; mean(temp(1:10,:))];
            dat_part2_I1 = [dat_part2_I1; mean(temp(10:20,:))];
            dat_part2_I2 = [dat_part2_I2; mean(temp(30:40,:))];
            dat_part2_I3 = [dat_part2_I3; mean(temp(50:60,:))];
            dat_part2_I4 = [dat_part2_I4; mean(temp(70:80,:))];
            dat_part2_I5 = [dat_part2_I5; mean(temp(90:100,:))];
        end
                
        for iC = 3:6
            within_factor = [[dat_part1_I0(:,iC); dat_part2_I0(:,iC)] [dat_part1_I1(:,iC); dat_part2_I1(:,iC)]...
                [dat_part1_I2(:,iC); dat_part2_I2(:,iC)] [dat_part1_I3(:,iC); dat_part2_I3(:,iC)]...
                [dat_part1_I4(:,iC); dat_part2_I4(:,iC)] [dat_part1_I5(:,iC); dat_part2_I5(:,iC)]];
            between_factor = [repmat(1,length(part1),1); repmat(2,length(part2),1)];

            t = table([part1; part2],within_factor(:,1),within_factor(:,2),within_factor(:,3),within_factor(:,4),...
            within_factor(:,5),within_factor(:,6),between_factor,'VariableNames',{'Part','MF1','MF20','MF40',...
            'MF60','MF80','MF100','Group'});
            t.Group = categorical(t.Group);
            Meas = table([1 2 3 4 5 6]','VariableNames',{'Measurement'});
            rm = fitrm(t,'MF1-MF100~Group','WithinDesign',Meas);
            ranovatbl = ranova(rm)
            anovatbl = anova(rm)
            pause()
            
%             % comparaison de groupe à chaque instant
%             for iT = 1:6
%                 t = table(within_factor(:,iT),between_factor,'VariableNames',{'XSENS','Group'});
%                 t.Group = categorical(t.Group);
%                 [P,ANOVATAB,STATS] = anova1(within_factor(:,iT),between_factor)
%                 P_Group(:,iT) = P;
%                 [QGroup] = mafdr(P_Group','BHFDR',true);
% %                 ranovatbl = ranova(rm);
% %                 posthoc = multcompare(rm,'Measurement');
% %                 PostHocResults = posthoc([1 7 13 19 25],:);
% %                 [QTime] = mafdr(table2array(PostHocResults(:,5)),'BHFDR',true);
% %                 computeCohen_d(table2array(t(1:STATS.n(1),1)),table2array(t(STATS.n(2):end,1)),'independent')
%                 pause()
%             end
%             close all
            
            % Comparaison par groupe pour l'effet de fatigue
            for iG = 1:2
                if iG == 1
                    Colors = [0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410];
                    figure(1);subplot(2,2,iC-2);UnivarScatter(within_factor(1:length(part1),:), 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');  
                    title([Feature '__' XSENS '__' Membres{iC-2} '__G1']);
                    t = table([part1],within_factor(1:length(part1),1),within_factor(1:length(part1),2),...
                        within_factor(1:length(part1),3),within_factor(1:length(part1),4),within_factor(1:length(part1),5),...
                        within_factor(1:length(part1),6),between_factor(1:length(part1)),'VariableNames',{'Part','MF1',...
                        'MF20','MF40','MF60','MF80','MF100','Group'});
                else
                    Colors = [0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980];
                    figure(2);subplot(2,2,iC-2);UnivarScatter(within_factor(length(part1)+1:end,:), 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');  
                    title([Feature '_' XSENS '__' Membres{iC-2} '_G2']);
                    t = table([part2],within_factor(length(part1)+1:end,1),within_factor(length(part1)+1:end,2),...
                        within_factor(length(part1)+1:end,3),within_factor(length(part1)+1:end,4),within_factor(length(part1)+1:end,5),...
                        within_factor(length(part1)+1:end,6),between_factor(length(part1)+1:end),'VariableNames',{'Part','MF1',...
                        'MF20','MF40','MF60','MF80','MF100','Group'});
                end

                t.Group = categorical(t.Group);
                Meas = table([1 2 3 4 5 6]','VariableNames',{'Measurement'});
                rm = fitrm(t,'MF1-MF100~Part','WithinDesign',Meas);
                ranovatbl = ranova(rm);
                posthoc = multcompare(rm,'Measurement');
                PostHocResults = posthoc([1 7 13 19 25],:);
                PostHocResults_General = posthoc([5],:);
                [QTime] = mafdr(table2array(PostHocResults(:,5)),'BHFDR',true);
                P_General = table2array(PostHocResults_General(:,5));
                if P_General<0.05
                    hold on
                    plot(3.4,1.1,'k*')
                    plot(3.6,1.1,'k*')
                end
                for iQ = 1:length(QTime)
                    if QTime(iQ)<0.05
                        hold on
                        plot(iQ + 0.5,1.05,'k*')
                    end
                end
                computeCohen_d(table2array(t(1:end,6)),table2array(t(1:end,7)),'paired')
                pause()
            end
        end
        close all
    end
end


%% Autres
% Muscles = {...
%     'm1_2';'m2_3';'m3_4';'m4_5';'m5_6';'m6_7';'m8_9';'m9_10';'m10_11';'m11_12';'m12_13';'m13_14';'m15_16';'m16_17';'m17_18';...
%     'm18_19';'m19_20';'m20_21';'m22_23';'m23_24';'m24_25';'m25_26';'m26_27';'m27_28';'m29_30';'m30_31';'m31_32';'m32_33';...
%     'm33_34';'m34_35';'m36_37';'m37_38';'m38_39';'m39_40';'m40_41';'m41_42';'m43_44';'m44_45';'m45_46';'m46_47';'m47_48';...
%     'm48_49';'Biceps';'Triceps';'DeltAnt';'DeltMed';'SupTrap'};

% FoldersNames = dir(['J:\Piano_Fatigue\Data_Exported\Tables_XSENS_LMM\']);
% Filenames = {...
%     FoldersNames(3:end).name...
%    } ;
% save(['FileNames_Feature_XSENS.mat'],'Filenames')

% EMG median freq
%     cd(['J:\Piano_Fatigue\Data_Exported\TFR_MedianFreq_Li'])
%     load([FilenamesEMG{1,iP}])
%     VarInc = 1;
%     Sig = [];
%     for iSig = [42 41 40 36 35 34 30 29 28]
%         Sig(:,VarInc) = TFR_MedianFreq.MedianFreq.(Muscles{iSig});
%         VarInc = VarInc +1 ;
%     end
%     Cy = round(cycles(iP).seq*1024);
% %             plot(Sig(:,1))
% %             line([Cy(:,1),Cy(:,1)]', repmat(ylim,length(Cy(:,1)),1)','color','red')
% %             line([Cy(:,end),Cy(:,end)]', repmat(ylim,length(Cy(:,end)),1)','color','red')
% %             pause()
%     SIG = [];
%     for iCycles = 1:size(Cy,1)-1
%         SIG(iCycles,:) = mean(Sig(Cy(iCycles,1):Cy(iCycles,end),:));
%     end
%     SIG = mean(SIG,2);
%     normalize_SIG = interp1(1:length(SIG),SIG,linspace(1,length(SIG),100))';

% Regression

%         [MoyRPE,CIrpe,STDRPE] = grpstats(dat_part1(:,2),repmat([1:100],1,49),{'mean','meanci','std'});
%         plot_ci(1:length(MoyRPE),[MoyRPE,MoyRPE-STDRPE,MoyRPE+STDRPE],'PatchColor', 'b', 'PatchAlpha', 0.3, 'MainLineColor', 'b','LineWidth', 0.1, 'LineColor', 'w')
%         for iC = 3:size(dat,2)
%             [Moy,CI,STD] = grpstats(dat_part1(:,iC),repmat([1:100],1,49),{'mean','meanci','std'});
% %             scatter(MoyRPE,Moy)
%             figure(1); subplot(2,2,iC-2); plot_ci(1:length(Moy),[Moy,Moy-STD,Moy+STD],'PatchColor', 'b', 'PatchAlpha', 0.3, 'MainLineColor', 'b','LineWidth', 0.1, 'LineColor', 'w')
%             title(['Mean ' Feature ' ' XSENS])
%             mdl = fitlm(Moy,MoyRPE);
%             legend(['R² = ' num2str(mdl.Rsquared.Adjusted)])
%             figure(2); subplot(2,2,iC-2); plot_ci(1:length(Moy),[STD],'PatchColor', 'b', 'PatchAlpha', 0.3, 'MainLineColor', 'b','LineWidth', 0.1, 'LineColor', 'w')
%             title(['STD ' Feature ' ' XSENS])
%             mdl = fitlm(STD,MoyRPE);
%             legend(['R² = ' num2str(mdl.Rsquared.Adjusted)])
%             figure(3); subplot(2,2,iC-2); UnivarScatter([dat_part1_Debut(:,iC) dat_part1_Fin(:,iC)], 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black')  
%             title(['Mean 10 cycles ' Feature ' ' XSENS])
%             [h,p,ci,stats]=ttest(dat_part1_Debut(:,iC),dat_part1_Fin(:,iC),'Tail','both','Alpha',0.05);
%             legend(['p = ' num2str(p)])
%         end
% %         DAT = [DAT, dat(:,3:end)];
% %         mdl = fitlm(dat_part1(:,3:end),dat_part1(:,2))
%         %mdl = stepwiselm(dat(:,3:end),dat(:,2),'constant','Upper','linear','PEnter',0.05)
%         disp([Feature '_' XSENS])
%         pause()
%         cd(['J:\Piano_Fatigue\Data_Exported\Tables_XSENS_LMM'])
%     save(['Data_' Feature '_' XSENS '.mat'],['Data_' Feature '_' XSENS])
%         save(['datpart1_' Feature '_' XSENS '.mat'],'dat_part1')

% mdl = fitlm(DAT(:,3:end),dat(:,2))
% mdl = stepwiselm(DAT(:,3:end),dat(:,2),'constant','Upper','linear','PEnter',0.05)
% 
% corr(dat_part1(:,5),dat_part1(:,2),'type','Spearman')
% scatter(dat_part1(:,2),dat_part1(:,5))
% tbl = table(dat_part1(:,2),dat_part1(:,3),dat_part1(:,4),dat_part1(:,5),dat_part1(:,6),'VariableNames',{'RPE','Shoulder','Arm','ForeArm','Hand'});
% lm = fitlm(tbl,'RPE~Arm+ForeArm+Hand')
% 
% figure;
% for iP = 1:length(part1)
%     temp = dat_part1(dat_part1(:,1)==part1(iP),:);
%     scatter(temp(:,2),temp(:,5))
% %     hold on
%     pause()
% end



