clear ;
close all ;
clc ;

Muscles = {'Dent1';'Dent2';'DeltA';'DeltM';'DeltP'; ...
    'Bi';'Tri';'TrapInf';'TrapSup';'TrapMed';};

% Défini un vecteur avec le nom de tous tes features sans '_'
Features_Names = {'MedianFreq';...
    'SpectralEntropy';'ApEntropy';'FuzzyEntropy';'MSEntropy';'SampleEntropy';...
    'FractalHiguchi';'FractalKatz';'FractalHurst';...
    'CorrDim'};...Rajoute les autres

% Défini le nom des colonnes
VarNames = {};
for iFeatures_Names = 1:length(Features_Names)
    RepFeatureName = repmat({Features_Names{iFeatures_Names}},1,10);
    for iM = 1:length(Muscles)
        RepFeatureName{iM} = [RepFeatureName{iM} '_' Muscles{iM}];
    end
    VarNames = {VarNames{:}, RepFeatureName{:}};
end
        
% Extrait les données normalisées de variation dans les tables
Part = [];
X = [];
Y = [];
Xnorm = [];
Ynorm = [];
for iP = 1:24
    if iP == 18
    else
        % Ça c'était pour extraire le RPE et ne pas charger les données brutes
        % à chaque fois (gain de temps)
        %     cd(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\EMG_Clean'])
        %     load(['P' num2str(iP) '_Pointage.mat'])
        %     Borg = EMG.Borg(:,2);
        %     cd(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\Borg_Scores'])
        %     save(['P' num2str(iP) '_Pointage_Borg.mat'],'Borg')
        % end
        % Crée un vecteur Part pour savoir à quel participant correspond chaque
        % ligne de X et Y (PLS)
        
        % Importe le RPE dans une matrice
        load(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\Borg_Scores\P' num2str(iP) '_Pointage_Borg.mat']);
        BorgZ = (Borg-mean(Borg))/(std(Borg));
        Y = [Y; Borg];
        Ynorm = [Ynorm; BorgZ];
        part = repmat(iP,length(Borg),1);
        Part = [Part; part];
        
        % Importe tes features EMG
        load(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\Results\TFR_MedianFreq\Instant_RPE\val_MedianFreq_P' num2str(iP) '.mat'])
        load(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\Results\TFR_SpectralEntropy\Instant_RPE\val_SpectralEntropy_P' num2str(iP) '.mat'])
        load(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\Results\EntropyMeasures\ApEntropy\Instant_RPE\val_ApEnt_P' num2str(iP) '.mat'])
        load(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\Results\EntropyMeasures\FuzzyEntropy\Instant_RPE\val_FuzzyEnt_P' num2str(iP) '.mat'])
        load(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\Results\EntropyMeasures\MSEntropy\Instant_RPE\val_MSEnt_P' num2str(iP) '.mat'])
        load(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\Results\EntropyMeasures\SampleEntropy\Instant_RPE\val_SamplEnt_P' num2str(iP) '.mat'])
        load(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\Results\Fractals\Higuchi\Instant_RPE\val_Hig_P' num2str(iP) '.mat'])
        load(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\Results\Fractals\Katz\Instant_RPE\val_Katz_P' num2str(iP) '.mat'])
        load(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\Results\Fractals\MFDFA\Instant_RPE\val_HurstExp_P' num2str(iP) '.mat'])
        load(['J:\IRSST_Fatigue\Pointage_repetitif_EMG\Results\Correlation\CorrelationDimension\Instant_RPE\val_CD_P' num2str(iP) '.mat'])
        ...Rajoute les autres

        val_Brute = {val_MedianFreq val_SpectralEntropy val_ApEnt val_FuzzyEnt val_MSEnt val_SamplEnt val_Hig...
                val_Katz val_HurstExp val_CD};
                ...Rajoute les autres

        % Calcule les valeurs Z
        val_Z = {};
        for iVal = 1:length(val_Brute)
            for iM = 1:length(Muscles)
                val_Z{iVal}(:,iM) = (val_Brute{iVal}(:,iM)-mean(val_Brute{iVal}(:,iM)))/(std(val_Brute{iVal}(:,iM)));
            end
        end
        
        % Concatene les features d'un participant dans le même ordre que
        % Features_Names
        temp = []; temp_norm = [];
        for iVal = 1:length(val_Brute)
            temp = [temp val_Brute{iVal}];
            temp_norm = [temp_norm val_Z{iVal}];
        end
    
        % Concatene les participants
        X = [X; temp];
        Xnorm = [Xnorm; temp_norm];
    end
end

%% Choix à faire sur l'analyse (tous les muscles ensemble ou muscle par muscle)
% Si tu analyses tous les muscles ensembles il faudra supprimer les sujets
% 9 et 13 car les NaN ne sont pas admis dans l'algorithme
X_PLS = X; Y_PLS = Y; Xnorm_PLS = Xnorm; Ynorm_PLS = Ynorm; Part_PLS = Part;
X_PLS(Part_PLS==9,:) = [];
Y_PLS(Part_PLS==9,:) = [];
Xnorm_PLS(Part_PLS==9,:) = [];
Ynorm_PLS(Part_PLS==9,:) = [];
Part_PLS(Part_PLS==9) = [];

X_PLS(Part_PLS==13,:) = [];
Y_PLS(Part_PLS==13,:) = [];
Xnorm_PLS(Part_PLS==13,:) = [];
Ynorm_PLS(Part_PLS==13,:) = [];
Part_PLS(Part_PLS==13) = [];
Participants = [1:8,10:12,14:17,19:24]';

% Si tu analyse muscle par muscle il faudra seulement supprimer les
% colonnes 6:8 du sujet 9
iM = 4;
X_PLS = X; Y_PLS = Y; Xnorm_PLS = Xnorm; Ynorm_PLS = Ynorm; Part_PLS = Part;
X_PLS = X_PLS(:,iM:10:size(X_PLS,2));
Xnorm_PLS = Xnorm_PLS(:,iM:10:size(Xnorm_PLS,2));
Participants = [1:17,19:24]';
if iM == 6 | iM == 7 | iM == 8
    X_PLS(Part==9,:) = [];
    Y_PLS(Part==9,:) = [];
    Xnorm_PLS(Part==9,:) = [];
    Ynorm_PLS(Part==9,:) = [];
    Part_PLS(Part==9) = [];
    Participants = [1:8,10:17,19:24]';
end

%% La première fois tu roules l'algorithme avec un max de ncomp
%% Define training set and test set
k = 5;
c = cvpartition(length(Participants),'KFold',k);
MABSERROR = [];
for ncomp = 1:size(Xnorm_PLS,2)
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
    Part_Training_Set = Participants(idx);
    Part_Test_Set = Participants(idx==0);
    
    X_Training_Set = Xnorm_PLS(ismember(Part_PLS,Part_Training_Set),:);
    X_Test_Set = Xnorm_PLS(ismember(Part_PLS,Part_Test_Set),:);
    
    Y_Training_Set = Ynorm_PLS(ismember(Part_PLS,Part_Training_Set),:);
    Y_Test_Set = Ynorm_PLS(ismember(Part_PLS,Part_Test_Set),:);
    
%     Row_Part_Training_Set = Part(ismember(Part,Part_Training_Set),1);
%     Row_Part_Test_Set = Part(ismember(Part,Part_Test_Set),1);
    
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
        Yobs = [Yobs; Y_PLS(Part_PLS==Part_Test_Set(iP))];
    end
    
    yfitnorm = [ones(size(X_Test_Set,1),1) X_Test_Set]*beta;
    %     yfitnorm = [ones(size(X_Test_Set,1),1) X_Test_Set(:,Releavant_Var)]*beta;
    yfit = (yfitnorm*std(Yobs))+mean(Yobs);
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
MABSERROR(ncomp) = mean(ABSERROR);
% mean(ABSERROR)
end

% La deuxième fois tu roules l'algorithme avec ncomp qui te donne la plus
% petite erreur de prédiction
ncomp = find(MABSERROR == min(MABSERROR));
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
    Part_Training_Set = Participants(idx);
    Part_Test_Set = Participants(idx==0);
    
    X_Training_Set = Xnorm_PLS(ismember(Part_PLS,Part_Training_Set),:);
    X_Test_Set = Xnorm_PLS(ismember(Part_PLS,Part_Test_Set),:);
    
    Y_Training_Set = Ynorm_PLS(ismember(Part_PLS,Part_Training_Set),:);
    Y_Test_Set = Ynorm_PLS(ismember(Part_PLS,Part_Test_Set),:);
    
%     Row_Part_Training_Set = Part(ismember(Part,Part_Training_Set),1);
%     Row_Part_Test_Set = Part(ismember(Part,Part_Test_Set),1);
    
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
        Yobs = [Yobs; Y_PLS(Part_PLS==Part_Test_Set(iP))];
    end
    
    yfitnorm = [ones(size(X_Test_Set,1),1) X_Test_Set]*beta;
    %     yfitnorm = [ones(size(X_Test_Set,1),1) X_Test_Set(:,Releavant_Var)]*beta;
    yfit = (yfitnorm*std(Yobs))+mean(Yobs);
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
mean(ABSERROR)

%% Quelques graph et stats
% Colors = [0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410];
Colors = [0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980];
UnivarScatter(VAREXPLAINED, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
xlabel('Components')
ylabel('% Var Exp in Y')
axis([0 ncomp+1 0 100])

% Colors = [0 0.4470 0.7410; 0 0.4470 0.7410];
Colors = [0.8500 0.3250 0.0980];
UnivarScatter([ABSERROR], 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
axis([0 2 0 2.6])

mean(VAREXPLAINED(:,end))
std(VAREXPLAINED(:,end))

mean(ABSERROR)
std(ABSERROR)
scatter(1:10,ABSERROR,'filled')
xlabel('Model')
ylabel('Absolute error')

%% Les variables les plus importantes
sum_INDVIP = zeros(size(X_Training_Set,2),1);
for ivar = 1:size(X_Training_Set,2)
    for ik = 1:5
        y = sum(INDVIP{ik}==ivar);
        sum_INDVIP(ivar,1) = sum_INDVIP(ivar,1) + y;
    end
end
Releavant_Var = find(sum_INDVIP > 3); % Le critère est arbritraire ! En gros sur 5 ittérations (k-folds avec k=5), combien de fois la variable a un VIP>1 
Mean_VIPSCORE = mean(VIPSCORE,2);
% Les 10 variables les plus importantes
% max10 = maxk(Mean_VIPSCORE,10);
% for iS = 1:10
%     idmax10(iS) = find(Mean_VIPSCORE == max10(iS));
% end
% VarNames(idmax10)

VarNames_Releavant_Var_dec = {};
VIPSCORE_Releavant_Var = VIPSCORE(Releavant_Var,:);
VarNames_Releavant_Var = VarNames(Releavant_Var);
for iV = 1:length(VarNames_Releavant_Var)
    temp = strsplit(VarNames_Releavant_Var{iV},'_');
    VarNames_Releavant_Var_dec = [VarNames_Releavant_Var_dec temp'];
end
% Features les plus importants
for iFeatures = 1:length(Features_Names)
    OccurenceFeatures(iFeatures) = sum(strcmp(VarNames_Releavant_Var_dec(1,:),Features_Names{iFeatures}));
    VIPFeatures{1,iFeatures} = mean(VIPSCORE_Releavant_Var(strcmp(VarNames_Releavant_Var_dec(1,:),Features_Names{iFeatures}),:),2);
end
Conc_VIPFeatures = padcat(VIPFeatures{1},VIPFeatures{2},VIPFeatures{3},VIPFeatures{4},VIPFeatures{5},VIPFeatures{6},VIPFeatures{7},VIPFeatures{8},VIPFeatures{9},VIPFeatures{10},VIPFeatures{11},VIPFeatures{12},VIPFeatures{13},VIPFeatures{14},VIPFeatures{15});
[p,tbl,stats] = anova1(Conc_VIPFeatures);
multcomp = multcompare(stats);
% Muscles les plus importants
for iM = 1:length(Muscles)
    OccurenceMembres(iM) = sum(count(VarNames_Releavant_Var_dec(2,:),Muscles{iM}));
    Conc_VIPMuscles{1,iM} = mean(VIPSCORE_Releavant_Var(strcmp(VarNames_Releavant_Var_dec(2,:),Muscles{iM}),:),2);
end
Conc_VIPMuscles = padcat(Conc_VIPMuscles{1},Conc_VIPMuscles{2},Conc_VIPMuscles{3},Conc_VIPMuscles{4},Conc_VIPMuscles{5},Conc_VIPMuscles{6},Conc_VIPMuscles{7},Conc_VIPMuscles{8},Conc_VIPMuscles{9},Conc_VIPMuscles{10});
[p,tbl,stats] = anova1(Conc_VIPMuscles);
multcomp = multcompare(stats);

Colors = repmat([0 0.4470 0.7410],length(Features_Names),1);
% Colors = repmat([0.8500 0.3250 0.0980],length(FeaturesVarNames),1);
figure; UnivarScatter(Conc_VIPFeatures, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
xlabel('Features')
ylabel('VIP')
axis([0 16 0.5 4])

Colors = repmat([0 0.4470 0.7410],length(Muscles),1);
% Colors = repmat([0.8500 0.3250 0.0980],length(XSENSsVarNames),1);
figure; UnivarScatter(Conc_VIPMuscles, 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
xlabel('Muscles')
ylabel('VIP')
axis([0 16 0.5 4])

%% La troisième fois tu roules l'algorithme avec seulement les variables importantes
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
    Part_Training_Set = Participants(idx);
    Part_Test_Set = Participants(idx==0);
    
    X_Training_Set = Xnorm_PLS(ismember(Part_PLS,Part_Training_Set),:);
    X_Test_Set = Xnorm_PLS(ismember(Part_PLS,Part_Test_Set),:);
    
    Y_Training_Set = Ynorm_PLS(ismember(Part_PLS,Part_Training_Set),:);
    Y_Test_Set = Ynorm_PLS(ismember(Part_PLS,Part_Test_Set),:);
    
%     Row_Part_Training_Set = Part(ismember(Part,Part_Training_Set),1);
%     Row_Part_Test_Set = Part(ismember(Part,Part_Test_Set),1);
    
    %% Apply PLSR on training set
%     ncomp = 2; % use only 10 components ? change this to desired number of components
%     [x_loadings, y_loadings, x_scores, y_scores, beta, pctvar, mse, stats] = plsregress(X_Training_Set, Y_Training_Set, ncomp);
        [x_loadings, y_loadings, x_scores, y_scores, beta, pctvar, mse, stats] = plsregress(X_Training_Set(:,Releavant_Var), Y_Training_Set, ncomp);
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
        Yobs = [Yobs; Y_PLS(Part_PLS==Part_Test_Set(iP))];
    end
    
%     yfitnorm = [ones(size(X_Test_Set,1),1) X_Test_Set]*beta;
    yfitnorm = [ones(size(X_Test_Set,1),1) X_Test_Set(:,Releavant_Var)]*beta;
    yfit = (yfitnorm*std(Yobs))+mean(Yobs);
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
mean(ABSERROR)


