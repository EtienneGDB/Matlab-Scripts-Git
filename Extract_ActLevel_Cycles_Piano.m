clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Entropy Fouaz'))

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
    FoldersNames(16:65).name...
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

Act_Cycles = [];
Act_Cycles_Interval = [];
for iM = 1:length(Muscles)
    Act_Cycles_Interval.(Muscles{iM}) = [0 0 0 0 0 0 0];
    Act_Cycles.(Muscles{iM}) = repmat(0,100,1)';
end
for iSubjects = 2:50
    %Load MVC Data
    cd(['J:\Piano_Fatigue\Data_Exported\EMG_Max_MVC'])
    load([MVCFilenames{1,iSubjects}]);
    
    % Load EMG data
    cd(['J:\Piano_Fatigue\Data_Exported\EMG_Ha'])
    load([HaFilenames{1,iSubjects}]);
    Data = EMG.data ;
    Freq = EMG.freq ;
        
    EMG = Data ;

        % Preprocessing
        % Freq : fréquence d'échantillonage de l'EMG
        [b,a] = butter(2,2*[10 400]/Freq) ; % Parametre du filtre BP 10-400 Hz
        EMGBP = filtfilt(b,a,EMG) ;
        EMGBS = EMGBP;
%         for iM = 1:47
%             plot(EMGBS(:,iM))
%             spec_fft(EMGBS(:,iM),Freq,1);
%             pause()
%         end
        
        % Filtre les données à 60Hz et 120Hz pour EMG_intra
        stop = [60 120 180];
        for iM = 1:length(Muscles)
            for iStop = 1:3
                [b,a] = butter(2,2*[stop(iStop)-1 stop(iStop)+1 ]/Freq,'stop') ;  % Parametre du filtre BS 55-65 Hz
                EMGBS(:,iM) = filtfilt(b,a,EMGBS(:,iM)) ;
            end
        end
        
        % Remove baseline
        EMGBL = EMGBS - repmat(mean(EMGBS),length(EMGBS),1) ; % Remove baseline (la fonction repmat donne une matrice de taille length(EMGBS) x 1 remplie de moyennes de EMGBS. On soustrait donc la moyenne de EMGBS à chacune de ses valeurs)
%         for iMuscles = 1:size(Muscles,1)
%             plot(EMGBL(:,iMuscles)) ; title(Muscles{iMuscles},'interpreter','none')
%             pause()
%         end
        
        %% Data analysis
        % EMG normalization
        Normalization = [];
        for iNormalization = 1:length(Muscles)
            Normalization (:,iNormalization) = EMGBL(:,iNormalization)./EMG_max(1,iNormalization)*100;
        end
%         for iMuscles = 1:size(Muscles,1)
%             subplot(2,1,1);plot(EMGBL(:,iMuscles)) ; title(Muscles{iMuscles},'interpreter','none')
%             subplot(2,1,2);plot(Normalization(:,iMuscles)) ; title(Muscles{iMuscles},'interpreter','none')
%             pause()
%         end

        %Reducing sample frequency to 1024 Hz
        Normalization2 = interp1(1:length(EMGBL),EMGBL,linspace(1,length(EMGBL),length(EMGBL)/2)) ;
%         subplot(2,1,1);plot(Normalization(:,1))
%         subplot(2,1,2);plot(Normalization2(:,1))

        % Compute EMG envelop
        [b,a] = butter(2,2*9/Freq) ; % Parametre du filtre Low pass 9 Hz
        EMG_envelop = filtfilt(b,a,abs(Normalization2)) ; % On ne garde que les valeurs absolues de EMGBL donc signal redressé au-dessus de 0
%         plot(EMG_envelop(:,1))
        
        %% EMG Median frequency
        %Libère mémoire
        clear Data EMG EMG_max EMGBL EMGBP EMGBS iHaFilenames iM iMuscles iMVCFilenames stop 
        clear functions
        
        cd(['J:\Piano_Fatigue\Data_Exported'])
%         load(['Cycle_Li_Beatrice.mat'])
%         Cy = cycles(iSubjects).seq*1024;
%         plot(Normalization2(:,1))
%         line([Cy(:,1) Cy(:,1)]', repmat(ylim,length(Cy(:,:)),1)','color','green')
%         line([Cy(:,2) Cy(:,2)]', repmat(ylim,length(Cy(:,:)),1)','color','red')

        load(['Cycle_Ha_Oceane.mat'])
        Cy = Cycle_Ha_Oceane(iSubjects).seq*1024;
%         plot(Normalization2(:,1))
%         line([Cy(:,1) Cy(:,1)]', repmat(ylim,length(Cy(:,:)),1)','color','green')
%         line([Cy(:,16) Cy(:,16)]', repmat(ylim,length(Cy(:,:)),1)','color','red')
        
        for iM = 1:47
            Normtemp = [];
            for iCy = 1:length(Cy)-2
                temp = EMG_envelop(round(Cy(iCy,1)):round(Cy(iCy,16)),iM);
                Normtemp(iCy,:) = interp1(1:length(temp),temp,linspace(1,length(temp),100));
            end
%             plot(Normtemp(iCy,:))
            MeanNormTemp = mean(Normtemp(:,:));
%             plot(MeanNormTemp)
            Act_Cycles.(Muscles{iM}) = [Act_Cycles.(Muscles{iM}); MeanNormTemp];
            I0 = mean(Act_Cycles.(Muscles{iM})(1:10));
            I1 = mean(Act_Cycles.(Muscles{iM})(11:20));
            I2 = mean(Act_Cycles.(Muscles{iM})(30:40));
            I3 = mean(Act_Cycles.(Muscles{iM})(50:60));
            I4 = mean(Act_Cycles.(Muscles{iM})(70:80));
            I5 = mean(Act_Cycles.(Muscles{iM})(90:100));
            
            Act_Cycles_Interval.(Muscles{iM}) = [Act_Cycles_Interval.(Muscles{iM}); iSubjects I0 I1 I2 I3 I4 I5];
        end
end
for iM = 1:length(Muscles)
    Act_Cycles.(Muscles{iM})(1,:) = [];
    Act_Cycles_Interval.(Muscles{iM})(1,:) = [];
end
cd(['J:\Piano_Fatigue\Data_Exported'])
% save('Act_Cycles_Interval_Li.mat','Act_Cycles_Interval')
% save('Act_Cycles_Ha.mat','Act_Cycles')

%% ANOVA repeated-measures
cd(['J:\Piano_Fatigue\Data_Exported'])
load('Act_Cycles_Interval_Li.mat')
for iM = 1:length(Muscles)
    for iInt = 2:7
        out = find(Act_Cycles_Interval.(Muscles{iM})(:,iInt) > median(Act_Cycles_Interval.(Muscles{iM})(:,iInt))+...
            1.5*iqr(Act_Cycles_Interval.(Muscles{iM})(:,iInt)));
        Act_Cycles_Interval.(Muscles{iM})(out,:) = [];
    end
end

RepeatedMeasuresResults = [];
GroupResults = [];
PostHocResults = [];
PostHocResults2 = [];
for iM = 1:length(Muscles)
    t = table(Act_Cycles_Interval.(Muscles{iM})(:,1),Act_Cycles_Interval.(Muscles{iM})(:,2),Act_Cycles_Interval.(Muscles{iM})(:,3),...
        Act_Cycles_Interval.(Muscles{iM})(:,4),Act_Cycles_Interval.(Muscles{iM})(:,5),Act_Cycles_Interval.(Muscles{iM})(:,6),...
        Act_Cycles_Interval.(Muscles{iM})(:,7),'VariableNames',{'Part','MFRef','MF10_20',...
        'MF30_40','MF50_60','MF70_80','MF90_100'});
    Meas = table([1 2 3 4 5 6]','VariableNames',{'Measurement'});
%     Meas2 = Meas;
%     Meas2.Measurement = categorical(Meas2.Measurement);
    rm = fitrm(t,'MFRef-MF90_100~Part','WithinDesign',Meas);
%     rm = fitrm(t,'MFRef-MF90_100~Group','WithinDesign',Meas);
    ranovatbl = ranova(rm);
    anovatbl = anova(rm);
    posthoc = multcompare(rm,'Measurement');
    RepeatedMeasuresResults.(Muscles{iM}) = ranovatbl;
    GroupResults.(Muscles{iM}) = anovatbl;
    PostHocResults.(Muscles{iM}) = posthoc([1 7 13 19 25],:);
    PostHocResults2.(Muscles{iM}) = posthoc([2 8 14 20],:);
end

Correction = [];
Correction2 = [];
Correction_interaction = [];
Correction_Group = [];
% Correction_DetailGroupResults = [];
for iM = 1:length(Muscles)
    Correction = [Correction; table2array(PostHocResults.(Muscles{iM})(:,5))];
    Correction2 = [Correction2; table2array(PostHocResults2.(Muscles{iM})(:,5))];
    Correction_interaction = [Correction_interaction; table2array(RepeatedMeasuresResults.(Muscles{iM})(1,5))]; %(1,5)_Time; (2,5)_Interaction
    Correction_Group = [Correction_Group; table2array(GroupResults.(Muscles{iM})(2,7))];
%     Correction_DetailGroupResults = [Correction_DetailGroupResults; DetailGroupResults.(Muscles{iM})];
end
[FDR,Qtime] = mafdr(Correction);
PCorrected = [];
seq = 1:5:length(Qtime);
for iM = 1:length(Muscles)
    PCorrected.(Muscles{iM}) = Qtime(seq(iM):seq(iM)+4);
end
[FDR,Qtime2] = mafdr(Correction2);
PCorrected2 = [];
seq = 1:4:length(Qtime2);
for iM = 1:length(Muscles)
    PCorrected2.(Muscles{iM}) = Qtime2(seq(iM):seq(iM)+3);
end
[Qinter] = mafdr(Correction_interaction,'BHFDR',true);
Qinter;
[Qgroup] = mafdr(Correction_Group,'BHFDR',true);
Qgroup;

    % Plot profil cycle
cd(['J:\Piano_Fatigue\Data_Exported'])
load('Act_Cycles_Li.mat')

for iM = 1:length(Muscles)
    for iInt = 1:100
        out = find(Act_Cycles.(Muscles{iM})(:,iInt) > median(Act_Cycles.(Muscles{iM})(:,iInt))+...
            1.5*iqr(Act_Cycles.(Muscles{iM})(:,iInt)));
        Act_Cycles.(Muscles{iM})(out,:) = [];
    end
end
    inc = 1;
    incM = 1;
    incm = 1;
    for iM = 42:-1:1
        Median_Act_Cycles.(Muscles{iM}) = median(Act_Cycles.(Muscles{iM}));
        IQR_Act_Cycles.(Muscles{iM}) = iqr(Act_Cycles.(Muscles{iM}));

        subplot(6,7,inc); plot(Median_Act_Cycles.(Muscles{iM}))
        axis([0 inf 0 450])
        hold on;
        line([0 100],[20 20],'color','green')
        if inc==36 || inc==37 || inc==38 || inc==39 || inc==40 || inc==41
            inc = incM + 1;
            incM = incM + 1;
        else inc = inc + 7;
        end
    end
    
    inc = 1;
    incM = 1;
    for iM = 1:42
        Median_Act_Cycles.(Muscles{iM}) = median(Act_Cycles.(Muscles{iM}));
        IQR_Act_Cycles.(Muscles{iM}) = iqr(Act_Cycles.(Muscles{iM}));

        subplot(6,7,inc); plot(Median_Act_Cycles.(Muscles{iM}))
        axis([0 inf 0 450])
        hold on;
        line([0 100],[20 20],'color','green')
        if inc==36 || inc==37 || inc==38 || inc==39 || inc==40 || inc==41
            inc = incM + 1;
            incM = incM + 1;
        else inc = inc + 7;
        end
    end
        
    % Plot map
    if graph == 1       
        %% Plot forearm map
        schem = [1:6;7:12;13:18;19:24;25:30;31:36;37:42]';
        
        Mean_Act_Int = [];
        for iInt = 1:6
            for iM = 1:length(Muscles)
                Mean_Act_Int(iM,iInt) = mean(Act_Cycles_Interval.(Muscles{iM})(:,iInt+1));
            end
        end
        iCol = 1;
        for iM = 1:6:42
            M1(:,iCol) = Mean_Act_Int(iM:iM+5,1);
            M2(:,iCol) = Mean_Act_Int(iM:iM+5,2);
            M3(:,iCol) = Mean_Act_Int(iM:iM+5,3);
            M4(:,iCol) = Mean_Act_Int(iM:iM+5,4);
            M5(:,iCol) = Mean_Act_Int(iM:iM+5,5);
            M6(:,iCol) = Mean_Act_Int(iM:iM+5,6);
            iCol = iCol +1;
        end
        
%         figure
        subplot(2,5,1);imagesc(M2-M1,[-120 250]) % [-15 15] Ha ; [-120 250] Li
        hold on
        for iM = 1:length(Muscles)-5
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(1) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end
        view([-180 90])

%         figure
        subplot(2,5,2);imagesc(M3-M1,[-120 250])
        hold on
        for iM = 1:length(Muscles)-5
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(2) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            elseif PCorrected2.(Muscles{iM})(1) < 0.05 & PCorrected.(Muscles{iM})(1) > 0.05
                plot(col-0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
                plot(col+0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end
        view([-180 90])

%         figure
        subplot(2,5,3);imagesc(M4-M1,[-120 250])
        hold on
        for iM = 1:length(Muscles)-5
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(3) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            elseif PCorrected2.(Muscles{iM})(2) < 0.05 & PCorrected.(Muscles{iM})(2) > 0.05
                plot(col-0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
                plot(col+0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end
        view([-180 90])

%         figure
        subplot(2,5,4);imagesc(M5-M1,[-120 250])
        hold on
        for iM = 1:length(Muscles)-5
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(4) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            elseif PCorrected2.(Muscles{iM})(3) < 0.05 & PCorrected.(Muscles{iM})(3) > 0.05
                plot(col-0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
                plot(col+0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end
        view([-180 90])

%         figure
        subplot(2,5,5);imagesc(M6-M1,[-120 250])
        hold on
        for iM = 1:length(Muscles)-5
            [row,col] = find(schem == iM);
            if PCorrected.(Muscles{iM})(5) < 0.05
                plot(col, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            elseif PCorrected2.(Muscles{iM})(4) < 0.05 & PCorrected.(Muscles{iM})(4) > 0.05
                plot(col-0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
                plot(col+0.05, row, 'black.', 'LineWidth', 2, 'MarkerSize', 10);
            end
        end
        view([-180 90])

        subplot(2,5,6);imagesc(M6-M1,[-120 250])
