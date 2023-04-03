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

for iSubjects = 34:50
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
%         for iM = 1
%             plot(EMGBS(:,iM))
%             [FFT, f, Mean_Freq, Med_Freq] = spec_fft(EMGBS(:,iM),Freq,0);
%             plot(f,FFT)
%             line([Mean_Freq Mean_Freq],[0 1000],'color','red')
%             line([Med_Freq Med_Freq],[0 1000],'color','green')
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
%         for iM = 1
%             [FFT, f, Mean_Freq, Med_Freq] = spec_fft(EMGBS(:,iM),Freq,0);
%             figure;
%             plot(f,FFT)
%             line([Mean_Freq Mean_Freq],[0 1000],'color','red')
%             line([Med_Freq Med_Freq],[0 1000],'color','green')
%         end

        
        % Remove baseline
        EMGBL = EMGBS - repmat(mean(EMGBS),length(EMGBS),1) ; % Remove baseline (la fonction repmat donne une matrice de taille length(EMGBS) x 1 remplie de moyennes de EMGBS. On soustrait donc la moyenne de EMGBS à chacune de ses valeurs)
%         for iMuscles = 1:size(Muscles,1)
%             plot(EMGBL(:,iMuscles)) ; title(Muscles{iMuscles},'interpreter','none')
%             pause()
%         end
        
        %% Data analysis
        % EMG normalization
%         Normalization = [];
%         for iNormalization = 1:length(Muscles)
%             Normalization (:,iNormalization) = EMGBL(:,iNormalization)./EMG_max(1,iNormalization)*100;
%         end
%         for iMuscles = 1:size(Muscles,1)
%             subplot(2,1,1);plot(EMGBL(:,iMuscles)) ; title(Muscles{iMuscles},'interpreter','none')
%             subplot(2,1,2);plot(Normalization(:,iMuscles)) ; title(Muscles{iMuscles},'interpreter','none')
%             pause()
%         end

        %Reducing sample frequency to 1024 Hz
        Normalization2 = interp1(1:length(EMGBL),EMGBL,linspace(1,length(EMGBL),length(EMGBL)/2)) ;
%         subplot(2,1,1);plot(Normalization(:,1))
%         subplot(2,1,2);plot(Normalization2(:,:))
        
        %% EMG Median frequency
        %Libère mémoire
        clear Data EMG EMG_max EMGBL EMGBP EMGBS iHaFilenames iM iMuscles iMVCFilenames stop 
        clear functions
        
        % TFR
        FreqMin = 1 ;
        FreqMax = 400 ;
        Resolution = 1 ;... pas de 1Hz
        WaveNumber = 7 ;... nombre de pics dans l'ondelette de Morlet
        Args = WaveletParameters(FreqMin,FreqMax,Resolution,WaveNumber,Freq) ;
        FreqRange = length(FreqMin:Resolution:FreqMax) ;
        Nb_Interp_Pnts = 1 ;
                
        % Compute TFR
        % Cut movements
        TFR = [] ;
        for iM = 1:10%42%length(Muscles)
            disp('.')
            if sum(Normalization2(:,iM))~=0 & isnan(sum(Normalization2(:,iM)))==0
                tic
                [norm, Time, Wave_FreqS] = TimeFreqTransform(Normalization2(:,iM),Freq,Args,Nb_Interp_Pnts) ;
                toc
%                 TFR.TFR.(Muscles{iM}) = norm ;
                TFR.MedianFreq.(Muscles{iM})(:,:,:) = Compute_Median_Frequency(norm,Wave_FreqS) ;
                TFR.SpectralEntropy.(Muscles{iM})(:,:,:) = Compute_Spectral_Entropy(norm, Wave_FreqS) ;
            else
                TFR.TFR.(Muscles{iM})(:,:,:) = nan(400,500) ;
                TFR.MedianFreq.(Muscles{iM})(:,:,:) = nan(500,1)  ;
                TFR.SpectralEntropy.(Muscles{iM})(:,:,:) = nan(500,1)  ;
            end
        end
        
%         TFR.freq = Wave_FreqS ;
%         TFR.muscles = Muscles ;
%         MedianFreq.muscles = Muscles ;
        
        TFR_MedianFreq = {};
        TFR_MedianFreq.(['MedianFreq']) = TFR.MedianFreq ;
        
        TFR_SpectralEntropy = {};
        TFR_SpectralEntropy.(['SpectralEntropy']) = TFR.SpectralEntropy ;

        % Plot TFR/Median_Freq/FFT
%         for iM=1:length(Muscles)
%             figure(2);
%             subplot(2,1,1) ; imagesc(norm); title([HaFilenames{2,iSubjects} '__' Muscles{iM}])
% %             subplot(3,1,2) ; plot(TFR.MedianFreq.(Muscles{iM})(:,:),'.')
%             subplot(2,1,2) ; plot(Normalization(:,iM))
% %             pause
% %             close all
%         end

        % Data saving
%         cd(['J:\Piano_Fatigue\Data_Exported\TFR'])
%         save(['CleanData_TFR_' (HaFilenames{2,iSubjects}) '_.mat'],'TFR')

        cd(['J:\Piano_Fatigue\Data_Exported\TFR_MedianFreq'])
        save(['TFR_MedianFreq_' (HaFilenames{2,iSubjects}) '.mat'],'TFR_MedianFreq')

        cd(['J:\Piano_Fatigue\Data_Exported\TFR_SpectralEntropy'])
        save(['TFR_SpectralEntropy_' (HaFilenames{2,iSubjects}) '.mat'],'TFR_SpectralEntropy')

        TFR = [];
        TFR_MedianFreq = [];
        TFR_SpectralEntropy = [];
        EMG = [];
        EMGBL = [];
        EMGBP = [];
        EMGBS = [];
        Movement = [];
        Normalization = [];
        clear functions
end
    
% xxx = struct2cell(TFR_MedianFreq.MedianFreq);
% xxx2 = struct2cell(TFR.MedianFreq);
% yyy = {xxx{1:8},xxx2{1},xxx{9:end}};
% TFR_MedianFreq.MedianFreq = cell2struct(yyy',Muscles',1);
% 
% xxx = struct2cell(TFR_SpectralEntropy.SpectralEntropy);
% xxx2 = struct2cell(TFR.SpectralEntropy);
% yyy = {xxx{1:8},xxx2{1},xxx{9:end}};
% TFR_SpectralEntropy.SpectralEntropy = cell2struct(yyy',Muscles',1);

