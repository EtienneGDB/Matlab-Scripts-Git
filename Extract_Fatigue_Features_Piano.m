clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Sample entropy Fouaz'))

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

%iSubjects=1;
for iSubjects = 45:50
    %Load MVC Data
    cd(['J:\Piano_Fatigue\Data_Exported\EMG_Max_MVC'])
    load([MVCFilenames{1,iSubjects}]);
    
    % Load EMG data
    cd(['J:\Piano_Fatigue\Data_Exported\EMG_Li'])
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
        Normalization2 = interp1(1:length(Normalization),Normalization,linspace(1,length(Normalization),length(Normalization)/2)) ;
%         subplot(2,1,1);plot(Normalization(:,1))
%         subplot(2,1,2);plot(Normalization2(:,1))
        
        % Dection of activity
        cd(['J:\Piano_Fatigue\Data_Exported'])
        load(['cycles_Li.mat'])
        participants(iSubjects).t_do = participants(iSubjects).t_do*(Freq/2);

%         plot(Normalization2(:,1))
%         line([participants(iSubjects).t_do,participants(iSubjects).t_do]', repmat(ylim,length(participants(iSubjects).t_do),1)','color','red')
        
        % Compute EMG envelop
        [b,a] = butter(2,2*9/Freq) ; % Parametre du filtre Low pass 9 Hz
        EMG_envelop = filtfilt(b,a,abs(Normalization2)) ; % On ne garde que les valeurs absolues de EMGBL donc signal redressé au-dessus de 0
%         plot(EMG_envelop(:,1))
%         line([participants(iSubjects).t_do,participants(iSubjects).t_do]', repmat(ylim,length(participants(iSubjects).t_do),1)','color','red')

        n=2000;
        Z = zeros(n/2);
        signal = [Z(:,1); EMG_envelop(:,1); Z(:,1)];
        Env = [];
        for i = 501:length(signal)-500
            Env(i,1) = mean(signal(i-499:i+499,1));
        end
%         plot(EMG_envelop(:,1))
%         hold on
%         plot(Env)
        
        % Compute Activity (feature)
        Activity = [];
        for iM = 1:length(Muscles)
            for iSeg = 1:(length(participants(iSubjects).t_do)-1)
                Activity.(Muscles{iM})(:,iSeg) = var(EMG_envelop(participants(iSubjects).t_do(iSeg):participants(iSubjects).t_do(iSeg+1),iM));
            end     
        end
        
        % Compute Mobility (feature)
        Mobility = [];
        for iM = 1:length(Muscles)
            % Derive Signal
            axT = [1:length(EMG_envelop(:,iM))]/Freq;
            Derivative = diff(EMG_envelop(:,iM))./diff(axT(1:2));
            
            for iSeg = 1:(length(participants(iSubjects).t_do)-1)
                Mobility.(Muscles{iM})(:,iSeg) = var(Derivative(participants(iSubjects).t_do(iSeg):participants(iSubjects).t_do(iSeg+1)))/var(EMG_envelop(participants(iSubjects).t_do(iSeg):participants(iSubjects).t_do(iSeg+1),iM));
            end     
        end
        
        % EMG activation level
        Amplitude = [];
        for iM = 1:length(Muscles)
            for iSeg = 1:(length(participants(iSubjects).t_do)-1)
                Amplitude.(Muscles{iM})(:,iSeg) = mean(EMG_envelop(participants(iSubjects).t_do(iSeg):participants(iSubjects).t_do(iSeg+1),iM));
            end
        end
        
        % Compute Sample Entropy (feature)
        SampleEntropy = [];
        for iM = 1:length(Muscles)
            for iSeg = 1:(length(participants(iSubjects).t_do)-1)
                SampleEntropy.(Muscles{iM})(:,iSeg) = sampenc(Normalization2(participants(iSubjects).t_do(iSeg):participants(iSubjects).t_do(iSeg+1),iM),2,0.2);
            end
        end

        cd(['J:\Piano_Fatigue\Data_Exported\Activity_Li'])
        save(['Activity' (HaFilenames{2,iSubjects}) '.mat'],'Activity')

        cd(['J:\Piano_Fatigue\Data_Exported\Mobility_Li'])
        save(['Mobility' (HaFilenames{2,iSubjects}) '.mat'],'Mobility')
        
        cd(['J:\Piano_Fatigue\Data_Exported\Amplitude_Li'])
        save(['Amplitude_' (HaFilenames{2,iSubjects}) '.mat'],'Amplitude')
        
        cd(['J:\Piano_Fatigue\Data_Exported\SampleEntropy_Li'])
        save(['SampleEntropy_' (HaFilenames{2,iSubjects}) '.mat'],'SampleEntropy')

end