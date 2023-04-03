clear ;
close all ;
clc ;

addpath E:\Bureau\Etienne\MATLAB\Functions
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\btk'))
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))

addpath H:\Bureau\Etienne\MATLAB\Functions
addpath(genpath('H:\Bureau\Etienne\MATLAB\Functions\btk'))
addpath(genpath('H:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))

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
    
% Define muscles used
% Column 1 = muscle name, colonne 2 = location in EMG matrix
%warning('vérifier que le numéro correspond bien pour le pec')
Muscles = {...
    'm1_2';'m2_3';'m3_4';'m4_5';'m5_6';'m6_7';'m8_9';'m9_10';'m10_11';'m11_12';'m12_13';'m13_14';'m15_16';'m16_17';'m17_18';...
    'm18_19';'m19_20';'m20_21';'m22_23';'m23_24';'m24_25';'m25_26';'m26_27';'m27_28';'m29_30';'m30_31';'m31_32';'m32_33';...
    'm33_34';'m34_35';'m36_37';'m37_38';'m38_39';'m39_40';'m40_41';'m41_42';'m43_44';'m44_45';'m45_46';'m46_47';'m47_48';...
    'm48_49';'Biceps';'Triceps';'DeltAnt';'DeltMed';'SupTrap'};

%iSubjects=1;
figure
figure
for iSubjects = 7:50
    disp(iSubjects)
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
        inc = 1;
        for iM = 1%:7:47
%             subplot(2,3,iM) ; plot([0:length(EMGBS(:,iM))-1]/Freq,EMGBS(:,iM))
            [FFT, f] = spec_fft(EMGBS(:,iM),Freq,1);title(Muscles{iM},'interpreter','none')
%             subplot(2,3,inc) ; 
%             plot(f,FFT)
%             inc = inc + 1;
%             pause()
%             close all
        end
        pause
end
iM=26
            [FFT, f] = spec_fft(EMGBS(:,iM),Freq,0);title(Muscles{iM},'interpreter','none')
            subplot(2,3,4) ; plot(f,FFT)
iM=41
            [FFT, f] = spec_fft(EMGBS(:,iM),Freq,0);title(Muscles{iM},'interpreter','none')
            subplot(2,3,6) ; plot(f,FFT)

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
%             plot([0:length(EMGBS(:,iM))-1]/Freq,EMGBL(:,iMuscles)) ; title(Muscles{iMuscles},'interpreter','none')
%             pause()
%         end
        
        % Spectre de fréquence raw
%         for iM = 1:47
%             figure(1) ; plot([0:length(EMGBL(:,iM))-1]/Freq,EMGBL(:,iM)) ; title(['Subject',num2str(iSubjects),'_',Muscles{iM}],'interpreter','none')
%             [FFT,PSD,f] = spec_fft(EMGBL(:,iM),Freq,0);
%             figure(2) ; plot(f,FFT) ; title(['FFT']) ; xlabel('Frequency (Hz)') ; ylabel('FFT')
%             pause()
%         end
inc=1;
        for iM = 1:7:47
%             subplot(2,3,iM) ; plot([0:length(EMGBS(:,iM))-1]/Freq,EMGBS(:,iM))
            [FFT, f] = spec_fft(EMGBL(:,iM),Freq,0);title(Muscles{iM},'interpreter','none')
            subplot(2,3,inc) ; plot(f,FFT)
            inc = inc + 1;
%             pause()
%             close all
        end
        iM=26
            [FFT, f] = spec_fft(EMGBL(:,iM),Freq,0);title(Muscles{iM},'interpreter','none')
            subplot(2,3,4) ; plot(f,FFT)
iM=41
            [FFT, f] = spec_fft(EMGBL(:,iM),Freq,0);title(Muscles{iM},'interpreter','none')
            subplot(2,3,6) ; plot(f,FFT)

        pause
    end
    close all