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
for iSubjects = 1:50
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

        % EMG normalization
        Normalization = [];
        for iNormalization = 1:length(Muscles)
            Normalization (:,iNormalization) = EMGBL(:,iNormalization)./EMG_max(1,iNormalization)*100;
        end
        
        % Reducing sample frequency to 1024 Hz
        Normalization2 = interp1(1:length(Normalization),Normalization,linspace(1,length(Normalization),length(Normalization)/2)) ;

        % Dection of activity
        cd(['J:\Piano_Fatigue\Data_Exported'])
%         load(['cycles_Li.mat'])
%         load(['Cycle_Li_Beatrice.mat'])
%         cycles(12).seq = cycles(12).seq;
        load(['Cycle_Ha_Oceane.mat'])
        Cycle_Ha_Oceane(iSubjects).seq = Cycle_Ha_Oceane(iSubjects).seq*60;
%         participants(iSubjects).t_do = participants(iSubjects).t_do*(Freq/2);

        for i = 1%:7:length(Muscles)
            t = 0:1/2048:(length(EMGBL)-1)/2048;
            plot(t,EMGBL(:,i)) ; title(['Subject',num2str(iSubjects),'_',Muscles{i}],'interpreter','none')
%             line([participants(iSubjects).t_do,participants(iSubjects).t_do]', repmat(ylim,length(participants(iSubjects).t_do),1)','color','green')
            % Li_Beatrice
%             line([cycles(12).seq(:,1),cycles(12).seq(:,1)]', repmat(ylim,length(cycles(12).seq(:,1)),1)','color','red')
%             line([cycles(12).seq(:,2),cycles(12).seq(:,2)]', repmat(ylim,length(cycles(12).seq(:,2)),1)','color','red')
            % Ha_Oceane
            line([Cycle_Ha_Oceane(21).seq(:,1),Cycle_Ha_Oceane(21).seq(:,1)]', repmat(ylim,length(Cycle_Ha_Oceane(21).seq(:,1)),1)','color','red')
            line([Cycle_Ha_Oceane(21).seq(:,16),Cycle_Ha_Oceane(21).seq(:,16)]', repmat(ylim,length(Cycle_Ha_Oceane(21).seq(:,16)),1)','color','red')
            pause
        end
        
    end
    close all
    