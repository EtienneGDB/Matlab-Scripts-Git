clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Entropy Fouaz'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\boxcount'))

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
%         subplot(2,1,2);plot(Normalization2(:,:))
        
        %% Box-counting computation
        BCC = [];
        Seg = 1:1024*10:length(Normalization2);
        for iMuscles = 42:size(Muscles,1)
            inc = 1;
            for iSeg = 1:length(Seg)-1
        
            % Box-counting computation
            [n, r] = boxcount(Normalization2(Seg(iSeg):Seg(iSeg+1)-1,iMuscles),'slope');

            % log-log plot of the box counting
    %         loglog(r, n,'bo-', r, (r/r(end)).^(-2), 'r--')
    %         xlabel('r')
    %         ylabel('n(r)')
    %         legend('actual box-count','space-filling box-count');
            N = log(n)';
            R = log(r.^-1)';
%             scatter(R,N);
%             hold on;

            % Linear regression computation
            ini = 1; %initial box size
            fin = length(n)-1; %final box size
            Rr = R(ini:fin); %R limited to the selected range of boxes
            Nr = N(ini:fin); %R limited to the selected range of boxes
        
            x = [ones(length(Rr),1) Rr]; %Adds a column of ones to Rr
            b = x\Nr; %b(1) is the y-intercept and b(2) is the slope of the line
            y = x*b; %linear regression
            corr = 1-sum((Nr-y).^2)/sum((Nr-mean(Nr)).^2); %correlation

            BCC.(Muscles{iMuscles,1})(inc) = b(2);
            % Plots the regression line and the fractal dimension results
%             plot(Rr,y);
%             xlabel('log(1/r)');
%             ylabel('log(n)');
%             title(sprintf('Fractal Dimension: %f, Correlation: %f',b(2),corr));
            inc = inc+1;
            
            end
            plot(BCC.(Muscles{iMuscles,1}))
        end
        
        cd(['H:\Bureau\Etienne\Extracted data\Fatigue\FractalDimension'])
        save(['FractalDimension' (FileNames{2,iFiles}) '_' (Subjects{iSubjects}) '.mat'],'BCC')
        clear functions

    end
end