close all ;
clear all ;
clc ;

addpath E:\Bureau\Etienne\MATLAB\Functions
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\btk'))
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))

addpath H:\Bureau\Etienne\MATLAB\Functions
addpath(genpath('H:\Bureau\Etienne\MATLAB\Functions\btk'))
addpath(genpath('H:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))

FoldersNames = dir('J:\Piano_Fatigue\Data_Exported\EMG_MVC');
MVCFilenames = {...
    FoldersNames(3:52).name...
   } ;

Muscles = {...
    '1-2';'2-3';'3-4';'4-5';'5-6';'6-7';'8-9';'9-10';'10-11';'11-12';'12-13';'13-14';'15-16';'16-17';'17-18';...
    '18-19';'19-20';'20-21';'22-23';'23-24';'24-25';'25-26';'26-27';'27-28';'29-30';'30-31';'31-32';'32-33';...
    '33-34';'34-35';'36-37';'37-38';'38-39';'39-40';'40-41';'41-42';'43-44';'44-45';'45-46';'46-47';'47-48';...
    '48-49';'Biceps';'Triceps';'DeltAnt';'DeltMed';'SupTrap'};

for iSubjects = 1:50
    cd(['J:\Piano_Fatigue\Data_Exported\EMG_MVC'])
    load([MVCFilenames{iSubjects}]);
    Freq = MVC(1).freq;

    EMG1 = transpose(MVC(1).data);
    EMG1(EMG1(:,1)==0,:) = [];
    EMG2 = transpose(MVC(2).data);
    EMG2(EMG2(:,1)==0,:) = [];
%     plot(EMG1(:,4))

% Preprocessing
    [b,a] = butter(2,2*[10 400]/Freq) ; % Parametre du filtre BP (Band Pass) 10-400 Hz
    EMGBP1 = filtfilt(b,a,EMG1) ;
    EMGBP2 = filtfilt(b,a,EMG2) ;

    [b,a] = butter(2,2*[55 65]/Freq,'stop') ;  % Parametre du filtre BS (Band Stop) 55-65 Hz (prises électriques)
    EMGBS1 = filtfilt(b,a,EMGBP1) ;
    EMGBS2 = filtfilt(b,a,EMGBP2) ;
    EMGBL1 = EMGBS1 - repmat(mean(EMGBS1),length(EMGBS1),1) ; % Remove baseline (la fonction repmat donne une matrice de taille length(EMGBS) x 1 remplie de moyennes de EMGBS. On soustrait donc la moyenne de EMGBS à chacune de ses valeurs)
    EMGBL2 = EMGBS2 - repmat(mean(EMGBS2),length(EMGBS2),1) ; % Remove baseline (la fonction repmat donne une matrice de taille length(EMGBS) x 1 remplie de moyennes de EMGBS. On soustrait donc la moyenne de EMGBS à chacune de ses valeurs)

%     for iM = 1:47
%         plot(EMGBL1(:,iM)) ; title (Muscles{iM}, 'interpreter', 'none') ;
%         hold on
%         plot(EMGBL2(:,iM)) ; title (Muscles{iM}, 'interpreter', 'none') ;
%         pause()
%     end
    
    %Concatenate signals
    EMG_allMVC = [EMGBL1 ; EMGBL2];
%     for iM = 1:47
%         plot(EMG_allMVC(:,iM)) ; title (Muscles{iM}, 'interpreter', 'none') ;
%         pause()
%     end
    
    % Data analysis
    % EMG activation level
    [b,a] = butter(2,2*9/Freq) ; % Parametre du filtre Low pass 9 Hz
    EMG_MVCenvelop = filtfilt(b,a,abs(EMG_allMVC)) ; %abs : on ne garde que les valeurs absolues de EMGBL donc signal redressé au-dessus de 0
    
%     for iM = 1:47
%         plot(EMG_MVCenvelop(:,iM)) ; title (Muscles{iM}, 'interpreter', 'none') ;
%         pause()
%     end
    
    %Sort EMG values from highest to lowest (descending order)
    EMG_sort = sort(EMG_MVCenvelop,'descend') ;
%     plot(EMG_sort(:,1))
%     title('Enveloppe MVC ordonnée'); xlabel('frame'); ylabel('EMG (mV)')
    
    %Mean EMG over 2000 Hz for each muscle
    EMG_sort = EMG_sort(1:Freq,:) ;
    EMG_max = mean(EMG_sort,1);
    
    % Data saving
    cd(['J:\Piano_Fatigue\Data_Exported\EMG_Max_MVC'])
    save(['EMG_Max_' MVCFilenames{iSubjects}],'EMG_max')
    clear functions
end
