clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))

% addpath H:\Bureau\Etienne\MATLAB\Functions
% addpath(genpath('H:\Bureau\Etienne\MATLAB\Functions\btk'))
% addpath(genpath('H:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))

FoldersNames = dir('F:\Data\IRSST\RAW');
Subjects = {...
    FoldersNames(4:36).name...
   } ;
Subjects(5) = [];
Subjects(5) = [];

cd(['H:\Projet_ExpertsNovices\excel'])
[num, txt, raw] = xlsread('participants.xlsx');
partic = participants ;

% Crée une 2eme ligne avec l'expertise du participant
Subjects{2,1}=[];
Subjects{3,1}=[];
for iSubjects = 1:length(Subjects)
    Subjects{2,iSubjects} = txt{find(contains(partic,Subjects(1,iSubjects)))+3,10};
    if Subjects{2,iSubjects} == 'Novice'
        Subjects{3,iSubjects} = 0;
    else
        Subjects{3,iSubjects} = 1;
    end
end

TimeData = [];
varInc = 1;
for iSubjects = 1:length(Subjects)
    % Plus de place sur serveur F donc les données ne sont pas au même
    % endroit pour tout le monde
    if ismember(iSubjects,[2:16 18:19 21:24 26:length(Subjects)])
        FolderContent = dir(['F:\Data\IRSST\RAW\' Subjects{1,iSubjects} '\fatigue']);
        FileNames = {...
        FolderContent(3:length(FolderContent)).name...
        } ;
    end
    if ismember(iSubjects,[1 17 19 22 25 28])
        if iSubjects == 1
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-09-29\aled\fatigue'];
        end
        if iSubjects == 17
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-11-27\jono\fatigue'];
        end
        if iSubjects == 19
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-11-22\marc\fatigue'];
        end
        if iSubjects == 22
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-11-10\nicl\fatigue'];
        end        
        if iSubjects == 25
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-12-12\pasd\fatigue'];
        end
        if iSubjects == 28
            pathData = ['H:\Projet_ExpertsNovices\data\raw\2017-09-28\samc\fatigue'];
        end
        FolderContent = dir(pathData);
        FileNames = {...
        FolderContent(3:length(FolderContent)).name...
        } ;
        for iFileNames = 1:length(FileNames)
            if iFileNames >= length(FileNames)
                while strcmp(FileNames{length(FileNames)}(length(FileNames{length(FileNames)})-2:length(FileNames{length(FileNames)})), 'c3d') == 0
                    FileNames(length(FileNames)) = [];
                end
                break
            end
            while strcmp(FileNames{iFileNames}(length(FileNames{iFileNames})-2:length(FileNames{iFileNames})), 'c3d') == 0
                FileNames(iFileNames) = [];
                if iFileNames >= length(FileNames)
                    break
                end
            end
        end
    end

    while (str2num(FileNames{1}(1)) == 1) | (str2num(FileNames{1}(1)) == 2)
        FileNames(1) = [];
    end
    
    %Crée une deuxième ligne avec le nom du fichier sans .c3d
    FileNames{2,1}=[];
    for iFileNames = 1:length(FileNames)
        FileNames{2,iFileNames} = erase(FileNames{1,iFileNames},'.c3d');
    end
    
    for iFiles = 1:length(FileNames)
        cd(['H:\Bureau\Etienne\Extracted data\Fatigue\Signal Segments'])
        load(['Seg_' (FileNames{2,iFiles}) '_' (Subjects{1,iSubjects}) '.mat'])

        SegNonAct = [];
        for iSeg = 1:length(Seg)-1
            SegNonAct(iSeg,1) = Seg(iSeg,2);
            SegNonAct(iSeg,2) = Seg(iSeg+1,1);
        end
        
        if FileNames{2,iFiles}(1) == 'l'
            C1 = 1;
        elseif FileNames{2,iFiles}(1) == 'm'
            C1 = 2;
        elseif FileNames{2,iFiles}(1) == 's'
            C1 = 3;
        end
        
        if FileNames{2,iFiles}(3) == 'd'
            C2 = 1;
        elseif FileNames{2,iFiles}(3) == 'u'
            C2 = 2;
        end
        
        for iSeg = 1:length(SegNonAct)
            TimeData(varInc,1) = iSubjects;
            TimeData(varInc,2) = Subjects{3,iSubjects};
            TimeData(varInc,3) = C1;
            TimeData(varInc,4) = C2;
            TimeData(varInc,5) = str2num(FileNames{2,iFiles}(5));
            TimeData(varInc,6) = iSeg;
            TimeData(varInc,7) = (SegNonAct(iSeg,2)-SegNonAct(iSeg,1))/2000;
            varInc = varInc + 1;
        end
    end
end

VarNames = {'Participant','Expertise','BoxSize','UpDown','Trial','Seg','Time'};
save('H:\Bureau\Etienne\Extracted data\Time_NonActivation.mat','TimeData')
save('H:\Bureau\Etienne\Extracted data\VarNames_LMM.mat','VarNames')


boxplot(TimeData(:,7),TimeData(:,6))





