clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))

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

% Define muscles used
% Column 1 = muscle name, colonne 2 = location in EMG matrix
%warning('vérifier que le numéro correspond bien pour le pec')
Muscles = {...
    'deltant',  'deltant_EMG1';...
    'deltmed',  'deltmed_EMG2';...
    'deltpost',    'deltpost_EMG3';...
    'biceps',    'biceps_EMG4';...
    'triceps',  'triceps_EMG5';...
    'uptrap',  'uptrap_EMG6';...
    'pect',  'pect_IM_EMG7';...
    'ssp',  'ssp_EMG8';... supra epineux
    'isp',  'isp_EMG9';... infra ep.
    'subs',       'subs_EMG10';... subscapulaire
    } ;

Trials = {...
    'l_d_';...
    'l_u_';...
    'm_d_';...
    'm_u_';...
    's_d_';...
    's_u_';...
    };

% Load TableData_LMM
load('H:\Bureau\Etienne\Extracted data\TableData_LMM.mat')

% Choix boîte
DAT = TableData;
for iM = 1:length(Muscles)
%   DAT.(Muscles{iM}) = DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,3)==1,:);...Grosses boîtes
%    DAT.(Muscles{iM}) = DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,3)==2,:);...Moyennes boîtes
    DAT.(Muscles{iM}) = DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,3)==3,:);...Petites boîtes
    
    DAT.(Muscles{iM}) = DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,5)<5,:);...4 Premiers essais
end

% Mean Activation level
MeanActLev = [];
for iM = 1:length(Muscles)
    MeanActLev(1,iM) = mean(DAT.(Muscles{iM})(:,7));
    MeanActLev(2,iM) = std(DAT.(Muscles{iM})(:,7));
end

% Order DAT per trial including Up/Down
DAT_Trial = DAT;
for iM = 1:length(Muscles)
    varInc = 1;
    DAT_Trial.(Muscles{iM}) = NaN;
    for iP = 1:31
        Vp = 1;
        for iT = 1:4
            for iUD = 2:-1:1
                V = DAT.(Muscles{iM})(DAT.(Muscles{iM})(:,1)==iP & DAT.(Muscles{iM})(:,5)==iT & DAT.(Muscles{iM})(:,4)==iUD,:);
                for iCol = 1:14
                    DAT_Trial.(Muscles{iM})(varInc,iCol) = mean(V(:,iCol));
                    DAT_Trial.(Muscles{iM})(varInc,6) = Vp;
                end
                Vp=Vp+1;
                varInc=varInc+1;
            end
        end
    end
end
for iM = 1:length(Muscles)
    DAT_Trial.(Muscles{iM}) = DAT_Trial.(Muscles{iM})(isfinite(DAT_Trial.(Muscles{iM})(:,1)),:);
end

%
ES = [];
for iM = 1:length(Muscles)
    for iSubjects = 1:31
        ES(:,) = 
        DAT_Trial.(Muscles{iM})(DAT_Trial.(Muscles{iM})(:,1)==iSubjects,:)
    end
end










