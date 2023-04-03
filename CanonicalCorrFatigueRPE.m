clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))

%FoldersNames = dir('F:\Data\IRSST\RAW');
FoldersNames = dir('J:\IRSST_temp\RAW');
Subjects = {...
    FoldersNames(4:36).name...
   } ;
Subjects(5) = [];
Subjects(5) = [];

cd(['Z:\Projet_ExpertsNovices\excel'])
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
TableData = TableData_FracDim;
for iM = 1:length(Muscles)
    TableData.(Muscles{iM}) = TableData.(Muscles{iM})(:,1:14)
end

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

% Load RPE Data
cd(['Z:\Projet_ExpertsNovices\excel'])
[num, txt, raw] = xlsread('rpe.xlsx');
RPE = num(2:59,6:14);
Sbj = transpose(txt(3:2:59,2));
Sbj{2,1} = [];
Sbj1 = [];
for i = 1:29
    for j = 1:31
        if contains(Subjects{1,j},Sbj{1,i})
            Sbj1(i) = j;
        end
    end
end
Sbj1 = Sbj1';
Sbj1 = repmat(Sbj1,1,2)';
Sbj1 = Sbj1(:)';

% Table RPE
TableRPE = [transpose(Sbj1) RPE];
TableRPE(TableRPE(:,1)==0,:) = [];
TableRPE = sortrows(TableRPE,1);
TableRPE(11:12,:) = [];...2x sujet 6
TableRPE(27:28,:) = [];...2x sujet 16

RPE_Shoulder=[];
varInc=1;
for i = 2:2:52 %1:2:51 pour shoulder & 2:2:52 pour overall
    for j = 2:9
        RPE_Shoulder(varInc,1) = TableRPE(i,1);
        RPE_Shoulder(varInc,2) = TableRPE(i,j);
        varInc=varInc+1;
    end
end
% RPE_Shoulder(RPE_Shoulder(:,1)==10,:)=[];

% Spearman correlations
% iMus = 7;
% DAT_Shoulder = DAT_Trial;
% Vemg = [7:11 13];
% 
% for iMus = 1:length(Muscles)
%     RPE_SH = RPE_Shoulder;
%     Del = [5 11 14 18 25];...missing RPE data
%     DAT_Shoulder.(Muscles{iMus})(ismember(DAT_Shoulder.(Muscles{iMus})(:,1),Del),:)=[];
% 
%     X1 = unique(DAT_Shoulder.(Muscles{iMus})(:,1));
%     X2 = unique(RPE_SH(:,1));
%     N1 = numel(X1);
%     N2 = numel(X2);
%     count1 = zeros(N1,1);
%     count2 = zeros(N2,1);
%     for k = 1:N1
%         count1(k) = sum(DAT_Shoulder.(Muscles{iMus})(:,1)==X1(k));
%     end
%     CT1 = [ X1(:) count1 ];
%     for k = 1:N2
%         count2(k) = sum(RPE_SH(:,1)==X2(k));
%     end
%     CT2 = [ X2(:) count2 ];
%     Vsync = setdiff(CT2,CT1);
%     if length(Vsync)~=0
%         for iVsync = 1:length(Vsync)
%             DAT_Shoulder.(Muscles{iMus})(DAT_Shoulder.(Muscles{iMus})(:,1)==Vsync(iVsync),:)=[];
%             RPE_SH(RPE_SH(:,1)==Vsync(iVsync),:)=[];
%         end
%     end
%     X = unique(DAT_Shoulder.(Muscles{iMus})(:,1));
%     N = numel(X);
%     count = zeros(N,1);
%     for k = 1:N
%         count(k) = sum(DAT_Shoulder.(Muscles{iMus})(:,1)==X(k));
%     end
%     CT = [ X(:) count ];
%     
%     for iP = 1:length(CT)
%         if CT(iP,2)<8
%             DAT_Shoulder.(Muscles{iMus})(DAT_Shoulder.(Muscles{iMus})(:,1)==CT(iP,1),:)=[];
%             RPE_SH(RPE_SH(:,1)==CT(iP,1),:)=[];
%         end
%     end
%     % Spearman correlations
%     for iVemg = 1:6
%         [rho,pval] = corr(RPE_SH(:,2),DAT_Shoulder.(Muscles{iMus})(:,Vemg(iVemg)),'Type','Spearman');
%         RHO(iMus,iVemg) = rho;
%         PVAL(iMus,iVemg) = pval;
%         plot(RPE_SH(:,2),DAT_Shoulder.(Muscles{iMus})(:,Vemg(iVemg)),'.')
%         pause
%     end
%     % Regression
%     [b,bint,r,rint,stats] = regress(RPE_SH(:,2),DAT_Shoulder.(Muscles{iMus})(:,Vemg));
%     aaa = DAT_Shoulder.(Muscles{iMus})(:,[7:11 13]);
%     for i = 1:length(aaa)
%         RPE_Shoulder(i,3) = aaa(i,1)*b(1)+aaa(i,2)*b(2)+aaa(i,3)*b(3)+aaa(i,4)*b(4)+aaa(i,5)*b(5)+aaa(i,6)*b(6);
%     end
%     plot(RPE_Shoulder(:,3),RPE_Shoulder(:,2),'.')
%     pause
% end

% [A,B,r,U,V,stats] = canoncorr(DAT_Shoulder.(Muscles{iMus})(:,[7:11 13]),RPE_Shoulder(:,2));

% Regression
DAT_Shoulder = DAT_Trial;
DAT_SH = [];
Vemg = [7:11 13];

% Armonisation DAT_Shoulder
% Étape 1-3
MinSize = 1000;
for iMus = 1:length(Muscles)
    Sz = size(DAT_Shoulder.(Muscles{iMus}));
    if Sz(1) < MinSize
        MinSize = Sz(1);
        IdRef = iMus;
    end
end
XRef = unique(DAT_Shoulder.(Muscles{IdRef})(:,1));
NRef = numel(XRef);
countRef = zeros(NRef,1);
for k = 1:NRef
    countRef(k) = sum(DAT_Shoulder.(Muscles{IdRef})(:,1)==XRef(k));
end
CTRef = [ XRef(:) countRef ];
    
for iMus = 1:length(Muscles)
    X2 = unique(DAT_Shoulder.(Muscles{iMus})(:,1));
    N2 = numel(X2);
    count2 = zeros(N2,1);
    for k = 1:N2
        count2(k) = sum(DAT_Shoulder.(Muscles{iMus})(:,1)==X2(k));
    end
    CT2 = [ X2(:) count2 ];
    Vsync = setdiff(CT2,CTRef);
    if length(Vsync)~=0
        for iVsync = 1:length(Vsync)
            DAT_Shoulder.(Muscles{iMus})(DAT_Shoulder.(Muscles{iMus})(:,1)==Vsync(iVsync),:)=[];
        end
    end
end

% Étape 2-4
for iMus = 1:length(Muscles)
    RPE_SH = RPE_Shoulder;
    Del = [5 11 14 18 25];...missing RPE data
    DAT_Shoulder.(Muscles{iMus})(ismember(DAT_Shoulder.(Muscles{iMus})(:,1),Del),:)=[];

    X1 = unique(DAT_Shoulder.(Muscles{iMus})(:,1));
    X2 = unique(RPE_SH(:,1));
    N1 = numel(X1);
    N2 = numel(X2);
    count1 = zeros(N1,1);
    count2 = zeros(N2,1);
    for k = 1:N1
        count1(k) = sum(DAT_Shoulder.(Muscles{iMus})(:,1)==X1(k));
    end
    CT1 = [ X1(:) count1 ];
    for k = 1:N2
        count2(k) = sum(RPE_SH(:,1)==X2(k));
    end
    CT2 = [ X2(:) count2 ];
    Vsync = setdiff(CT2,CT1);
    if length(Vsync)~=0
        for iVsync = 1:length(Vsync)
            DAT_Shoulder.(Muscles{iMus})(DAT_Shoulder.(Muscles{iMus})(:,1)==Vsync(iVsync),:)=[];
            RPE_SH(RPE_SH(:,1)==Vsync(iVsync),:)=[];
        end
    end
    X = unique(DAT_Shoulder.(Muscles{iMus})(:,1));
    N = numel(X);
    count = zeros(N,1);
    for k = 1:N
        count(k) = sum(DAT_Shoulder.(Muscles{iMus})(:,1)==X(k));
    end
    CT = [ X(:) count ];
    
    for iP = 1:length(CT)
        if CT(iP,2)<8
            DAT_Shoulder.(Muscles{iMus})(DAT_Shoulder.(Muscles{iMus})(:,1)==CT(iP,1),:)=[];
            RPE_SH(RPE_SH(:,1)==CT(iP,1),:)=[];
        end
    end
end

% Suite
DAT_SH = [];
for iMus = 1:length(Muscles)
    DAT_SH = [DAT_SH DAT_Shoulder.(Muscles{iMus})(:,[7])];
end
% DAT_SH = [ones(152,1) DAT_SH];

[b,bint,r,rint,stats] = regress(RPE_SH(:,2),DAT_SH(:,:));
stats
[b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(DAT_SH(:,:),RPE_SH(:,2),'penter',0.05,'premove',0.1);

mdl = fitlm(DAT_SH(:,:),RPE_SH(:,2))
mdl = stepwiselm(DAT_SH(:,:),RPE_SH(:,2),'constant','Upper','linear','PEnter',0.05)

aaa = DAT_Shoulder.(Muscles{iMus})(:,[7:11 13]);
for i = 1:length(aaa)
    RPE_Shoulder(i,3) = aaa(i,1)*b(1)+aaa(i,2)*b(2)+aaa(i,3)*b(3)+aaa(i,4)*b(4)+aaa(i,5)*b(5)+aaa(i,6)*b(6);
end
plot(RPE_Shoulder(:,3),RPE_Shoulder(:,2),'.')



