clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Entropy Fouaz'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\UnivarScatter-master'))

load('J:\Violon_SAB\results\Kin\allSubjects\Fatigue\Data')
Subjects = fatigueData.Subjects;
fatigueData = rmfield(fatigueData,'Subjects');
Variables = fieldnames(fatigueData);
Conditions = [1 2];

Interet = {'mean'};
VarInc = 1;
for iInteret = 1%:length(Interet)
    for iVariables = 1:length(Variables)
        if iVariables==1 | iVariables == 2 | iVariables == 3
            for iSubjects = 1:length(Subjects)
                for iConditions = 1:length(Conditions)
                    temp = fatigueData.(Variables{iVariables}).Subjects.All(~isnan(fatigueData.(Variables{iVariables}).Subjects.All(:,iSubjects,iConditions)),iSubjects,iConditions);
                    First(iSubjects,iConditions) = eval([Interet{iInteret} '(temp(1:5))']);
                    Last(iSubjects,iConditions) = eval([Interet{iInteret} '(temp(end-5:end-1))']);
                end
            end
            Y = [First(:,1); First(:,2); Last(:,1); Last(:,2)]; % Variable dépendante
            S = repmat([1:length(Subjects)],1,4)'; % Sujets
            F1 = [repmat(1,length(Subjects)*2,1); repmat(2,length(Subjects)*2,1)]; % Facteur temps
            F2 = [repmat(1,length(Subjects),1); repmat(2,length(Subjects),1);...
                repmat(1,length(Subjects),1); repmat(2,length(Subjects),1)]; % Facteur condition
            FactNames = {'Time', 'Condition'};
            
            test.(Interet{iInteret}).(Variables{iVariables}) = rm_anova2(Y,S,F1,F2,FactNames);
            
%             Colors = [0 0.4470 0.7410; 0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980];
%             UnivarScatter([First(:,1) First(:,2) Last(:,1) Last(:,2)], 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
%             xticklabels({'Without First' 'Without Last' 'With First' 'With Last'})
%             ylabel([Interet{iInteret} ' ' Variables{iVariables}])
%             if test.(Interet{iInteret}).(Variables{iVariables}){4,6} < 0.05
%                 hold on
%                 plot(2.5, max(Y), '*', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
%             end
%             if test.(Interet{iInteret}).(Variables{iVariables}){2,6} < 0.05
%                 hold on
%                 plot(1.5, max(Y)-std(Y), '+', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
%                 plot(3.5, max(Y)-std(Y), '+', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
%             end
%             if test.(Interet{iInteret}).(Variables{iVariables}){3,6} < 0.05
%                 hold on
%                 plot(2.5, max(Y)-std(Y), '-o', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
%             end
            
            Ycol = [First(:,1) First(:,2) Last(:,1) Last(:,2)]; % Variable dépendante
            TableResult(VarInc,:) = table(cellstr(Variables{iVariables}),cellstr(Variables{iVariables}),cellstr(Variables{iVariables}),cellstr(Variables{iVariables}),mean(Ycol),mean(Ycol(:,1:2),'all'),mean(Ycol(:,3:4),'all'),...
                mean(Ycol(:,1:3),'all'),mean(Ycol(:,2:4),'all'),test.(Interet{iInteret}).(Variables{iVariables}){2,6},...
                test.(Interet{iInteret}).(Variables{iVariables}){3,6},test.(Interet{iInteret}).(Variables{iVariables}){4,6},...
                'VariableNames',{'Variable','DoF','Position','RawCalculation','mean_F1_F2_L1_L2','mean_First','mean_Last','mean_Without',...
                'mean_With','p_Time','p_Condition','p_Interaction'});
            VarInc = VarInc+1;

%             pause()
%             close all
        else
            DOF = fieldnames(fatigueData.(Variables{iVariables}));
            for iDOF = 1:length(DOF)
                Position = fieldnames(fatigueData.(Variables{iVariables}).(DOF{iDOF}));
                for iPosition = 2%:length(Position)
                    RawCalculation = fieldnames(fatigueData.(Variables{iVariables}).(DOF{iDOF}).(Position{iPosition}));
                    for iRawCalc = 1%:length(RawCalculation)
                        for iSubjects = 1:length(Subjects)
                            for iConditions = 1:length(Conditions)
                                temp = fatigueData.(Variables{iVariables}).(DOF{iDOF}).(Position{iPosition}).(RawCalculation{iRawCalc}).Subjects.All(~isnan(fatigueData.(Variables{iVariables}).(DOF{iDOF}).(Position{iPosition}).(RawCalculation{iRawCalc}).Subjects.All(:,iSubjects,iConditions)),iSubjects,iConditions);
                                if ~isempty(temp)
                                    First(iSubjects,iConditions) = eval([Interet{iInteret} '(temp(1:5))']);
                                    Last(iSubjects,iConditions) = eval([Interet{iInteret} '(temp(end-5:end-1))']);
                                end
                            end
                        end
                    end
                    Y = [First(:,1); First(:,2); Last(:,1); Last(:,2)]; % Variable dépendante
                    S = repmat([1:length(Subjects)],1,4)'; % Sujets
                    F1 = [repmat(1,length(Subjects)*2,1); repmat(2,length(Subjects)*2,1)]; % Facteur temps
                    F2 = [repmat(1,length(Subjects),1); repmat(2,length(Subjects),1);...
                        repmat(1,length(Subjects),1); repmat(2,length(Subjects),1)]; % Facteur condition
                    FactNames = {'Time', 'Condition'};
                    
                    test.(Interet{iInteret}).(Variables{iVariables}).(DOF{iDOF}).(Position{iPosition}).(RawCalculation{iRawCalc}) = rm_anova2(Y,S,F1,F2,FactNames);
                    
%                     Colors = [0 0.4470 0.7410; 0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980];
%                     UnivarScatter([First(:,1) First(:,2) Last(:,1) Last(:,2)], 'LineWidth', 1.5, 'Whiskers', 'box', 'MarkerFaceColor', Colors, 'MeanColor', 'w', 'MarkerEdgeColor', 'none', 'SEMColor', 'black', 'StdColor', 'black');
%                     xticklabels({'Without First' 'Without Last' 'With First' 'With Last'})
%                     ylabel([Interet{iInteret} ' ' Variables{iVariables} ' ' DOF{iDOF} ' ' RawCalculation{iRawCalc}])
%                     if test.(Interet{iInteret}).(Variables{iVariables}).(DOF{iDOF}).(Position{iPosition}).(RawCalculation{iRawCalc}){4,6} < 0.05
%                         hold on
%                         plot(2.5, max(Y), '*', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
%                     end
%                     if test.(Interet{iInteret}).(Variables{iVariables}).(DOF{iDOF}).(Position{iPosition}).(RawCalculation{iRawCalc}){2,6} < 0.05
%                         hold on
%                         plot(1.5, max(Y)-std(Y), '+', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
%                         plot(3.5, max(Y)-std(Y), '+', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
%                     end
%                     if test.(Interet{iInteret}).(Variables{iVariables}).(DOF{iDOF}).(Position{iPosition}).(RawCalculation{iRawCalc}){3,6} < 0.05
%                         hold on
%                         plot(2.5, max(Y)-std(Y), '-o', 'color', 'black', 'LineWidth', 1, 'MarkerSize', 10);
%                     end
                    
                    Ycol = [First(:,1) First(:,2) Last(:,1) Last(:,2)]; % Variable dépendante
                    TableResult(VarInc,:) = table(cellstr(Variables{iVariables}),cellstr(DOF{iDOF}),...
                        cellstr(Position{iPosition}),cellstr(RawCalculation{iRawCalc}),mean(Ycol),...
                        mean(Ycol(:,1:2),'all'),mean(Ycol(:,3:4),'all'),mean(Ycol(:,1:3),'all'),mean(Ycol(:,2:4),'all'),...
                        test.(Interet{iInteret}).(Variables{iVariables}).(DOF{iDOF}).(Position{iPosition}).(RawCalculation{iRawCalc}){2,6},...
                        test.(Interet{iInteret}).(Variables{iVariables}).(DOF{iDOF}).(Position{iPosition}).(RawCalculation{iRawCalc}){3,6},...
                        test.(Interet{iInteret}).(Variables{iVariables}).(DOF{iDOF}).(Position{iPosition}).(RawCalculation{iRawCalc}){4,6},...
                        'VariableNames',{'Variable','DoF','Position','RawCalculation','mean_F1_F2_L1_L2','mean_First','mean_Last','mean_Without',...
                        'mean_With','p_Time','p_Condition','p_Interaction'});
                    VarInc = VarInc+1;

%                     pause()
%                     close all
                end
            end
        end
        
    end
end
