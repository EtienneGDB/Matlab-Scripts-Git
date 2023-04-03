clear all; 
close all; 
clc;

addpath(genpath('Z:\Bureau\Clara\UdeM\S2M\projetSAB\Analyse_Donnees\drawingFigures\'));
load('W:\Violon_SAB\results\Kin\allSubjects\Fatigue\Data.mat')


Analyses = {'highVSlowRatio', 'elbowROM_R', 'shiftSpeed', 'ROM'}; %{'highVSlowRatio', 'elbowROM_R', 'shiftSpeed', 'Time2Fatigue'};
yLabels = {'High:low position time ratio', 'Righ elbow range of motion (°)', 'Left elbow angular speeed (°/s)'}; %{'High:low position time ratio', 'Righ elbow range of motion (°)', 'Left elbow angular speeed (°/s)', 'Time (minutes)'};
Titles = {'High:Low position ratio with no DAS', 'Righ elbow range of motion with no DAS', 'Left elbow shift speed with no DAS'}; %{'High:Low position ratio with no DAS', 'Righ elbow range of motion with no DAS', 'Left elbow shift speed with no DAS', 'Time to fatigue'};
Blocks = {'First Block', 'Last Block'};
Conditions = {'No-DAS Condition', 'DAS Condition'};
Colors = [rgb('DodgerBlue'); rgb('Orange')];
legendLabels = {'No-DAS Condition', 'DAS Condition'};
% trunkDOF = {'Extension', 'lateralFlexion_R', 'Circumduction_L'};
% pelvisDOF = {'x', 'y', 'z'};
Stats = {'Median', 'Mean', 'Min', 'Max'}; 

for iAnalysis = 4:4+10
    if iAnalysis == 2 || iAnalysis == 3
        fatigueData.(Analyses{iAnalysis}).Subjects.firstBlock = fatigueData.(Analyses{iAnalysis}).Subjects.firstBlock * 180/pi;
        fatigueData.(Analyses{iAnalysis}).Subjects.lastBlock = fatigueData.(Analyses{iAnalysis}).Subjects.lastBlock * 180/pi;
    end
    
    figure; ax = gca; 
    for iSubject = 1:14
        
        for iTesting = 1:2
            i = 1;
            
            if iAnalysis == 1 || iAnalysis == 2 || iAnalysis == 3
                scatter(1, fatigueData.(Analyses{iAnalysis}).Subjects.firstBlock(iSubject,iTesting), 80, Colors(iTesting,:), 'filled'); title(Titles{iAnalysis}); hold on;
                scatter(2, fatigueData.(Analyses{iAnalysis}).Subjects.lastBlock(iSubject,iTesting),  80, Colors(iTesting,:), 'filled'); title(Titles{iAnalysis}); hold on;
                line([1 2], [fatigueData.(Analyses{iAnalysis}).Subjects.firstBlock(iSubject,iTesting) fatigueData.(Analyses{iAnalysis}).Subjects.lastBlock(iSubject,iTesting)], 'Color', Colors(iTesting,:))
                
            elseif iAnalysis >= 4
                nDOF = 0;
                k = strfind(alphaRef(1,:), alphaRef(1,i));
                for i=1:length(k)
                    if k{i} == 1
                        nDOF = nDOF + 1;
                    end                        
                end
                
                for iDOF = 1: nDOF
                    for iStat = 1:length(Stats)
                        data(:,:,iTesting) = horzcat(fatigueData.([alphaRef{1,iDOF} Analyses{iAnalysis}]).(trunkDOF{iDOF}).(Stats{iStat}).Subjects.firstBlock(:,iTesting) * 180/pi, fatigueData.(Analyses{iAnalysis}).(trunkDOF{iDOF}).(Stats{iStat}).Subjects.lastBlock(:,iTesting) * 180/pi);
                        subplot(3,4,i);
                        scatter([1 2], [data(iSubject,1,iTesting) data(iSubject,2,iTesting)], 80, Colors(iTesting,:), 'filled');  title(strcat({'Trunk '}, (trunkDOF{iDOF}), {' '} ,(Stats{iStat}))); hold on;
                        line([1 2], data(iSubject,:,iTesting) , 'Color', Colors(iTesting,:))
                        xlim([0 3])
                        xticks([1 2])
                        xticklabels(Blocks)
                        ylabel('Degrees (°)');
                        i = i+1;
                    end
                end
                
%             elseif iAnalysis == 5
%                 for iDOF = 1:length(pelvisDOF)
%                     for iStat = 1:length(Stats)
%                         dataPelvis(:,:,iTesting) = horzcat(fatigueData.(Analyses{iAnalysis}).(pelvisDOF{iDOF}).(Stats{iStat}).Subjects.firstBlock(:,iTesting) * 180/pi, fatigueData.(Analyses{iAnalysis}).(pelvisDOF{iDOF}).(Stats{iStat}).Subjects.lastBlock(:,iTesting) * 180/pi);
%                         subplot(3,4,i);
%                         scatter([1 2], [dataPelvis(iSubject,1,iTesting) dataPelvis(iSubject,2,iTesting)], 80, Colors(iTesting,:), 'filled');  title(strcat({'Pelvis '}, (pelvisDOF{iDOF}), {' '} ,(Stats{iStat}))); hold on;
%                         line([1 2], dataPelvis(iSubject,:,iTesting) , 'Color', Colors(iTesting,:))
%                         xlim([0 3])
%                         xticks([1 2])
%                         xticklabels(Blocks)
%                         ylabel('Degrees (°)');
%                         i = i+1;
%                     end
%                 end
            end
            
            %             xlim([0 3])
            %             xticks([1 2])
            %             set(ax, 'xticklabel', Blocks);
            if iAnalysis == 1 || iAnalysis == 2 || iAnalysis == 3
                ax.YLabel.String = yLabels{iAnalysis};
            end
            
        end
        
    end
    
end

% % Time to fatigue plot
% figure; ax = gca; 
% for iSubject = 1:15
%     for iCondition = 1:2
%         
%         scatter(iCondition, squeeze(fatigueData.Time2Fatigue.Subjects.All(1,iSubject,iCondition))', 80, Colors(iCondition,:), 'filled'); title(Titles{4}); hold on;
%         line([1 1.5], [squeeze(fatigueData.Time2Fatigue.Subjects.All(1,iSubject,1))' (squeeze(fatigueData.Time2Fatigue.Subjects.All(1,iSubject,1))' + ((squeeze(fatigueData.Time2Fatigue.Subjects.All(1,iSubject,2))' - squeeze(fatigueData.Time2Fatigue.Subjects.All(1,iSubject,1))') / 2))], 'Color', Colors(1,:))
%         line([1.5 2], [(squeeze(fatigueData.Time2Fatigue.Subjects.All(1,iSubject,1))' + ((squeeze(fatigueData.Time2Fatigue.Subjects.All(1,iSubject,2))' - squeeze(fatigueData.Time2Fatigue.Subjects.All(1,iSubject,1))') / 2)) squeeze(fatigueData.Time2Fatigue.Subjects.All(1,iSubject,2))'], 'Color', Colors(2,:))
% 
%         alpha(.5)
%         xlim([0 3])
%         ylim([0 31])
%         xticks([1 2])
%         set(ax, 'xticklabel', Conditions);
%         ax.YLabel.String = yLabels{iAnalysis};
%         
%     end
% end
% [p, h, stats] = signrank(squeeze(fatigueData.Time2Fatigue.Subjects.All(1,:,1))', squeeze(fatigueData.Time2Fatigue.Subjects.All(1,:,2))');