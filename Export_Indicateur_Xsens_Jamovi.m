% plot the indicateur of the parameter wanted of all subject regarding the
% task

clear all;close all;clc;

Subjects = {'P01' 'P02' 'P03' 'P04' 'P05' 'P06' 'P07' 'P08' 'P09' 'P10' 'P11' 'P12' 'P13' 'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20' 'P21'...
    'P22' 'P23' 'P24' 'P25' 'P26' 'P27' 'P28' 'P29' 'P30' 'P31'};

cd(['\\10.89.24.15\q\IRSST_DavidsTea'])
load(['SegmentName_XSENS.mat'])

parametre  = {'Acceleration','Angular_Velocity','Jerk','Module_Acceleration',...
    'Module_Angular_Velocity','Module_Jerk'};

indicateur = {'MedianFreq','SpectralEntropy','PowerLF','PowerHF','PowerTot', ...
    'PeakPower','PeakPower_Freq'};

indicateur_time = {'Peak', 'Mean', 'STD', 'Median', 'IQR', 'Percentile10',...
    'Percentile25', 'Percentile75', 'Percentile90'};

Task =  'SP'; % RPT / Work / SP
showplot = false;

ipar = 1;
iind = 3;
coordonnees= 'x'; % x, y ou z
if coordonnees == 'x'
    start = 1;
elseif coordonnees == 'y'
    start = 32;
else
    start = 63;
end


TFR_Features = true;
Time_Features = true;


%% TFR_Features
if TFR_Features
    MatJamovi = [];
    Pvalues = [];
    Axes = {'x', 'y', 'z'};
    for ipar = 1:length(parametre)
        for iind = 1:length(indicateur)
            for iaxe = 1:3
                coordonnees = Axes{iaxe};
                if coordonnees == 'x'
                    start = 1;
                elseif coordonnees == 'y'
                    start = 32;
                else
                    start = 63;
                end
                
                cd(['\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted\TFR_XSENS\' parametre{ipar}])
                load(['TFR_' Task '_' parametre{ipar} '_' indicateur{iind} '.mat' ])
                Mat_Load = eval(['Mat_' parametre{ipar} '_' indicateur{iind}]);
                
                Norm = Mat_Load;
                if  ipar > 3
                    interval = 1:31;
                else
                    interval = 1:93;
                end
                
                % Normalization of the indicator chosen between trial by subtracting the
                % the trial by the first one
                for iSuj = interval
                    for iS = 1:length(SegmentName)
                        FirstTrial = 2;
                        
                        
                        for iTrial = 2:11
                            if (isnan(Norm.(SegmentName{iS,1})(iSuj,iTrial))) & (iTrial == 2)
                                inc = iTrial + 1;
                                while isnan(Norm.(SegmentName{iS,1})(iSuj,inc))
                                    if inc == 11
                                        break
                                    end
                                    inc = inc + 1;
                                end
                                FirstTrial = inc;
                            end
                            
                            NormPlot.(SegmentName{iS,1})(iSuj,iTrial) = Mat_Load.(SegmentName{iS,1})(iSuj,FirstTrial) - Norm.(SegmentName{iS,1})(iSuj,iTrial);
                            Norm.(SegmentName{iS,1})(iSuj,iTrial) = Norm.(SegmentName{iS,1})(iSuj,iTrial);
                            
                        end
                    end
                end
                
                
                for iS = 1:length(SegmentName)
                    vect = [1:2  4:6 8:13 15:26 28:31];
                    if ipar == 4 | ipar == 5 | ipar == 6
                        DataMat.(parametre{ipar}).(indicateur{iind}).(SegmentName{iS,1}) = Norm.(SegmentName{iS,1});
                        
                                        % box plot of all the subject by Segment
%                                         for iS = 1:length(SegmentName)
                                            % 2 figures
                                            if showplot == true
                                            if iS <10
                                                figure(1)
                                                sgtitle([parametre{ipar} ' ' indicateur{iind} ' ' Task])
                                                subplot(3,3,iS);
                                                plotver2 = NormPlot.(SegmentName{iS,1});
                                                plotver2 = plotver2(:,2:end);
                                                boxplot(plotver2,'DataLim',[min(prctile(plotver2,10)) max(prctile(plotver2,90))])
                                                line([1 10],[0 0],'color','black')
                                                title(SegmentName{iS},'interpreter','none')
                                            else
                                                figure(2)
                                                sgtitle([parametre{ipar} ' ' indicateur{iind} ' ' Task])
                                                subplot(3,3,iS-9);
                                                plotver2 = NormPlot.(SegmentName{iS,1});
                                                plotver2 = plotver2(:,2:end);
                                                boxplot(plotver2,'DataLim',[min(prctile(plotver2,10)) max(prctile(plotver2,90))])
                                                line([1 10],[0 0],'color','black')
                                                title(SegmentName{iS},'interpreter','none')
                                            end
                                            end
%                                         end
                        
                    else
                        DataMat.(parametre{ipar}).(indicateur{iind}).(SegmentName{iS,1}).(Axes{iaxe}) = Norm.(SegmentName{iS,1})(start:start+30,:);
                        %                 DataMat = DataMat(start:start+30,:);
                        
                                        % box plot of all the subject by Segment
%                                         for iS = 1:length(SegmentName)
                                            % 2 figures
                                            if showplot == true
                                            if iS <10
                                                figure(1)
                                                sgtitle([parametre{ipar} ' ' coordonnees ' ' indicateur{iind} ' ' Task])
                                                subplot(3,3,iS);
                                                plotver2 = NormPlot.(SegmentName{iS,1});
                                                plotver2 = plotver2(start:start+30,2:end);
                                                boxplot(plotver2,'DataLim',[min(prctile(plotver2,10)) max(prctile(plotver2,90))])
                                                line([1 10],[0 0],'color','black')
                                                title(SegmentName{iS},'interpreter','none')
                                            else
                                                figure(2)
                                                sgtitle([parametre{ipar} ' ' coordonnees ' ' indicateur{iind} ' ' Task])
                                                subplot(3,3,iS-9);
                                                plotver2 = NormPlot.(SegmentName{iS,1});
                                                plotver2 = plotver2(start:start+30,2:end);
                                                boxplot(plotver2,'DataLim',[min(prctile(plotver2,10)) max(prctile(plotver2,90))])
                                                line([1 10],[0 0],'color','black')
                                                title(SegmentName{iS},'interpreter','none')
                                            end
                                            end
%                                         end
                        
                    end
                end
                % pause()
                % clf(figure(1))
                % clf(figure(2))
                
                %---Export Jamovi---
                for iSegment = 1:length(SegmentName)
                    inc = 1;
                    for iSubjects = 1:31
                        for iTrial = 1:10
                            
                            if ipar == 4 | ipar == 5 | ipar == 6
                                MatJamovi.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1})(inc,1) = iSubjects;
                                MatJamovi.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1})(inc,2) = iTrial;
                                MatJamovi.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1})(inc,3) = ...
                                    DataMat.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1})(iSubjects, iTrial+1);
                            else
                                MatJamovi.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1}).(Axes{iaxe})(inc,1) = iSubjects;
                                MatJamovi.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1}).(Axes{iaxe})(inc,2) = iTrial;
                                MatJamovi.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1}).(Axes{iaxe})(inc,3) = ...
                                    DataMat.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1}).(Axes{iaxe})(iSubjects, iTrial+1);
                            end
                            inc = inc +1;
                        end
                    end
                end
                
                %---Linear Mixed Model---
                for iSegment = 3:8%length(SegmentName)
                    if ipar == 4 | ipar == 5 | ipar == 6
                        tbl = table(MatJamovi.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1})(:,1),...
                            categorical(MatJamovi.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1})(:,2)),...
                            MatJamovi.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1})(:,3),...
                            'VariableNames',{'Subject','Trial','Measure'});
                    else
                        tbl = table(MatJamovi.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1}).(Axes{iaxe})(:,1),...
                            categorical(MatJamovi.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1}).(Axes{iaxe})(:,2)),...
                            MatJamovi.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1}).(Axes{iaxe})(:,3),...
                            'VariableNames',{'Subject','Trial','Measure'});
                    end
                    lme = fitlme(tbl,'Measure~Trial+(1|Subject)');
                    for i = 2:10
                        if double(lme.Coefficients(i,6)) < 0.05
                            Pvalues.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1}).(Axes{iaxe}) = double(lme.Coefficients(i,6));
                            break
                        else
                            Pvalues.(parametre{ipar}).(indicateur{iind}).(SegmentName{iSegment,1}).(Axes{iaxe}) = 1;
                        end
                    end
                end
                
                
            end
        end
    end
    
    cd(['\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted'])
    save(['MatJamovi_', Task, '.mat'],'MatJamovi')
    save(['Pvalues_LM_', Task, '.mat'],'Pvalues')
end


%% Time_Features
if Time_Features
    MatJamovi = [];
    Pvalues = [];
    Axes = {'x', 'y', 'z'};
    cd(['\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted\Time_Domain_XSENS'])
    load(['TimeXSENS_' Task '.mat'])

    for ipar = 1:length(parametre)
        for iind_time = 1:length(indicateur_time)
            for iaxe = 1:3
                coordonnees = Axes{iaxe};
                if coordonnees == 'x'
                    start = 1;
                elseif coordonnees == 'y'
                    start = 32;
                else
                    start = 63;
                end
                
                % cd(['\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted\TFR_XSENS\' parametre{ipar}])
                % load(['TFR_' Task '_' parametre{ipar} '_' indicateur{iind} '.mat' ])
                % Mat_Load = eval(['Mat_' parametre{ipar} '_' indicateur{iind}]);
                
                Mat_Load = Struc_Indicateur_Time.(parametre{ipar}).(indicateur_time{iind_time});

                Norm = Mat_Load;
                if  ipar > 3
                    interval = 1:31;
                else
                    interval = 1:93;
                end
                
                % Normalization of the indicator chosen between trial by subtracting the
                % the trial by the first one
                for iSuj = interval
                    for iS = 1:length(SegmentName)
                        FirstTrial = 2;
                        
                        
                        for iTrial = 2:11
                            if (isnan(Norm.(SegmentName{iS,1})(iSuj,iTrial))) & (iTrial == 2)
                                inc = iTrial + 1;
                                while isnan(Norm.(SegmentName{iS,1})(iSuj,inc))
                                    if inc == 11
                                        break
                                    end
                                    inc = inc + 1;
                                end
                                FirstTrial = inc;
                            end
                            
                            NormPlot.(SegmentName{iS,1})(iSuj,iTrial) = Mat_Load.(SegmentName{iS,1})(iSuj,FirstTrial) - Norm.(SegmentName{iS,1})(iSuj,iTrial);
                            Norm.(SegmentName{iS,1})(iSuj,iTrial) = Norm.(SegmentName{iS,1})(iSuj,iTrial);
                            
                        end
                    end
                end
                
                
                for iS = 1:length(SegmentName)
                    vect = [1:2  4:6 8:13 15:26 28:31];
                    if ipar == 4 | ipar == 5 | ipar == 6
                        DataMat.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iS,1}) = Norm.(SegmentName{iS,1});
                        
                                        % box plot of all the subject by Segment
                                        % for iS = 1:length(SegmentName)
                                            % 2 figures
                                            if showplot == true
                                            if iS <10
                                                figure(1)
                                                sgtitle([parametre{ipar} ' ' indicateur_time{iind_time} ' ' Task])
                                                subplot(3,3,iS);
                                                plotver2 = NormPlot.(SegmentName{iS,1});
                                                plotver2 = plotver2(:,2:end);
                                                boxplot(plotver2,'DataLim',[min(prctile(plotver2,10)) max(prctile(plotver2,90))])
                                                line([1 10],[0 0],'color','black')
                                                title(SegmentName{iS},'interpreter','none')
                                            else
                                                figure(2)
                                                sgtitle([parametre{ipar} ' ' indicateur_time{iind_time} ' ' Task])
                                                subplot(3,3,iS-9);
                                                plotver2 = NormPlot.(SegmentName{iS,1});
                                                plotver2 = plotver2(:,2:end);
                                                boxplot(plotver2,'DataLim',[min(prctile(plotver2,10)) max(prctile(plotver2,90))])
                                                line([1 10],[0 0],'color','black')
                                                title(SegmentName{iS},'interpreter','none')
                                            end
                                            end
                                        % end
                        
                    else
                        DataMat.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iS,1}).(Axes{iaxe}) = Norm.(SegmentName{iS,1})(start:start+30,:);
                        %                 DataMat = DataMat(start:start+30,:);
                        
                                        % box plot of all the subject by Segment
                                        % for iS = 1:length(SegmentName)
                                            % 2 figures
                                            if showplot == true
                                            if iS <10
                                                figure(1)
                                                sgtitle([parametre{ipar} ' ' coordonnees ' ' indicateur_time{iind_time} ' ' Task])
                                                subplot(3,3,iS);
                                                plotver2 = NormPlot.(SegmentName{iS,1});
                                                plotver2 = plotver2(start:start+30,2:end);
                                                boxplot(plotver2,'DataLim',[min(prctile(plotver2,10)) max(prctile(plotver2,90))])
                                                line([1 10],[0 0],'color','black')
                                                title(SegmentName{iS},'interpreter','none')
                                            else
                                                figure(2)
                                                sgtitle([parametre{ipar} ' ' coordonnees ' ' indicateur_time{iind_time} ' ' Task])
                                                subplot(3,3,iS-9);
                                                plotver2 = NormPlot.(SegmentName{iS,1});
                                                plotver2 = plotver2(start:start+30,2:end);
                                                boxplot(plotver2,'DataLim',[min(prctile(plotver2,10)) max(prctile(plotver2,90))])
                                                line([1 10],[0 0],'color','black')
                                                title(SegmentName{iS},'interpreter','none')
                                            end
                                            end
                                        % end
                        
                    end
                end
                % pause()
                % clf(figure(1))
                % clf(figure(2))

                %---Export Jamovi---
                for iSegment = 1:length(SegmentName)
                    inc = 1;
                    for iSubjects = 1:31
                        for iTrial = 1:10
                            
                            if ipar == 4 | ipar == 5 | ipar == 6
                                MatJamovi.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1})(inc,1) = iSubjects;
                                MatJamovi.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1})(inc,2) = iTrial;
                                MatJamovi.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1})(inc,3) = ...
                                    DataMat.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1})(iSubjects, iTrial+1);
                            else
                                MatJamovi.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1}).(Axes{iaxe})(inc,1) = iSubjects;
                                MatJamovi.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1}).(Axes{iaxe})(inc,2) = iTrial;
                                MatJamovi.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1}).(Axes{iaxe})(inc,3) = ...
                                    DataMat.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1}).(Axes{iaxe})(iSubjects, iTrial+1);
                            end
                            inc = inc +1;
                        end
                    end
                end
                
                %---Linear Mixed Model---
                for iSegment = 3:8%length(SegmentName)
                    if ipar == 4 | ipar == 5 | ipar == 6
                        tbl = table(MatJamovi.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1})(:,1),...
                            categorical(MatJamovi.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1})(:,2)),...
                            MatJamovi.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1})(:,3),...
                            'VariableNames',{'Subject','Trial','Measure'});
                    else
                        tbl = table(MatJamovi.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1}).(Axes{iaxe})(:,1),...
                            categorical(MatJamovi.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1}).(Axes{iaxe})(:,2)),...
                            MatJamovi.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1}).(Axes{iaxe})(:,3),...
                            'VariableNames',{'Subject','Trial','Measure'});
                    end
                    lme = fitlme(tbl,'Measure~Trial+(1|Subject)');
                    for i = 2:10
                        if double(lme.Coefficients(i,6)) < 0.05
                            Pvalues.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1}).(Axes{iaxe}) = double(lme.Coefficients(i,6));
                            break
                        else
                            Pvalues.(parametre{ipar}).(indicateur_time{iind_time}).(SegmentName{iSegment,1}).(Axes{iaxe}) = 1;
                        end
                    end
                end
                
                
            end
        end
    end
    
    cd(['\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted'])
    save(['MatJamovi_Time_', Task, '.mat'],'MatJamovi')
    save(['Pvalues_LM_Time_', Task, '.mat'],'Pvalues')
    
end








