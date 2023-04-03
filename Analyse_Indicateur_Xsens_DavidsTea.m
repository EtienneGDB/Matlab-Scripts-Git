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

%%---chose parameter before running---
Task =  'SP'; % RPT / Work / SP

ipar = 1; % parametre
iind = 7; % indicateur
coordonnees= 'x'; % x, y ou z 
%%------------------------------------

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

for iSegment = 1:length(SegmentName)
    for iSubject = 1:length(Mat_Load.(SegmentName{iSegment,1}))
        if isnan(Mat_Load.(SegmentName{iSegment,1})(iSubject,2))
            Mat_Load.(SegmentName{iSegment,1})(iSubject,2:4) = Mat_Load.(SegmentName{iSegment,1})(iSubject,3:5);
            Mat_Load.(SegmentName{iSegment,1})(iSubject,5) = NaN;
        end
    end
end

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
%             if (isnan(Norm.(SegmentName{iS,1})(iSuj,iTrial))) & (iTrial == 2)
%                 inc = iTrial + 1;
%                 while isnan(Norm.(SegmentName{iS,1})(iSuj,inc))
%                     if inc == 11
%                         break
%                     end
%                     inc = inc + 1;
%                 end
%                 FirstTrial = inc;
%             end

            Norm.(SegmentName{iS,1})(iSuj,iTrial) = Mat_Load.(SegmentName{iS,1})(iSuj,FirstTrial) - Norm.(SegmentName{iS,1})(iSuj,iTrial);

        end
    end
end
 vect = [1:2  4:6 8:13 15:26 28:31];
if ipar == 4 | ipar == 5 | ipar == 6
    % box plot of all the subject by Segment
    for iS = 1:length(SegmentName)
        % 2 figures
        if iS <10
            figure(1)
            sgtitle([parametre{ipar} ' ' indicateur{iind} ' ' Task])
            subplot(3,3,iS);
            plotver2 = Norm.(SegmentName{iS,1});
            plotver2 = plotver2(:,3:end);
            boxplot(plotver2)
            line([1 10],[0 0],'color','black')
            title(SegmentName{iS},'interpreter','none') 
        else
            figure(2)
            sgtitle([parametre{ipar} ' ' indicateur{iind} ' ' Task])
            subplot(3,3,iS-9);
            plotver2 = Norm.(SegmentName{iS,1});
            plotver2 = plotver2(:,3:end);
            boxplot(plotver2)
            line([1 10],[0 0],'color','black')
            title(SegmentName{iS},'interpreter','none') 
        end
    end
    
else

    % box plot of all the subject by Segment
    for iS = 1:length(SegmentName)
        % 2 figures
        if iS <10
            figure(1)
            sgtitle([parametre{ipar} ' ' coordonnees ' ' indicateur{iind} ' ' Task])
            subplot(3,3,iS);
            plotver2 = Norm.(SegmentName{iS,1});
            plotver2 = plotver2(start:start+30,3:end);
            boxplot(plotver2)
            line([1 10],[0 0],'color','black')
            title(SegmentName{iS},'interpreter','none') 
        else
            figure(2)
            sgtitle([parametre{ipar} ' ' coordonnees ' ' indicateur{iind} ' ' Task])
            subplot(3,3,iS-9);
            plotver2 = Norm.(SegmentName{iS,1});
            plotver2 = plotver2(start:start+30,3:end);
            boxplot(plotver2)
            line([1 10],[0 0],'color','black')
            title(SegmentName{iS},'interpreter','none') 
        end
    end
    
end


