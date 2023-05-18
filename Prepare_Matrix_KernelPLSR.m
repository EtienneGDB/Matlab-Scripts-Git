clear all;
close all;
clc;

Subjects = {'P01' 'P02' 'P03' 'P04' 'P05' 'P06' 'P07' 'P08' 'P09' 'P10' 'P11' 'P12' 'P13' 'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20' 'P21'...
    'P22' 'P23' 'P24' 'P25' 'P26' 'P27' 'P28' 'P29' 'P30' 'P31'};

parametre  = {'Acceleration','Angular_Velocity','Jerk','Module_Acceleration',...
    'Module_Angular_Velocity','Module_Jerk'};

indicateur = {'MedianFreq','SpectralEntropy','PowerLF','PowerHF','PowerTot', ...
    'PeakPower','PeakPower_Freq'};

indicateur_time = {'Peak', 'Mean', 'STD', 'Median', 'IQR', 'Percentile10',...
    'Percentile25', 'Percentile75', 'Percentile90'};

Seg_of_interest = {'Pelvis', 'T8', 'Head', 'RightShoulder', 'RightUpperArm',...
    'RightForeArm', 'RightHand', 'LeftShoulder', 'LeftUpperArm', 'LeftForeArm', 'LeftHand'};

Axes = {'x', 'y', 'z'};



Tasks =  {'SP', 'RPT', 'Work'};
for iTask = 1:length(Tasks)
    Task = Tasks{iTask};


    cd(['\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted'])
    spec_features = load(['MatJamovi_' Task '.mat']);
    time_features = load(['MatJamovi_Time_' Task '.mat']);

    % X_MatKernel
    VarNames = {};
    X_MatKernel = [];
    inc = 1;
    % spec_features.MatJamovi.Acceleration.MedianFreq.Pelvis.x(:,1:2);
    for ipar = 1:length(parametre)
        for iind = 1:length(indicateur)
            for iSeg = 1:length(Seg_of_interest)

                if ipar == 4 | ipar == 5 | ipar == 6
                    X_MatKernel(:,inc) = spec_features.MatJamovi.(parametre{ipar}).(indicateur{iind}).(Seg_of_interest{iSeg})(:,end);
                    VarNames{inc} = [(indicateur{iind}) ' ' (parametre{ipar}) ' ' (Seg_of_interest{iSeg})];
                    inc = inc + 1;

                else

                    for iaxe = 1:3
                        X_MatKernel(:,inc) = spec_features.MatJamovi.(parametre{ipar}).(indicateur{iind}).(Seg_of_interest{iSeg}).(Axes{iaxe})(:,end);
                        VarNames{inc} = [(indicateur{iind}) ' ' (parametre{ipar}) ' ' (Seg_of_interest{iSeg}) ' ' (Axes{iaxe})];
                        inc = inc + 1;
                    end
                    X_MatKernel(:,inc) = mean(X_MatKernel(:,inc-3:inc-1),2);
                    VarNames{inc} = [(indicateur{iind}) ' ' (parametre{ipar}) ' ' (Seg_of_interest{iSeg}) ' MeanXYZ'];
                    inc = inc + 1;

                end

            end
        end
    end

    for ipar = 1:length(parametre)
        for iind = 1:length(indicateur_time)
            for iSeg = 1:length(Seg_of_interest)

                if ipar == 4 | ipar == 5 | ipar == 6
                    X_MatKernel(:,inc) = time_features.MatJamovi.(parametre{ipar}).(indicateur_time{iind}).(Seg_of_interest{iSeg})(:,end);
                    VarNames{inc} = [(indicateur_time{iind}) ' ' (parametre{ipar}) ' ' (Seg_of_interest{iSeg})];
                    inc = inc + 1;

                else

                    for iaxe = 1:3
                        X_MatKernel(:,inc) = time_features.MatJamovi.(parametre{ipar}).(indicateur_time{iind}).(Seg_of_interest{iSeg}).(Axes{iaxe})(:,end);
                        VarNames{inc} = [(indicateur_time{iind}) ' ' (parametre{ipar}) ' ' (Seg_of_interest{iSeg}) ' ' (Axes{iaxe})];
                        inc = inc + 1;
                    end
                    X_MatKernel(:,inc) = mean(X_MatKernel(:,inc-3:inc-1),2);
                    VarNames{inc} = [(indicateur_time{iind}) ' ' (parametre{ipar}) ' ' (Seg_of_interest{iSeg}) ' MeanXYZ'];
                    inc = inc + 1;

                end

            end
        end
    end


    % Y_MatKernel
    Subjects = {'P1' 'P2' 'P3' 'P4' 'P5' 'P6' 'P7' 'P8' 'P9' 'P10' 'P11' 'P12' 'P13' 'P14' 'P15' 'P16' 'P17' 'P18' 'P19' 'P20' 'P21'...
        'P22' 'P23' 'P24' 'P25' 'P26' 'P27' 'P28' 'P29' 'P30' 'P31'};

    Variables = {'Effort', 'Arm pain', 'Shoulder pain', 'Neck pain', 'Cycle/s'} ;
    Y_MatKernel = [];
    inc = 1;
    for iS = 1:31
        CR100Table = [];
        cd('\\10.89.24.15\q\IRSST_DavidsTea\CR100')
        load(['CR100_' Subjects{iS} '.mat'])
        % Y_MatKernel = [Y_MatKernel; CR100Table(:,:)];
        Y_MatKernel = [Y_MatKernel; [CR100Table(:,1) mean(CR100Table(:,2:4),2)]];
    end

    cd('\\10.89.24.15\q\IRSST_DavidsTea')
    load('Longueur_15_cycles.mat')
    temp = [];
    for iS = 1:31
        temp = [temp; MatLengthWT(iS,:)'];
    end
    Y_MatKernel(:,size(Y_MatKernel,2)+1) = temp;

    Lign_Participants = spec_features.MatJamovi.Acceleration.MedianFreq.Pelvis.x(:,1);
    Lign_Trial = repmat([1:10]',31,1);

    Lign_Participants(isnan(X_MatKernel(:,1)),:) = [];
    Lign_Trial(isnan(X_MatKernel(:,1)),:) = [];

    Y_MatKernel(isnan(X_MatKernel(:,1)),:) = [];
    X_MatKernel(isnan(X_MatKernel(:,1)),:) = [];

    % Save MatKernel
    cd(['\\10.89.24.15\q\IRSST_DavidsTea\Data_extracted'])
    save(['X_MatKernel_' Task '.mat'],'X_MatKernel');
    save(['Y_MatKernel_' Task '.mat'],'Y_MatKernel');
    save('VarNames.mat',"VarNames");
    save(['Lign_Participants_' Task '.mat'],"Lign_Participants");
    save(['Lign_Trial_' Task '.mat'],"Lign_Trial");

end

RMSE_CrossVal = aaa;
VarExplained_LOO = aaa;
save(['RMSE_CrossVal_LOO_' Task '_pain.mat'],"RMSE_CrossVal");
save(['VarExplained_LOO_' Task '_pain.mat'],"VarExplained_LOO");






