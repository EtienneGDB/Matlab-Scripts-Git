% Save a matrix of one indicators, one parameter and one body segment in a structure 
% with all the body segment   

clear ;
close all ;
clc ;

addpath \\10.89.24.15\E\Bureau\Etienne\MATLAB\Functions
addpath(genpath('\\10.89.24.15\E\Bureau\Etienne\MATLAB\Functions\btk'))
addpath(genpath('\\10.89.24.15\E\Bureau\Etienne\MATLAB\Functions\wtc-r16'))
addpath(genpath('\\10.89.24.15\E\Bureau\Etienne\MATLAB\Functions\Matlab-Functions-Git'))
addpath(genpath('\\10.89.24.15\E\Bureau\Etienne\MATLAB\Functions\FFT'))
addpath(genpath('\\10.89.24.15\E\Bureau\Etienne\MATLAB\Functions\Entropy Fouaz'))

FoldersNames = dir('\\10.89.24.15\q\IRSST_DavidsTea\Raw_Data');
Subjects = {...
    FoldersNames(3:length(FoldersNames)).name... % May change if others files different from subject were add
   } ;

Task =  'SP';  % 'Work' / 'RPT' / 'SP'
tic
parametre  = {'Acceleration','Angular_Velocity','Jerk','Module_Acceleration',...
    'Module_Angular_Velocity','Module_Jerk'}; 

indicateur = {'MedianFreq','SpectralEntropy','PowerLF','PowerHF','PowerTot', ...
    'PeakPower','PeakPower_Freq'};

cd('\\10.89.24.15\q\IRSST_DavidsTea')
load('SegmentName_XSENS.mat')

for iSubjects = 1%:length(Subjects)
    Task = convertCharsToStrings(Task);     

    if Task == 'SP'  
        cd(['\\10.89.24.15\q\IRSST_DavidsTea\Data_exported\XSENS\Matlab\Filtered\Static_Pose'])
    elseif Task == 'Work'
        cd(['\\10.89.24.15\q\IRSST_DavidsTea\Data_exported\XSENS\Matlab\Filtered\Working_Task'])
    elseif Task == 'RPT'
        cd(['\\10.89.24.15\q\IRSST_DavidsTea\Data_exported\XSENS\Matlab\Filtered\RPT'])
    end
    
    %% Load files
    Task = convertStringsToChars(Task);
    load(['filtredXSENS_' Task '_' Subjects{iSubjects} '.mat'])
    TempMVNBL = eval(['MVNBL_' Task]);
    
    % To only keep the trial in the strucutre
    FN = fieldnames(TempMVNBL);
    NewFN = {};
    jFN = 1;
    for iFN = 1:length(FN)
        if strcmp(FN{iFN}(1), 'T') == 1
            NewFN(jFN,1) = FN(iFN);
            jFN = jFN+1;
        end
    end
    
    Features_XSENS = [];
    for iTrial = 1:length(NewFN)
        Trial = NewFN{iTrial};
        
            Freq = 60 ;
            dt=1/Freq ;
            t=(0:length(TempMVNBL.(Trial).AngularVelocity.Pelvis  )-1)*dt ;
            Local_Angular_Velocity = [];
            Local_Acceleration = [];
            if sum(TempMVNBL.(Trial).SegmentAcceleration.Pelvis(:,1))~=0 & nansum(TempMVNBL.(Trial).SegmentAcceleration.Pelvis(:,1))~=0
            % don't do it if the data are full of NaN or full of 0


            iCol = 1:3:51;
            for iSegment = 1:length(fieldnames(TempMVNBL.Trial1.SegmentAcceleration))
                Quat = TempMVNBL.(Trial).Orientation.(SegmentName{iSegment})  ;
                temp = LocalRef(TempMVNBL.(Trial).SegmentAcceleration.(SegmentName{iSegment})  ,Quat);
                Local_Acceleration.(SegmentName{iSegment}) = [temp];
                temp = LocalRef(TempMVNBL.(Trial).AngularVelocity.(SegmentName{iSegment}),Quat);
                Local_Angular_Velocity.(SegmentName{iSegment}) = [temp];

    %             figure; subplot(3,1,1); plot(Local_Acceleration(:,iCol(iSegment)))
    %             subplot(3,1,2); plot(Local_Acceleration(:,iCol(iSegment)+1))
    %             subplot(3,1,3); plot(Local_Acceleration(:,iCol(iSegment)+2))
    %             
    %             figure; subplot(3,1,1); plot(TempMVNBL.(Trial).SegmentAcceleration.(SegmentName{iSegment})(:,1))
    %             subplot(3,1,2); plot(TempMVNBL.(Trial).SegmentAcceleration.(SegmentName{iSegment})(:,2))
    %             subplot(3,1,3); plot(TempMVNBL.(Trial).SegmentAcceleration.(SegmentName{iSegment})(:,3))
    %             pause()
    %             close all

            end
            % Already filter in a previous script 
            Angular_Velocity_Filtered = Local_Angular_Velocity;
            Acceleration_Filtered = Local_Acceleration;

    %         [b,a] = butter(2,2*[0.5 29]/Freq) ; % Parametre du filtre BP 0.5-30 Hz
    %         Angular_Velocity_Filtered = filtfilt(b,a,Local_Angular_Velocity) ;
    %         Acceleration_Filtered = filtfilt(b,a,Local_Acceleration) ;

                % Jerk calculation
            Jerk = [];
            for iSegment = 1:length(fieldnames(TempMVNBL.Trial1.SegmentAcceleration))
                for iC = 1:size(Acceleration_Filtered.(SegmentName{iSegment}),2)
                     Jerk.(SegmentName{iSegment})(:,iC) = diff(Acceleration_Filtered.(SegmentName{iSegment})(:,iC))./diff(t');
                end
            end

    %         %     figure
    %         subplot(2,2,1); plot(Acceleration_Filtered(:,4))
    %         [FFT,f] = spec_fft(Acceleration_Filtered(:,4),Freq,0);
    %         subplot(2,2,2); plot(f',FFT)
    %         subplot(2,2,3); plot(Jerk(:,4))
    %         [FFT,f] = spec_fft(Jerk(:,4),Freq,0);
    %         subplot(2,2,4); plot(f',FFT)

            % Module
            Module_Angular_Velocity = [];
            Module_Acceleration = [];
            Module_Jerk = [];
            for iSegment = 1:length(fieldnames(TempMVNBL.Trial1.SegmentAcceleration))
                Module_Angular_Velocity.(SegmentName{iSegment}) = Module(Angular_Velocity_Filtered.(SegmentName{iSegment}));
                Module_Acceleration.(SegmentName{iSegment}) = Module(Acceleration_Filtered.(SegmentName{iSegment}));
                Module_Jerk.(SegmentName{iSegment}) = Module(Jerk.(SegmentName{iSegment}));
        %         subplot(4,1,1) ; plot(t,Module_Acceleration(:,VarInc))
        %         subplot(4,1,2) ; plot(t,Acceleration_Filtered(:,icol))
        %         subplot(4,1,3) ; plot(Acceleration_Filtered(:,icol+1))
        %         subplot(4,1,4) ; plot(Acceleration_Filtered(:,icol+2))
        %         pause()
            end

    %% EMG Median frequency
            f = 0.05:0.05:24.05;...Range de fréqence de la décomposition
            ond = ["cmor8-1"];
            scale = Freq*centfrq(ond)./f;

            % Acceleration
            disp([num2str(iSubjects) '.Acceleration'])
            for iSegment = 1:length(fieldnames(TempMVNBL.Trial1.SegmentAcceleration))
                for iData = 1:3
                    [norm, Period] = cwt(Acceleration_Filtered.(SegmentName{iSegment})(:,iData),scale,ond,dt);
                    Features_XSENS.(Trial).Acceleration.(SegmentName{iSegment}).MedianFreq(:,iData) = Compute_Median_Frequency(abs(norm),Period) ;
                    PSD = abs(norm).^2;
                    Features_XSENS.(Trial).Acceleration.(SegmentName{iSegment}).SpectralEntropy(:,iData) = Compute_Spectral_Entropy(PSD, Period) ;
                    [Features_XSENS.(Trial).Acceleration.(SegmentName{iSegment}).PowerLF(:,iData), Features_XSENS.(Trial).Acceleration.(SegmentName{iSegment}).PowerHF(:,iData), Features_XSENS.(Trial).Acceleration.(SegmentName{iSegment}).PowerTot(:,iData), Features_XSENS.(Trial).Acceleration.(SegmentName{iSegment}).PeakPower(:,iData), Features_XSENS.(Trial).Acceleration.(SegmentName{iSegment}).PeakPower_Freq(:,iData)] = Compute_Power(PSD, Period) ;
        %             pcolor(t,f,abs(norm)); shading interp
        %             pause()
                end
            end
            
            % Angular Velocity
            disp([num2str(iSubjects) '.Angular Velocity'])
            for iSegment = 1:length(fieldnames(TempMVNBL.Trial1.SegmentAcceleration))
                for iData = 1:3
                    [norm, Period] = cwt(Angular_Velocity_Filtered.(SegmentName{iSegment})(:,iData),scale,ond,dt);
                    Features_XSENS.(Trial).Angular_Velocity.(SegmentName{iSegment}).MedianFreq(:,iData) = Compute_Median_Frequency(abs(norm),Period) ;
                    PSD = abs(norm).^2;
                    Features_XSENS.(Trial).Angular_Velocity.(SegmentName{iSegment}).SpectralEntropy(:,iData) = Compute_Spectral_Entropy(PSD, Period) ;
                    [Features_XSENS.(Trial).Angular_Velocity.(SegmentName{iSegment}).PowerLF(:,iData), Features_XSENS.(Trial).Angular_Velocity.(SegmentName{iSegment}).PowerHF(:,iData), Features_XSENS.(Trial).Angular_Velocity.(SegmentName{iSegment}).PowerTot(:,iData), Features_XSENS.(Trial).Angular_Velocity.(SegmentName{iSegment}).PeakPower(:,iData), Features_XSENS.(Trial).Angular_Velocity.(SegmentName{iSegment}).PeakPower_Freq(:,iData)] = Compute_Power(PSD, Period) ;
        %             pcolor(t,f,abs(norm)); shading interp
        %             pause()
                end
            end
            
            % Jerk
            disp([num2str(iSubjects) '.Jerk'])
            for iSegment = 1:length(fieldnames(TempMVNBL.Trial1.SegmentAcceleration))
                for iData = 1:3
                    [norm, Period] = cwt(Jerk.(SegmentName{iSegment})(:,iData),scale,ond,dt);
                    Features_XSENS.(Trial).Jerk.(SegmentName{iSegment}).MedianFreq(:,iData) = Compute_Median_Frequency(abs(norm),Period) ;
                    PSD = abs(norm).^2;
                    Features_XSENS.(Trial).Jerk.(SegmentName{iSegment}).SpectralEntropy(:,iData) = Compute_Spectral_Entropy(PSD, Period) ;
                    [Features_XSENS.(Trial).Jerk.(SegmentName{iSegment}).PowerLF(:,iData), Features_XSENS.(Trial).Jerk.(SegmentName{iSegment}).PowerHF(:,iData), Features_XSENS.(Trial).Jerk.(SegmentName{iSegment}).PowerTot(:,iData), Features_XSENS.(Trial).Jerk.(SegmentName{iSegment}).PeakPower(:,iData), Features_XSENS.(Trial).Jerk.(SegmentName{iSegment}).PeakPower_Freq(:,iData)] = Compute_Power(PSD, Period) ;
        %             pcolor(t(1:length(t)-1),f,abs(norm)); shading interp
        %             pause()
                end
            end

            % Module Acceleration
            disp([num2str(iSubjects) '.Module Acceleration'])
            for iSegment = 1:length(fieldnames(TempMVNBL.Trial1.SegmentAcceleration))
                for iData = 1
                    [norm, Period] = cwt(Module_Acceleration.(SegmentName{iSegment})(:,iData),scale,ond,dt);
                    Features_XSENS.(Trial).Module_Acceleration.(SegmentName{iSegment}).MedianFreq(:,iData) = Compute_Median_Frequency(abs(norm),Period) ;
                    PSD = abs(norm).^2;
                    Features_XSENS.(Trial).Module_Acceleration.(SegmentName{iSegment}).SpectralEntropy(:,iData) = Compute_Spectral_Entropy(PSD, Period) ;
                    [Features_XSENS.(Trial).Module_Acceleration.(SegmentName{iSegment}).PowerLF(:,iData), Features_XSENS.(Trial).Module_Acceleration.(SegmentName{iSegment}).PowerHF(:,iData), Features_XSENS.(Trial).Module_Acceleration.(SegmentName{iSegment}).PowerTot(:,iData), Features_XSENS.(Trial).Module_Acceleration.(SegmentName{iSegment}).PeakPower(:,iData), Features_XSENS.(Trial).Module_Acceleration.(SegmentName{iSegment}).PeakPower_Freq(:,iData)] = Compute_Power(PSD, Period) ;
        %             pcolor(t,f,abs(norm)); shading interp
        %             pause()
                end
            end
            
            % Module Angular Velocity
            disp([num2str(iSubjects) '.Module Angular Velocity'])
            for iSegment = 1:length(fieldnames(TempMVNBL.Trial1.SegmentAcceleration))
                for iData = 1
                    [norm, Period] = cwt(Module_Angular_Velocity.(SegmentName{iSegment})(:,iData),scale,ond,dt);
                    Features_XSENS.(Trial).Module_Angular_Velocity.(SegmentName{iSegment}).MedianFreq(:,iData) = Compute_Median_Frequency(abs(norm),Period) ;
                    PSD = abs(norm).^2;
                    Features_XSENS.(Trial).Module_Angular_Velocity.(SegmentName{iSegment}).SpectralEntropy(:,iData) = Compute_Spectral_Entropy(PSD, Period) ;
                    [Features_XSENS.(Trial).Module_Angular_Velocity.(SegmentName{iSegment}).PowerLF(:,iData), Features_XSENS.(Trial).Module_Angular_Velocity.(SegmentName{iSegment}).PowerHF(:,iData), Features_XSENS.(Trial).Module_Angular_Velocity.(SegmentName{iSegment}).PowerTot(:,iData), Features_XSENS.(Trial).Module_Angular_Velocity.(SegmentName{iSegment}).PeakPower(:,iData), Features_XSENS.(Trial).Module_Angular_Velocity.(SegmentName{iSegment}).PeakPower_Freq(:,iData)] = Compute_Power(PSD, Period) ;
        %             pcolor(t,f,abs(norm)); shading interp
        %             pause()
                end
            end
          
            % Module Jerk
            disp([num2str(iSubjects) '.Module Jerk'])
            for iSegment = 1:length(fieldnames(TempMVNBL.Trial1.SegmentAcceleration))
                for iData = 1
                    [norm, Period] = cwt(Module_Jerk.(SegmentName{iSegment})(:,iData),scale,ond,dt);
                    Features_XSENS.(Trial).Module_Jerk.(SegmentName{iSegment}).MedianFreq(:,iData) = Compute_Median_Frequency(abs(norm),Period) ;
                    PSD = abs(norm).^2;
                    Features_XSENS.(Trial).Module_Jerk.(SegmentName{iSegment}).SpectralEntropy(:,iData) = Compute_Spectral_Entropy(PSD, Period) ;
                    [Features_XSENS.(Trial).Module_Jerk.(SegmentName{iSegment}).PowerLF(:,iData), Features_XSENS.(Trial).Module_Jerk.(SegmentName{iSegment}).PowerHF(:,iData), Features_XSENS.(Trial).Module_Jerk.(SegmentName{iSegment}).PowerTot(:,iData), Features_XSENS.(Trial).Module_Jerk.(SegmentName{iSegment}).PeakPower(:,iData), Features_XSENS.(Trial).Module_Jerk.(SegmentName{iSegment}).PeakPower_Freq(:,iData)] = Compute_Power(PSD, Period) ;
        %             pcolor(t,f,abs(norm)); shading interp
        %             pause()
                end
            end
            
        else                
                
                for iind = 1:length(parametre)
                    for iind2 = 1:length(indicateur)
                        cpt = 1;
                            if iind > 3
                                Features_XSENS.(Trial).(parametre{iind}).(indicateur{iind2}) =   NaN(500,17) ;                        
                                cpt = cpt+1;
                                continue
                            end
                            Features_XSENS.(Trial).(parametre{iind}).(indicateur{iind2}) =   NaN(500,51) ;
                            cpt = cpt+1;
                    end
                end
        end
        
  
     
    end
    save(['TFR_' Task '_' (Subjects{iSubjects}) '.mat'], 'Features_XSENS')
end



