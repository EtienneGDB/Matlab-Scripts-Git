clear ;
close all ;
clc ;

% addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
% addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
% addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
% addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
% addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
% addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Entropy Fouaz'))

addpath E:\Bureau\Etienne\MATLAB\Functions
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\btk'))
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\Matlab-Functions-Git'))
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\FFT'))
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\Entropy Fouaz'))

FoldersNames = dir('V:\Piano_Fatigue\Data_Compiled');
Filenames = {...
    FoldersNames(4:end).name...
   } ;

%Crée une deuxième ligne avec le nom du fichier sans .mat
Filenames{2,1}=[];
for iFilenames = 1:length(Filenames)
    Filenames{2,iFilenames} = erase(Filenames{1,iFilenames},'.mat');
end
    
cd(['V:\Piano_Fatigue\Data_Exported'])
load('cycles_Li.mat');  % Charge les cycle de base de Valentin pour garder juste le 1er et dernier

PartInc = 18;
for iSubjects = 36:2:100
    % Load XSENS data
    cd(['J:\Piano_Fatigue\Data_Compiled'])
    load([Filenames{1,iSubjects}]);
    Angular_Velocity = table2array(data.angularVelocity(:,[4:6 13:15 19:21 22:33])) ; % [13:15 19:21] pour L5, T8, head... & 22:33 pour Rshoulder à Rhand
    Acceleration = table2array(data.acceleration(:,[4:6 13:15 19:21 22:33])) ;
    Quaternion = table2array(data.orientation(:,[5:8 17:20 25:44])) ;
    Freq = 60 ;
    dt=1/Freq ;
    t=(0:length(Angular_Velocity)-1)*dt ;
    cycles_Valentin = participants(PartInc).t_do;
    
    Local_Angular_Velocity = [];
    Local_Acceleration = [];
    iCol = 1:3:size(Angular_Velocity,2);
    iCol_Quat = 1:4:size(Quaternion,2);
    for iSegment = 1:length(iCol_Quat)
        Quat = Quaternion(:,iCol_Quat(iSegment):iCol_Quat(iSegment)+3);
        temp = LocalRef(Acceleration(:,iCol(iSegment):iCol(iSegment)+2),Quat);
        Local_Acceleration(:,iCol(iSegment):iCol(iSegment)+2) = temp;
        temp = LocalRef(Angular_Velocity(:,iCol(iSegment):iCol(iSegment)+2),Quat);
        Local_Angular_Velocity(:,iCol(iSegment):iCol(iSegment)+2) = temp;
                
%         figure; subplot(3,1,1); plot(t,Acceleration(:,iCol(iSegment)))
%         subplot(3,1,2); plot(Acceleration(:,iCol(iSegment)+1))
%         subplot(3,1,3); plot(Acceleration(:,iCol(iSegment)+2))
% 
%         figure; subplot(3,1,1); plot(Local_Acceleration(:,iCol(iSegment)))
%         subplot(3,1,2); plot(Local_Acceleration(:,iCol(iSegment)+1))
%         subplot(3,1,3); plot(Local_Acceleration(:,iCol(iSegment)+2))
%         pause()
%         close all
    end
    
        % Preprocessing
        % Freq : fréquence d'échantillonage de l'EMG
        [b,a] = butter(2,2*[0.5 29]/Freq) ; % Parametre du filtre BP 0.5-30 Hz
        Angular_Velocity_Filtered = filtfilt(b,a,Local_Angular_Velocity) ;
        Acceleration_Filtered = filtfilt(b,a,Local_Acceleration) ;
%         cycles = round(data.cycles);
%         for iM = 1:12
%             plot([1:length(Acceleration_Filtered(:,iM))]/60,Acceleration_Filtered(:,iM))
%             hold on
%             line([cycles_Valentin(:,1),cycles_Valentin(:,1)]', repmat(ylim,length(cycles_Valentin(:,1)),1)','color','red')
%             line([cycles(:,1),cycles(:,1)]', repmat(ylim,length(cycles(:,1)),1)','color','red')
%             line([cycles(:,end),cycles(:,end)]', repmat(ylim,length(cycles(:,end)),1)','color','red')
%             spec_fft(Acceleration_Filtered(:,iM),Freq,1);
%             pause()
%         end

    % Jerk calculation
    Jerk = [];
    for iC = 1:size(Acceleration_Filtered,2)
        Jerk(:,iC) = diff(Acceleration_Filtered(:,iC))./diff(t');
    end
%     figure
%     subplot(2,2,1); plot(Acceleration_Filtered(:,1))
%     [FFT,f] = spec_fft(Acceleration_Filtered(:,1),Freq,0);
%     subplot(2,2,2); plot(f',FFT)
%     subplot(2,2,3); plot(Jerk(:,1))
%     [FFT,f] = spec_fft(Jerk(:,1),Freq,0);
%     subplot(2,2,4); plot(f',FFT)

    % Module
    Module_Angular_Velocity = [];
    Module_Acceleration = [];
    Module_Jerk = [];
    VarInc = 1;
    for icol = 1:3:size(Angular_Velocity,2)
        Module_Angular_Velocity(:,VarInc) = Module(Angular_Velocity_Filtered(:,icol:icol+2));
        Module_Acceleration(:,VarInc) = Module(Acceleration_Filtered(:,icol:icol+2));
        Module_Jerk(:,VarInc) = Module(Jerk(:,icol:icol+2));
%         subplot(4,1,1) ; plot(t,Module_Acceleration(:,VarInc))
%         subplot(4,1,2) ; plot(t,Acceleration_Filtered(:,icol))
%         subplot(4,1,3) ; plot(Acceleration_Filtered(:,icol+1))
%         subplot(4,1,4) ; plot(Acceleration_Filtered(:,icol+2))
%         pause()
        VarInc = VarInc+1 ;
    end

% %% EMG Median frequency
%         f = 0.05:0.05:24.05;...Range de fréqence de la décomposition
%         ond = ["cmor8-1"];
%         scale = Freq*centfrq(ond)./f;
%         Features_XSENS = [];
%         % Angular Velocity
%         disp([num2str(iSubjects) '.Angular Velocity'])
%         for iData = 1:size(Angular_Velocity,2)
%             [norm, Period] = cwt(Angular_Velocity_Filtered(:,iData),scale,ond,dt);
%             Features_XSENS.MedianFreq.Angular_Velocity(:,iData) = Compute_Median_Frequency(abs(norm),Period) ;
%             Features_XSENS.SpectralEntropy.Angular_Velocity(:,iData) = Compute_Spectral_Entropy(abs(norm), Period) ;
%             [Features_XSENS.PowerLF.Angular_Velocity(:,iData), Features_XSENS.PowerHF.Angular_Velocity(:,iData), Features_XSENS.PowerTot.Angular_Velocity(:,iData), Features_XSENS.PeakPower.Angular_Velocity(:,iData), Features_XSENS.PeakPower_Freq.Angular_Velocity(:,iData)] = Compute_Power(abs(norm), Period) ;
% %             pcolor(t,f,abs(norm)); shading interp
% %             pause()
%         end
%         
%         % Acceleration
%         disp([num2str(iSubjects) '.Acceleration'])
%         for iData = 1:size(Acceleration,2)
%             [norm, Period] = cwt(Acceleration_Filtered(:,iData),scale,ond,dt);
%             Features_XSENS.MedianFreq.Acceleration(:,iData) = Compute_Median_Frequency(abs(norm),Period) ;
%             Features_XSENS.SpectralEntropy.Acceleration(:,iData) = Compute_Spectral_Entropy(abs(norm), Period) ;
%             [Features_XSENS.PowerLF.Acceleration(:,iData), Features_XSENS.PowerHF.Acceleration(:,iData), Features_XSENS.PowerTot.Acceleration(:,iData), Features_XSENS.PeakPower.Acceleration(:,iData), Features_XSENS.PeakPower_Freq.Acceleration(:,iData)] = Compute_Power(abs(norm), Period) ;
% %             pcolor(t,f,abs(norm)); shading interp
% %             pause()
%         end
%         
%         % Jerk
%         disp([num2str(iSubjects) '.Jerk'])
%         for iData = 1:size(Jerk,2)
%             [norm, Period] = cwt(Jerk(:,iData),scale,ond,dt);
%             Features_XSENS.MedianFreq.Jerk(:,iData) = Compute_Median_Frequency(abs(norm),Period) ;
%             Features_XSENS.SpectralEntropy.Jerk(:,iData) = Compute_Spectral_Entropy(abs(norm), Period) ;
%             [Features_XSENS.PowerLF.Jerk(:,iData), Features_XSENS.PowerHF.Jerk(:,iData), Features_XSENS.PowerTot.Jerk(:,iData), Features_XSENS.PeakPower.Jerk(:,iData), Features_XSENS.PeakPower_Freq.Jerk(:,iData)] = Compute_Power(abs(norm), Period) ;
% %             pcolor(t(1:length(t)-1),f,abs(norm)); shading interp
% %             pause()
%         end
%         
%         % Module Angular Velocity
%         disp([num2str(iSubjects) '.Module Angular Velocity'])
%         for iData = 1:size(Module_Angular_Velocity,2)
%             [norm, Period] = cwt(Module_Angular_Velocity(:,iData),scale,ond,dt);
%             Features_XSENS.MedianFreq.Module_Angular_Velocity(:,iData) = Compute_Median_Frequency(abs(norm),Period) ;
%             Features_XSENS.SpectralEntropy.Module_Angular_Velocity(:,iData) = Compute_Spectral_Entropy(abs(norm), Period) ;
%             [Features_XSENS.PowerLF.Module_Angular_Velocity(:,iData), Features_XSENS.PowerHF.Module_Angular_Velocity(:,iData), Features_XSENS.PowerTot.Module_Angular_Velocity(:,iData), Features_XSENS.PeakPower.Module_Angular_Velocity(:,iData), Features_XSENS.PeakPower_Freq.Module_Angular_Velocity(:,iData)] = Compute_Power(abs(norm), Period) ;
% %             pcolor(t,f,abs(norm)); shading interp
% %             pause()
%         end
%         
%         % Module Acceleration
%         disp([num2str(iSubjects) '.Module Acceleration'])
%         for iData = 1:size(Module_Acceleration,2)
%             [norm, Period] = cwt(Module_Acceleration(:,iData),scale,ond,dt);
%             Features_XSENS.MedianFreq.Module_Acceleration(:,iData) = Compute_Median_Frequency(abs(norm),Period) ;
%             Features_XSENS.SpectralEntropy.Module_Acceleration(:,iData) = Compute_Spectral_Entropy(abs(norm), Period) ;
%             [Features_XSENS.PowerLF.Module_Acceleration(:,iData), Features_XSENS.PowerHF.Module_Acceleration(:,iData), Features_XSENS.PowerTot.Module_Acceleration(:,iData), Features_XSENS.PeakPower.Module_Acceleration(:,iData), Features_XSENS.PeakPower_Freq.Module_Acceleration(:,iData)] = Compute_Power(abs(norm), Period) ;
% %             pcolor(t,f,abs(norm)); shading interp
% %             pause()
%         end
%         
%         % Module Jerk
%         disp([num2str(iSubjects) '.Module Jerk'])
%         for iData = 1:size(Module_Jerk,2)
%             [norm, Period] = cwt(Module_Jerk(:,iData),scale,ond,dt);
%             Features_XSENS.MedianFreq.Module_Jerk(:,iData) = Compute_Median_Frequency(abs(norm),Period) ;
%             Features_XSENS.SpectralEntropy.Module_Jerk(:,iData) = Compute_Spectral_Entropy(abs(norm), Period) ;
%             [Features_XSENS.PowerLF.Module_Jerk(:,iData), Features_XSENS.PowerHF.Module_Jerk(:,iData), Features_XSENS.PowerTot.Module_Jerk(:,iData), Features_XSENS.PeakPower.Module_Jerk(:,iData), Features_XSENS.PeakPower_Freq.Module_Jerk(:,iData)] = Compute_Power(abs(norm), Period) ;
% %             pcolor(t(1:length(t)-1),f,abs(norm)); shading interp
% %             pause()
%         end
%         
%         % Features par Cycles
%         disp([num2str(iSubjects) '.Features par Cycles'])
%         cycles = round(data.cycles*60);
% %             plot(Acceleration_Filtered(:,1))
% %             line([cycles(:,1),cycles(:,1)]', repmat(ylim,length(cycles(:,1)),1)','color','red')
% %             line([cycles(:,end),cycles(:,end)]', repmat(ylim,length(cycles(:,end)),1)','color','red')
% %             pause()
%         Features_XSENS_Cycles = [];
%         for iCycles = 1:size(cycles,1)-1
%             Features_XSENS_Cycles.MedianFreq.Angular_Velocity(iCycles,:) = mean(Features_XSENS.MedianFreq.Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.MedianFreq.Acceleration(iCycles,:) = mean(Features_XSENS.MedianFreq.Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.MedianFreq.Jerk(iCycles,:) = mean(Features_XSENS.MedianFreq.Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.MedianFreq.Module_Angular_Velocity(iCycles,:) = mean(Features_XSENS.MedianFreq.Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.MedianFreq.Module_Acceleration(iCycles,:) = mean(Features_XSENS.MedianFreq.Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.MedianFreq.Module_Jerk(iCycles,:) = mean(Features_XSENS.MedianFreq.Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
% 
%             Features_XSENS_Cycles.SpectralEntropy.Angular_Velocity(iCycles,:) = mean(Features_XSENS.SpectralEntropy.Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.SpectralEntropy.Acceleration(iCycles,:) = mean(Features_XSENS.SpectralEntropy.Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.SpectralEntropy.Jerk(iCycles,:) = mean(Features_XSENS.SpectralEntropy.Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.SpectralEntropy.Module_Angular_Velocity(iCycles,:) = mean(Features_XSENS.SpectralEntropy.Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.SpectralEntropy.Module_Acceleration(iCycles,:) = mean(Features_XSENS.SpectralEntropy.Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.SpectralEntropy.Module_Jerk(iCycles,:) = mean(Features_XSENS.SpectralEntropy.Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
% 
%             Features_XSENS_Cycles.PowerLF.Angular_Velocity(iCycles,:) = mean(Features_XSENS.PowerLF.Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerLF.Acceleration(iCycles,:) = mean(Features_XSENS.PowerLF.Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerLF.Jerk(iCycles,:) = mean(Features_XSENS.PowerLF.Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerLF.Module_Angular_Velocity(iCycles,:) = mean(Features_XSENS.PowerLF.Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerLF.Module_Acceleration(iCycles,:) = mean(Features_XSENS.PowerLF.Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerLF.Module_Jerk(iCycles,:) = mean(Features_XSENS.PowerLF.Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
% 
%             Features_XSENS_Cycles.PowerHF.Angular_Velocity(iCycles,:) = mean(Features_XSENS.PowerHF.Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerHF.Acceleration(iCycles,:) = mean(Features_XSENS.PowerHF.Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerHF.Jerk(iCycles,:) = mean(Features_XSENS.PowerHF.Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerHF.Module_Angular_Velocity(iCycles,:) = mean(Features_XSENS.PowerHF.Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerHF.Module_Acceleration(iCycles,:) = mean(Features_XSENS.PowerHF.Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerHF.Module_Jerk(iCycles,:) = mean(Features_XSENS.PowerHF.Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
% 
%             Features_XSENS_Cycles.PowerTot.Angular_Velocity(iCycles,:) = mean(Features_XSENS.PowerTot.Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerTot.Acceleration(iCycles,:) = mean(Features_XSENS.PowerTot.Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerTot.Jerk(iCycles,:) = mean(Features_XSENS.PowerTot.Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerTot.Module_Angular_Velocity(iCycles,:) = mean(Features_XSENS.PowerTot.Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerTot.Module_Acceleration(iCycles,:) = mean(Features_XSENS.PowerTot.Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PowerTot.Module_Jerk(iCycles,:) = mean(Features_XSENS.PowerTot.Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
% 
%             Features_XSENS_Cycles.PeakPower.Angular_Velocity(iCycles,:) = mean(Features_XSENS.PeakPower.Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PeakPower.Acceleration(iCycles,:) = mean(Features_XSENS.PeakPower.Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PeakPower.Jerk(iCycles,:) = mean(Features_XSENS.PeakPower.Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PeakPower.Module_Angular_Velocity(iCycles,:) = mean(Features_XSENS.PeakPower.Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PeakPower.Module_Acceleration(iCycles,:) = mean(Features_XSENS.PeakPower.Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PeakPower.Module_Jerk(iCycles,:) = mean(Features_XSENS.PeakPower.Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
% 
%             Features_XSENS_Cycles.PeakPower_Freq.Angular_Velocity(iCycles,:) = mean(Features_XSENS.PeakPower_Freq.Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PeakPower_Freq.Acceleration(iCycles,:) = mean(Features_XSENS.PeakPower_Freq.Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PeakPower_Freq.Jerk(iCycles,:) = mean(Features_XSENS.PeakPower_Freq.Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PeakPower_Freq.Module_Angular_Velocity(iCycles,:) = mean(Features_XSENS.PeakPower_Freq.Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PeakPower_Freq.Module_Acceleration(iCycles,:) = mean(Features_XSENS.PeakPower_Freq.Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.PeakPower_Freq.Module_Jerk(iCycles,:) = mean(Features_XSENS.PeakPower_Freq.Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
% 
%             Features_XSENS_Cycles.Peak.Angular_Velocity(iCycles,:) = max(Angular_Velocity_Filtered(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Peak.Acceleration(iCycles,:) = max(Acceleration_Filtered(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Peak.Jerk(iCycles,:) = max(Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Peak.Module_Angular_Velocity(iCycles,:) = max(Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Peak.Module_Acceleration(iCycles,:) = max(Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Peak.Module_Jerk(iCycles,:) = max(Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
% 
%             Features_XSENS_Cycles.Mean.Angular_Velocity(iCycles,:) = mean(Angular_Velocity_Filtered(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Mean.Acceleration(iCycles,:) = mean(Acceleration_Filtered(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Mean.Jerk(iCycles,:) = mean(Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Mean.Module_Angular_Velocity(iCycles,:) = mean(Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Mean.Module_Acceleration(iCycles,:) = mean(Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Mean.Module_Jerk(iCycles,:) = mean(Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
% 
%             Features_XSENS_Cycles.Median.Angular_Velocity(iCycles,:) = median(Angular_Velocity_Filtered(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Median.Acceleration(iCycles,:) = median(Acceleration_Filtered(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Median.Jerk(iCycles,:) = median(Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Median.Module_Angular_Velocity(iCycles,:) = median(Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Median.Module_Acceleration(iCycles,:) = median(Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.Median.Module_Jerk(iCycles,:) = median(Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
% 
%             Features_XSENS_Cycles.STD.Angular_Velocity(iCycles,:) = std(Angular_Velocity_Filtered(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.STD.Acceleration(iCycles,:) = std(Acceleration_Filtered(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.STD.Jerk(iCycles,:) = std(Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.STD.Module_Angular_Velocity(iCycles,:) = std(Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.STD.Module_Acceleration(iCycles,:) = std(Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:));
%             Features_XSENS_Cycles.STD.Module_Jerk(iCycles,:) = std(Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:));
% 
%             Features_XSENS_Cycles.Percentile10.Angular_Velocity(iCycles,:) = prctile(Angular_Velocity_Filtered(cycles(iCycles,1):cycles(iCycles,end),:),10);
%             Features_XSENS_Cycles.Percentile10.Acceleration(iCycles,:) = prctile(Acceleration_Filtered(cycles(iCycles,1):cycles(iCycles,end),:),10);
%             Features_XSENS_Cycles.Percentile10.Jerk(iCycles,:) = prctile(Jerk(cycles(iCycles,1):cycles(iCycles,end),:),10);
%             Features_XSENS_Cycles.Percentile10.Module_Angular_Velocity(iCycles,:) = prctile(Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:),10);
%             Features_XSENS_Cycles.Percentile10.Module_Acceleration(iCycles,:) = prctile(Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:),10);
%             Features_XSENS_Cycles.Percentile10.Module_Jerk(iCycles,:) = prctile(Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:),10);
% 
%             Features_XSENS_Cycles.Percentile25.Angular_Velocity(iCycles,:) = prctile(Angular_Velocity_Filtered(cycles(iCycles,1):cycles(iCycles,end),:),25);
%             Features_XSENS_Cycles.Percentile25.Acceleration(iCycles,:) = prctile(Acceleration_Filtered(cycles(iCycles,1):cycles(iCycles,end),:),25);
%             Features_XSENS_Cycles.Percentile25.Jerk(iCycles,:) = prctile(Jerk(cycles(iCycles,1):cycles(iCycles,end),:),25);
%             Features_XSENS_Cycles.Percentile25.Module_Angular_Velocity(iCycles,:) = prctile(Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:),25);
%             Features_XSENS_Cycles.Percentile25.Module_Acceleration(iCycles,:) = prctile(Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:),25);
%             Features_XSENS_Cycles.Percentile25.Module_Jerk(iCycles,:) = prctile(Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:),25);
% 
%             Features_XSENS_Cycles.Percentile75.Angular_Velocity(iCycles,:) = prctile(Angular_Velocity_Filtered(cycles(iCycles,1):cycles(iCycles,end),:),75);
%             Features_XSENS_Cycles.Percentile75.Acceleration(iCycles,:) = prctile(Acceleration_Filtered(cycles(iCycles,1):cycles(iCycles,end),:),75);
%             Features_XSENS_Cycles.Percentile75.Jerk(iCycles,:) = prctile(Jerk(cycles(iCycles,1):cycles(iCycles,end),:),75);
%             Features_XSENS_Cycles.Percentile75.Module_Angular_Velocity(iCycles,:) = prctile(Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:),75);
%             Features_XSENS_Cycles.Percentile75.Module_Acceleration(iCycles,:) = prctile(Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:),75);
%             Features_XSENS_Cycles.Percentile75.Module_Jerk(iCycles,:) = prctile(Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:),75);
% 
%             Features_XSENS_Cycles.Percentile90.Angular_Velocity(iCycles,:) = prctile(Angular_Velocity_Filtered(cycles(iCycles,1):cycles(iCycles,end),:),90);
%             Features_XSENS_Cycles.Percentile90.Acceleration(iCycles,:) = prctile(Acceleration_Filtered(cycles(iCycles,1):cycles(iCycles,end),:),90);
%             Features_XSENS_Cycles.Percentile90.Jerk(iCycles,:) = prctile(Jerk(cycles(iCycles,1):cycles(iCycles,end),:),90);
%             Features_XSENS_Cycles.Percentile90.Module_Angular_Velocity(iCycles,:) = prctile(Module_Angular_Velocity(cycles(iCycles,1):cycles(iCycles,end),:),90);
%             Features_XSENS_Cycles.Percentile90.Module_Acceleration(iCycles,:) = prctile(Module_Acceleration(cycles(iCycles,1):cycles(iCycles,end),:),90);
%             Features_XSENS_Cycles.Percentile90.Module_Jerk(iCycles,:) = prctile(Module_Jerk(cycles(iCycles,1):cycles(iCycles,end),:),90);
%         end
%         
%         % Correct Features_XSENS en supprimant le début et la fin selon les cycles de Valentin pour
%         % garder uniquement les données correspondant à la tâche
%         porte = 1*60; % porte de 1sec
%         Indicators = fields(Features_XSENS);
%         for iIndicators = 1:length(fields(Features_XSENS))
%             Variables = fields(Features_XSENS.(Indicators{iIndicators}));
%             for iVariables = 1:length(fields(Features_XSENS.(Indicators{iIndicators})))
%                 O = zeros(porte,size(Features_XSENS.(Indicators{iIndicators}).(Variables{iVariables}),2));
%                 inc = 1;
%                 temp = [];
%                 for iData = round(cycles_Valentin(1)*60):length(O):round(cycles_Valentin(end)*60)
%                     temp(inc,:) = mean(Features_XSENS.(Indicators{iIndicators}).(Variables{iVariables})...
%                         (iData-length(O)/2:iData+length(O)/2,:));
%                     inc = inc + 1;
%                 end
%                 Features_XSENS.(Indicators{iIndicators}).(Variables{iVariables}) = temp;
%             end
%         end

%% Correction erreur
Task = 'Li'; % changer Li2 et Ha2 partout ('Matrix dimension must agree')
% FoldersNames = dir(['J:\Piano_Fatigue\Data_Exported\Features_XSENS_Cycles\' Task]);
FoldersNames = dir(['V:\Piano_Fatigue\Data_Exported\Features_XSENS\' Task]);
FilenamesXSENS = {...
    FoldersNames(3:end).name...
   } ;

%Crée une deuxième ligne avec le nom du fichier sans .mat
FilenamesXSENS{2,1}=[];
for iFilenames = 1:length(FilenamesXSENS)
    FilenamesXSENS{2,iFilenames} = erase(FilenamesXSENS{1,iFilenames},'.mat');
end
        
    cd(['V:\Piano_Fatigue\Data_Exported\Features_XSENS\' Task])
    load([FilenamesXSENS{1,PartInc}])
    % Run temporal features

        %% Temporal features
        porte = 1*60; % porte de 1sec
        % Angular Velocity
        disp([num2str(iSubjects) '.Angular Velocity'])
        O = zeros(porte,size(Angular_Velocity,2));
        inc = 1;
        for iData = round(cycles_Valentin(1)*60):length(O):round(cycles_Valentin(end)*60)
            Features_XSENS.Peak.Angular_Velocity(inc,:) = max(Angular_Velocity_Filtered(iData-length(O)/2:iData+length(O)/2,:));
            Features_XSENS.Mean.Angular_Velocity(inc,:) = mean(Angular_Velocity_Filtered(iData-length(O)/2:iData+length(O)/2,:));
            Features_XSENS.Median.Angular_Velocity(inc,:) = median(Angular_Velocity_Filtered(iData-length(O)/2:iData+length(O)/2,:));
            Features_XSENS.STD.Angular_Velocity(inc,:) = std(Angular_Velocity_Filtered(iData-length(O)/2:iData+length(O)/2,:));
            Features_XSENS.Percentile10.Angular_Velocity(inc,:) = prctile(Angular_Velocity_Filtered(iData-length(O)/2:iData+length(O)/2,:),10);
            Features_XSENS.Percentile25.Angular_Velocity(inc,:) = prctile(Angular_Velocity_Filtered(iData-length(O)/2:iData+length(O)/2,:),25);
            Features_XSENS.Percentile75.Angular_Velocity(inc,:) = prctile(Angular_Velocity_Filtered(iData-length(O)/2:iData+length(O)/2,:),75);
            Features_XSENS.Percentile90.Angular_Velocity(inc,:) = prctile(Angular_Velocity_Filtered(iData-length(O)/2:iData+length(O)/2,:),90);
            inc = inc + 1;
        end
        
        % Acceleration
        disp([num2str(iSubjects) '.Acceleration'])
        O = zeros(porte,size(Acceleration,2));
        inc = 1;
        for iData = round(cycles_Valentin(1)*60):length(O):round(cycles_Valentin(end)*60)
            Features_XSENS.Peak.Acceleration(inc,:) = max(Acceleration_Filtered(iData-length(O)/2:iData+length(O)/2,:));
            Features_XSENS.Mean.Acceleration(inc,:) = mean(Acceleration_Filtered(iData-length(O)/2:iData+length(O)/2,:));
            Features_XSENS.Median.Acceleration(inc,:) = median(Acceleration_Filtered(iData-length(O)/2:iData+length(O)/2,:));
            Features_XSENS.STD.Acceleration(inc,:) = std(Acceleration_Filtered(iData-length(O)/2:iData+length(O)/2,:));
            Features_XSENS.Percentile10.Acceleration(inc,:) = prctile(Acceleration_Filtered(iData-length(O)/2:iData+length(O)/2,:),10);
            Features_XSENS.Percentile25.Acceleration(inc,:) = prctile(Acceleration_Filtered(iData-length(O)/2:iData+length(O)/2,:),25);
            Features_XSENS.Percentile75.Acceleration(inc,:) = prctile(Acceleration_Filtered(iData-length(O)/2:iData+length(O)/2,:),75);
            Features_XSENS.Percentile90.Acceleration(inc,:) = prctile(Acceleration_Filtered(iData-length(O)/2:iData+length(O)/2,:),90);
            inc = inc + 1;
        end
        
%         % Jerk
%         disp([num2str(iSubjects) '.Jerk'])
%         O = zeros(porte,size(Jerk,2));
%         inc = 1;
%         for iData = round(cycles_Valentin(1)*60):length(O):round(cycles_Valentin(end)*60)
%             Features_XSENS.Peak.Jerk(inc,:) = max(Jerk(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.Mean.Jerk(inc,:) = mean(Jerk(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.Median.Jerk(inc,:) = median(Jerk(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.STD.Jerk(inc,:) = std(Jerk(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.Percentile10.Jerk(inc,:) = prctile(Jerk(iData-length(O)/2:iData+length(O)/2,:),10);
%             Features_XSENS.Percentile25.Jerk(inc,:) = prctile(Jerk(iData-length(O)/2:iData+length(O)/2,:),25);
%             Features_XSENS.Percentile75.Jerk(inc,:) = prctile(Jerk(iData-length(O)/2:iData+length(O)/2,:),75);
%             Features_XSENS.Percentile90.Jerk(inc,:) = prctile(Jerk(iData-length(O)/2:iData+length(O)/2,:),90);
%             inc = inc + 1;
%         end
%         
%         % Module Angular Velocity
%         disp([num2str(iSubjects) '.Module_Angular_Velocity'])
%         O = zeros(porte,size(Module_Angular_Velocity,2));
%         inc = 1;
%         for iData = round(cycles_Valentin(1)*60):length(O):round(cycles_Valentin(end)*60)
%             Features_XSENS.Peak.Module_Angular_Velocity(inc,:) = max(Module_Angular_Velocity(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.Mean.Module_Angular_Velocity(inc,:) = mean(Module_Angular_Velocity(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.Median.Module_Angular_Velocity(inc,:) = median(Module_Angular_Velocity(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.STD.Module_Angular_Velocity(inc,:) = std(Module_Angular_Velocity(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.Percentile10.Module_Angular_Velocity(inc,:) = prctile(Module_Angular_Velocity(iData-length(O)/2:iData+length(O)/2,:),10);
%             Features_XSENS.Percentile25.Module_Angular_Velocity(inc,:) = prctile(Module_Angular_Velocity(iData-length(O)/2:iData+length(O)/2,:),25);
%             Features_XSENS.Percentile75.Module_Angular_Velocity(inc,:) = prctile(Module_Angular_Velocity(iData-length(O)/2:iData+length(O)/2,:),75);
%             Features_XSENS.Percentile90.Module_Angular_Velocity(inc,:) = prctile(Module_Angular_Velocity(iData-length(O)/2:iData+length(O)/2,:),90);
%             inc = inc + 1;
%         end
%         
%         % Module Acceleration
%         disp([num2str(iSubjects) '.Module_Acceleration'])
%         O = zeros(porte,size(Module_Acceleration,2));
%         inc = 1;
%         for iData = round(cycles_Valentin(1)*60):length(O):round(cycles_Valentin(end)*60)
%             Features_XSENS.Peak.Module_Acceleration(inc,:) = max(Module_Acceleration(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.Mean.Module_Acceleration(inc,:) = mean(Module_Acceleration(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.Median.Module_Acceleration(inc,:) = median(Module_Acceleration(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.STD.Module_Acceleration(inc,:) = std(Module_Acceleration(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.Percentile10.Module_Acceleration(inc,:) = prctile(Module_Acceleration(iData-length(O)/2:iData+length(O)/2,:),10);
%             Features_XSENS.Percentile25.Module_Acceleration(inc,:) = prctile(Module_Acceleration(iData-length(O)/2:iData+length(O)/2,:),25);
%             Features_XSENS.Percentile75.Module_Acceleration(inc,:) = prctile(Module_Acceleration(iData-length(O)/2:iData+length(O)/2,:),75);
%             Features_XSENS.Percentile90.Module_Acceleration(inc,:) = prctile(Module_Acceleration(iData-length(O)/2:iData+length(O)/2,:),90);
%             inc = inc + 1;
%         end
%         
%         % Module Jerk
%         disp([num2str(iSubjects) '.Module_Jerk'])
%         O = zeros(porte,size(Module_Jerk,2));
%         inc = 1;
%         for iData = round(cycles_Valentin(1)*60):length(O):round(cycles_Valentin(end)*60)
%             Features_XSENS.Peak.Module_Jerk(inc,:) = max(Module_Jerk(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.Mean.Module_Jerk(inc,:) = mean(Module_Jerk(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.Median.Module_Jerk(inc,:) = median(Module_Jerk(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.STD.Module_Jerk(inc,:) = std(Module_Jerk(iData-length(O)/2:iData+length(O)/2,:));
%             Features_XSENS.Percentile10.Module_Jerk(inc,:) = prctile(Module_Jerk(iData-length(O)/2:iData+length(O)/2,:),10);
%             Features_XSENS.Percentile25.Module_Jerk(inc,:) = prctile(Module_Jerk(iData-length(O)/2:iData+length(O)/2,:),25);
%             Features_XSENS.Percentile75.Module_Jerk(inc,:) = prctile(Module_Jerk(iData-length(O)/2:iData+length(O)/2,:),75);
%             Features_XSENS.Percentile90.Module_Jerk(inc,:) = prctile(Module_Jerk(iData-length(O)/2:iData+length(O)/2,:),90);
%             inc = inc + 1;
%         end
        
        % Data saving
        cd(['V:\Piano_Fatigue\Data_Exported\Features_XSENS'])
        save(['Features3_XSENS_' (Filenames{2,iSubjects}) '.mat'],'Features_XSENS')

%         cd(['J:\Piano_Fatigue\Data_Exported\Features_XSENS_Cycles'])
%         save(['Features3_XSENS_Cycles_' (Filenames{2,iSubjects}) '.mat'],'Features_XSENS_Cycles')
        
        PartInc = PartInc + 1;
        clear functions
end
    