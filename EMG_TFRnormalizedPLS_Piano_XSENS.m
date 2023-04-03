clear ;
close all ;
clc ;

addpath C:\Users\p1098713\Documents\3.MATLAB\Fonctions
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\btk'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\wtc-r16'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Matlab-Functions-Git'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\FFT'))
addpath(genpath('C:\Users\p1098713\Documents\3.MATLAB\Fonctions\Entropy Fouaz'))

FoldersNames = dir('J:\Piano_Fatigue\Data_Compiled');
Filenames = {...
    FoldersNames(4:end).name...
   } ;

    %Crée une deuxième ligne avec le nom du fichier sans .mat
    Filenames{2,1}=[];
    for iFilenames = 1:length(Filenames)
        Filenames{2,iFilenames} = erase(Filenames{1,iFilenames},'.mat');
    end
    
for iSubjects = 2:2:100
    % Load XSENS data
    cd(['J:\Piano_Fatigue\Data_Compiled'])
    load([Filenames{1,iSubjects}]);
    Angular_Velocity = table2array(data.angularVelocity(:,[4:6 13:15 19:21 22:33])) ; % [13:15 19:21] pour L5, T8, head... & 22:33 pour Rshoulder à Rhand
    Acceleration = table2array(data.acceleration(:,[4:6 13:15 19:21 22:33])) ;
    Quaternion = table2array(data.orientation(:,[5:8 17:20 25:44])) ;
    Freq = 60 ;
    dt=1/Freq ;
    t=(0:length(Angular_Velocity)-1)*dt ;
    
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
%             line([cycles(:,1),cycles(:,1)]', repmat(ylim,length(cycles(:,1)),1)','color','red')
%             line([cycles(:,end),cycles(:,end)]', repmat(ylim,length(cycles(:,end)),1)','color','red')
%             spec_fft(Acceleration_Filtered(:,iM),Freq,1);
%             pause()
%         end

    % Module
    Module_Angular_Velocity = [];
    Module_Acceleration = [];
    VarInc = 1;
    for icol = 1:3:size(Angular_Velocity,2)
        Module_Angular_Velocity(:,VarInc) = Module(Angular_Velocity_Filtered(:,icol:icol+2));
        Module_Acceleration(:,VarInc) = Module(Acceleration_Filtered(:,icol:icol+2));
%         subplot(4,1,1) ; plot(t,Module_Acceleration(:,VarInc))
%         subplot(4,1,2) ; plot(t,Acceleration_Filtered(:,icol))
%         subplot(4,1,3) ; plot(Acceleration_Filtered(:,icol+1))
%         subplot(4,1,4) ; plot(Acceleration_Filtered(:,icol+2))
%         pause()
        VarInc = VarInc+1 ;
    end
    
%     Normalized_Acceleration_Filtered = interp1(1:length(Acceleration_Filtered),Acceleration_Filtered,linspace(1,length(Acceleration_Filtered),10000)) ;
%     for iC = 1:size(Acceleration_Filtered,2)
%         subplot(2,1,1); plot(Acceleration_Filtered(:,iC))
%         subplot(2,1,2); plot(Normalized_Acceleration_Filtered(:,iC))
%         pause()
%     end
    
%% EMG Median frequency
        f = 0.05:0.05:24.05;...Range de fréqence de la décomposition
        ond = ["cmor8-1"];
        scale = Freq*centfrq(ond)./f;
        TFR = [];
        % Angular Velocity
        disp([num2str(iSubjects) '.Angular Velocity'])
        for iData = 1:size(Angular_Velocity,2)
            V = ['V' num2str(iData)];
            [norm, Period] = cwt(Angular_Velocity_Filtered(:,iData),scale,ond,dt);
            Norm_norm = interp1(1:length(abs(norm)),abs(norm)',linspace(1,length(abs(norm)),1000))' ;
            Norm_t = interp1(1:length(t),t',linspace(1,length(t),1000))' ;
            TFR.AngularVelocity.(V) = Norm_norm;
%             figure; pcolor(t,f,abs(norm)); shading interp
%             figure; pcolor(Norm_t,f,Norm_norm); shading interp
        end
        
        % Acceleration
        disp([num2str(iSubjects) '.Acceleration'])
        for iData = 1:size(Acceleration,2)
            V = ['V' num2str(iData)];
            [norm, Period] = cwt(Acceleration_Filtered(:,iData),scale,ond,dt);
            Norm_norm = interp1(1:length(abs(norm)),abs(norm)',linspace(1,length(abs(norm)),1000))' ;
            Norm_t = interp1(1:length(t),t',linspace(1,length(t),1000))' ;
            TFR.Accelearation.(V) = Norm_norm;
        end
        
        % Module Angular Velocity
        disp([num2str(iSubjects) '.Module Angular Velocity'])
        for iData = 1:size(Module_Angular_Velocity,2)
            V = ['V' num2str(iData)];
            [norm, Period] = cwt(Module_Angular_Velocity(:,iData),scale,ond,dt);
            Norm_norm = interp1(1:length(abs(norm)),abs(norm)',linspace(1,length(abs(norm)),1000))' ;
            Norm_t = interp1(1:length(t),t',linspace(1,length(t),1000))' ;
            TFR.Module_AngularVelocity.(V) = Norm_norm;
        end
        
        % Module Acceleration
        disp([num2str(iSubjects) '.Module Acceleration'])
        for iData = 1:size(Module_Acceleration,2)
            V = ['V' num2str(iData)];
            [norm, Period] = cwt(Module_Acceleration(:,iData),scale,ond,dt);
            Norm_norm = interp1(1:length(abs(norm)),abs(norm)',linspace(1,length(abs(norm)),1000))' ;
            Norm_t = interp1(1:length(t),t',linspace(1,length(t),1000))' ;
            TFR.Module_Acceleration.(V) = Norm_norm;
        end
        TFR.Time = Norm_t;
        
        % Data saving
        cd(['J:\Piano_Fatigue\Data_Exported\TFRnormalized_XSENS'])
        save(['TFR_XSENS_' (Filenames{2,iSubjects}) '.mat'],'Features_XSENS')

        clear functions
end
    