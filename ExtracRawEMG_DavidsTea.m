clear ;
close all ;
clc ;

addpath E:\Bureau\Etienne\MATLAB\Functions
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\btk'))
addpath(genpath('E:\Bureau\Etienne\MATLAB\Functions\wtc-r16'))

Muscles = {...
    'DeltA',  'DeltA_IM_EMG1';...
    'DeltM',  'DeltM_IM_EMG2';...
    'DeltP',    'DeltP_IM_EMG3';...
    'Bi',    'Bi_IM_EMG4';...
    'Tri',  'Tri_IM_EMG5';...
    'TrapSup',  'TrapSup_IM_EMG6';...
    'TrapMed'  'TrapMed_IM_EMG7';...
    'TrapInf'  'TrapInf_IM_EMG8';...
    'Dent'  'Dent_IM_EMG9';...
    };

        cd(['J:\IRSST_DavidsTea\2022_07_28_Testing\Pilote2\EMG'])
        acq = btkReadAcquisition('Contractions_Alexis01.c3d') ;
        Data = btkGetAnalogs(acq) ;
        Freq = btkGetAnalogFrequency(acq) ;
        
        DataSelec = {};
        for iM = 1:length(Muscles)
            DataSelec.(Muscles{iM,1}) = Data.(Muscles{iM,2});
        end
        
        for iM = 1:length(Muscles)
            subplot(5,2,iM); 
            plot(DataSelec.(Muscles{iM}));
%             hold on
            title(Muscles{iM});
            axis([0 4e5 -1.7e-3 1.7e-3])
%             pause()
        end