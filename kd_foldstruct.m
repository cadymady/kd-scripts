%%% CHECK FOLDER STRUCTURE %%%
% Checks if folder contains subfolders: MPRAGE, PDT, SLM, REST,OTHER
% and if PDT, SLM, and REST contain enough EPI Images
% and if these folders contain the correct (amount) of files
%for NIFTI converted fMRI DATA - 

%   scrpt created for VPPG Data at Charité Berlin, August 2017

% Input args: cur_dir: current directory in which folder of subjects are
% stored
%             PDT: approximate number of MoCo files that should be in the
%             PDT, same for SLM, REST
% Output: subj: participant's ID

% Author: Katharina Dücker
% Date: 15.08.2017

%% ADD OPTIONS: PDT, REST, SLM ...

clear all
close all
clc
cur_dir = 'E:\data\';
PDT = 1500;
SLM = 700;

cd(cur_dir);    % current directory
all_subf = cellstr(ls('VPPG*'));  %list folders starting w vppg
check_folders = {'MPRAGE','OTHER','PDT','REST','SLM'};           % cell array that contains subfolder names

    for ii = 36:37 % for all subfolders
        cd([cur_dir,all_subf{ii},'\MRT\NIFTI\'])  % go into NIFTI folder
        folder_cons = cellstr(ls([cur_dir,all_subf{ii},'\MRT\NIFTI\']));        % list elements of current folder
        folder_cons = folder_cons(strncmp('.',folder_cons, 1) == 0);            % delete random point elements
        for nn = 1:length(check_folders)
            if numel(folder_cons) > 5
                subj = ii;
                error(['Subject ',all_subf{ii}, ' too many subfolders.'])
                
            elseif ismember(check_folders{nn}, folder_cons) == 0        % if subfolder is not in the folder list
                subj = ii;
                error(['Subject ',all_subf{ii}, ' is at least missing one subfolder. Please revise the folder structure.'])
           
                         
            elseif ismember(check_folders{nn}, folder_cons) == 1
                    cd([cur_dir,all_subf{ii},'\MRT\NIFTI\',check_folders{nn}])           % go into that folder
                
                %Check for MPRAGE if MPRAGE folder exists
                    if strcmp(check_folders{nn}, 'MPRAGE') == 1
                   MPRAGE = cellstr(ls('*MPRAGE*'));        % check if there is at least one folder that is called MPRAGE
                   switch numel(MPRAGE)
                       case 0
                       subj = ii;
                       error(['MPRAGE missing in Subject ',all_subf{ii}]) % otherwise: error message
                   end
               
                % check for PDT if there is enough files
                elseif strcmp(check_folders{nn}, 'PDT') == 1        % 
                    MOCO = cellstr(ls('*MoCoSeries*'));
                    EPI = cellstr(ls('*epi*'));
                    switch numel(MOCO) 
                        % if there is two MoCo folders
                        case 2
                            mtmp1 = dir([cd,'\',MOCO{1}]);           % help variable: all filenames in first moco
                            mtmp2 = dir([cd,'\',MOCO{2}]);           % all filenames in second
                            if length(mtmp1)+length(mtmp2) < PDT        % if there sum is smaller than the criterion (set in function)
                                subj = ii;
                                error(['Not enough MoCo files in PDT. Subject ', all_subf{ii}])
                            end
                        %if there is only one MoCo folder
                        case 1
                            mtmp = dir([cd,'\',MOCO{1}]);
                            if strcmp(MOCO{1},'') == 1
                                MOCO = [];
                                                           
                            elseif length(mtmp) < PDT   % if the folder contains less elements than set by criterion, produce an error
                                subj = ii;
                                error(['Not enough MoCo files in PDT. Subject ', all_subf{ii}])
                            end

                            % if there is no MoCos, check the above for epis
                        case 0
                            switch numel(EPI)
                                case 2
                                etmp1 = dir([cd,'\',EPI{1}]);           
                                etmp2 = dir([cd,'\',EPI{2}]);         
                            if length(etmp1)+length(etmp2) < PDT        % if there sum is smaller than the criterion (set in function)
                               subj = ii;
                                error(['Not enough epi files in PDT. Subject ', all_subf{ii}])
                            end
                        %if there is only one EPI folder
                            case 1
                            etmp = dir([cd,'\',EPI]);
                            if length(etmp) < PDT   % if the folder contains less elements than set by criterion, produce an error
                                subj = ii;
                                error(['Not enough EPI files in PDT. Subject ', all_subf{ii}])
                            end

                            % crash if there is neither MoCo nor EPI files
                            %% PROBABLE ERROR MESSAGE: OTHERWISE IF NONE OF THE IF STATEMENTS ABOVE HOLDS???
                                otherwise
                                    subj = ii;
                                    error(['No MoCo or EPI files in PDT folder, subject ', all_subf{ii}])
                            end
                    end


                    % Check if there are fieldmaps in REST

                     elseif strcmp(check_folders{nn}, 'REST') == 1
                         FIELD = cellstr(ls('*Fieldmap*'));
                         MOCO = cellstr(ls('*MoCoSeries*'));
                         EPI = cellstr(ls('*epi*'));
                         if numel(FIELD) < 2
                             subj = ii;
                             error(['Not enough fieldmaps in REST Subject ', all_subf{ii}])
                         end

                         switch numel(MOCO)
                             case 0
                                 switch numel(EPI)
                                     case 0
                                         subj = ii;
                                         error(['MoCo and EPI files in folder REST, Subject ', all_subf{ii}, ' missing.'])
                                 end
                         end

                         % check for number of elements in SLM
                     elseif strcmp(check_folders{nn}, 'SLM') == 1
                           MOCO = cellstr(ls('*MoCoSeries*'));
                            EPI = cellstr(ls('*epi*'));
                            switch numel(MOCO) 
                        % if there is two MoCo folders
                        case 2
                            mtmp1 = dir([cd,'\',MOCO{1}]);           % help variable: all filenames in first moco
                            mtmp2 = dir([cd,'\',MOCO{2}]);           % all filenames in second
                            if length(mtmp1)+length(mtmp2) < SLM        % if there sum is smaller than the criterion (set in function)
                               subj = ii;
                                error(['Not enough MoCo files in SLM. Subject ', all_subf{ii}])
                            end
                        %if there is only one MoCo folder
                        case 1
                            mtmp = dir([cd,'\',MOCO{1}]);
                             if strcmp(MOCO{1},'') == 1
                                MOCO = [];
                             elseif length(mtmp) < PDT   % if the folder contains less elements than set by criterion, produce an error
                                subj = ii;
                                error(['Not enough MoCo files in SLM. Subject ', all_subf{ii}])
                            end

                            % if there is no MoCos, check the above for epis
                        otherwise
                            switch numel(EPI)
                                case 2
                                etmp1 = dir([cd,'\',EPI{1}]);           
                                etmp2 = dir([cd,'\',EPI{2}]);         
                            if length(etmp1)+length(etmp2) < SLM        % if there sum is smaller than the criterion (set in function)
                               subj = ii;
                                error(['Not enough epi files in SLM. Subject ', all_subf{ii}])
                            end
                        %if there is only one EPI folder
                        case 1
                            etmp = dir([cd,'\',{EPI}]);
                            if length(etmp) < PDT   % if the folder contains less elements than set by criterion, produce an error
                                subj = ii;
                                error(['Not enough EPI files in SLM. Subject ', all_subf{ii}])
                            end

                            % crash if there is neither MoCo nor EPI files
                            %% PROBABLE ERROR MESSAGE: OTHERWISE IF NONE OF THE IF STATEMENTS ABOVE HOLDS???
                                otherwise
                                    subj = ii;
                                    error(['No MoCo or EPI files in SLM folder, subject ', all_subf{ii}])
                            end
                            end


                end
            end
        end
    end



