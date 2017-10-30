% wrapper function for preprocessing of PDT fMRI data only

% TEST RUN: copy one subject into a different directory!!!!

function [error_msg] = PDT_preprocess(cur_sub,which_folders,pm_defaults_file,preproc_tmpl,ow,run_pp,cdim_resl)

%% what to preprocess?

cd(cur_sub)
cur_home = [cur_sub,'\MRT\NIFTI\PDT'];         % added: current home is PDT folder
try
    cd(cur_home)
catch cur_error
    error_msg = cur_error;
    disp(error_msg)
    return
end
disp(['Now running ' cur_sub]);
folders_to_pp = {};
which_folders= {'*epi*', '*MoCo*'}


%% UPDATE MATLAB 2017 -> function 'contains' would be very useful here.... ?
% Does folder contain epi or moco files?

for ii = 1:length(which_folders)
    tmp = cellstr(ls(which_folders{ii}));
    if strcmp(tmp,'') == 1
        if ii == 1
            epi = 0;
            moco = 1;
        elseif ii == 2
            moco = 0;
            epi = 1;
        end
    end
    folders_to_pp = [folders_to_pp;cellstr(tmp)];
end
what_folders = [];
for ii = 1:length(folders_to_pp)
    if isdir(folders_to_pp{ii})
        what_folders = [what_folders,ii];
    end
end
folders_to_pp = folders_to_pp(what_folders);


%% what fieldmaps do we have?
% TODO here: Go into each folder (PDT, REST, SLM) and find fieldmaps there
cur_which_folder = {'*Fieldmap'};
folders_fm = {};
for ii = 1:length(cur_which_folder)
    tmp = cellstr(ls(cur_which_folder{ii}));
    folders_fm = [folders_fm;cellstr(tmp)];
end
what_folders = [];
for ii = 1:length(folders_fm)
    if isdir(folders_fm{ii})
        what_folders = [what_folders,ii];
    end
end
folders_fm = folders_fm(what_folders);

% Check dimensions and reslice if necessary
if cdim_resl
    all_fieldmaps = {};     % help variables: for fieldmaps and reference
    reference = {};
    P = {};     % list preprocessing folders
    cd(cur_home)
    for ii = 1:length(folders_to_pp)
        cd(folders_to_pp{ii})
        P{ii}= spm_select('FPList',pwd,'.nii');         % list contents of folders to pp
        VY{ii} = spm_vol(P{ii});
        cd(cur_home)
    end
    % check dimensions
    x_true = []; % help variables
    y_true = [];
    z_true = [];

for v = 1:size(VY,2)
    for ii = 1:length(VY{1,1}(:,1)) % logical matrices to check dimensions
    x_true{v}(ii) = (VY{1,v}(ii).dim(1) == 64);
    y_true{v}(ii) = (VY{1,v}(ii).dim(2) == 64);
    z_true{v}(ii) = (VY{1,v}(ii).dim(3) == 33);
    end
end

for v = 1:size(VY,2)        % if one dimension doesn't fit -> error
    if x_true{v} == 0 | y_true{v} == 0 | z_true{v} == 0 
    error('Please revise dimensions.')
    end
end
    F = {};         % fieldmap folder
    for jj = 1:length(folders_fm)
        cd(cur_home)
        cd(folders_fm{jj})
        F{jj} = spm_select('FPList',pwd,'.nii'); 
        FVY{jj} = spm_vol(F{jj});
    end
    
    % check if fieldmaps and epis are the same
    if (isequal(FVY{1,1}(1).mat, VY{1,1}(1).mat) == 0)     % if fieldmaps and epi are not equal
        reference = VY{1,1}(1).fname;           % take first epi as reference
        all_fieldmaps{1} = FVY{1}(1).fname;     % create fieldmap structure array
        all_fieldmaps{2} = FVY{2}(2).fname;
        all_fieldmaps{3} = FVY{1,2}(1).fname;
        all_fieldmaps{4} = FVY{1,2}(2).fname;
        agk_nii_in_new_space(reference,all_fieldmaps,all_fieldmaps)  % reslice
    end
end

%% what are the PDT folders and what are the SLM folders?
% checking if there were two runs or not in this subject, because
% depending on this we'll have to use a one session or two-session
% preprocessing batch template

% are there two or one runs in this subject?

cd(cur_sub)         % changed: cur_sub instead of cur_home -> cur_home is PDT folder
cd('Behav\PDT')     % changed -> before ..\..\Behav\PDT
%change: list P_* data and load the first one!
% BEHAV = cellstr(ls('P_*'));
% load(BEHAV{1});
load(ls('P_*'))
tworuns   = isfield(P.t,'triggerscannertpostpause');
%% batch preprocessing for all to-be-preprocessed folders
% now different: one for loop through each folder
% if tworuns then expecting both to be either epi or moco files in the PDT
% (!!!) experiment

% changed: go again into PDT folder (=home) we are in the subjects folder now


cd(cur_sub)
cd(cur_home)

% changed: two runs, if statement not in loop
if tworuns 
  load(preproc_tmpl{2})     % load the two session preproc template if there are two runs
else
    load(preproc_tmpl{1})
end
for jj = 1:length(folders_to_pp)        % loop through the folders to pp
        cd(cur_home) % first fo home again (every loop)
        cd(folders_to_pp{jj});  % then, go into subject folder

        % delete old preproc batch if ow == 1
        if ow == 1
            if exist('Preprocessing_pipeline_SPM12.mat')
                delete('Preprocessing_pipeline_SPM12.mat')
            end
        end
        
           if epi
            % here getting the fnames; try using faster function?
            % checking if already preprocessed, only if ow is not allowed
            
            if ow == 0
                list = cellstr(ls('swuaepi*'));
                if length(list)> 50
                    error_msg{jj} = ['The subject ' cur_sub ' series ' folders_to_pp{jj} '... has already been preprocessed!'];
                    disp(error_msg{jj})
                    disp('I will continue to next series.')
                    continue
                end
            end
            % checking if it is a too short measurement
            list = cellstr(ls('epi_*'));
            if length(list) < 250
                error_msg{jj} = ['The subject ' cur_sub ' series ' folders_to_pp{jj} '... is too short to be real. I ignore it!'];
                disp(error_msg{jj})
                disp('I will continue to next series.')
                continue
            end
                    %% Do the same for the moco files
           elseif moco          % if we have moco files
           if ow == 0
            list = cellstr(ls('swuaepi*'));
            if length(list)> 50
                error_msg{jj} = ['The subject ' cur_sub ' series ' folders_to_pp{jj} '... has already been preprocessed!'];
                disp(error_msg{jj})
                disp('I will continue to next series.')
                continue
            end
           end

        % checking if it is a too short measurement
        list = cellstr(ls('*MoCo*'));               % list moco files
        if length(list) < 250
            error_msg{jj} = ['The subject ' cur_sub ' series ' folders_to_pp{jj} '... is too short to be real. I ignore it!'];
            disp(error_msg{jj})
            disp('I will continue to next series.')
            continue
        end
     end
           
            %% MORE ELEGANT HERE: index in matlabbatch and overwrite fnames or concatenate second fnames to first fnames???
          
            % get scans
            if (jj == 1)
                cd(cur_home)
                cd(folders_to_pp{1});
                disp(['Specify slice timing of Session 1... ' jj])
                fnames           = spm_select('FPList',pwd,'^epi_');
                fnames_rel_1     = cellstr(spm_select('List',pwd,'^epi_'));
                fnames_rel_1_pwd = pwd;
                % fill batch
                matlabbatch{1,1}.spm.temporal.st.scans{1,1} = cellstr(fnames);
                first_epi = fnames(1,:);
            elseif (jj == 2)
                disp('Specify slice timing (Session 2)');
                cd(cur_home)
                cd(folders_to_pp{2});
                % fill batch
                fnames_2 = spm_select('FPList',pwd,'^epi');
                matlabbatch{1,1}.spm.temporal.st.scans{1,2} = cellstr(fnames_2);
                fnames_rel_2 = cellstr(spm_select('List',pwd,'^epi'));
                fnames_rel_2_pwd = pwd;
                scnd_epi = fnames_2(1,:);
                
            end

        
        %% CHANGE: always both fieldmaps in folder used, no matter if it is one or two runs
        % taken out of if statement that asks for one or two runs
       % select which fieldmap to use
       % folders to be picked
       switch numel(folders_fm)
           case 4 
               error('Specify the fieldmaps you want to use.')
           case 3
               error('Specify the fieldmaps you want to use.')
           otherwise
               fm_to_be_picked = 1;
        end
    
    %% now tworuns or onerun PDT is all decided; fill the rest of batch
    
    % fill batch with phase image
    disp('Specify field map (phase image)');
    cd(cur_home)
      try
        % change: cd to pwd -> did not work with cd (index exceeds matrix
        % dimensions)
        cd(folders_fm{fm_to_be_picked+1})       % fieldmap with phase image is the second fieldmap folder
    catch
       cd(folders_fm{2})   % in case only one fieldmap was recorded - should not be the case
    end
    matlabbatch{1,2}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = cellstr(spm_select('FPList',pwd,'^Field'));
    
    % fill batch with magnitude images
    disp('Specify field map (1 of 2 magnitude images)');
    cd(cur_home)
    try
        % again cd changed to pwd
        cd(folders_fm{fm_to_be_picked})
    catch
        cd(folders_fm{1}) % in case only one fieldmap was recorded
    end
    m_images = spm_select('FPList',pwd,'^Fieldmap');
    % fill batch
    matlabbatch{1,2}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = cellstr(m_images(1,:)); % will be the first one
    % fill batch
    matlabbatch{1,2}.spm.tools.fieldmap.calculatevdm.subj.session.epi = cellstr(first_epi);
    % specifiying the defaults file
    matlabbatch{1,2}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsfile{1,1} = pm_defaults_file;
    
    %% Find MPRAGE folder
    % get the coregistration source image
    cd(cur_sub)
    cd('MRT\NIFTI\MPRAGE')         % go back to the subjects folder
    try
        cd(ls('*MPRAGE'))       % go into MPRAGE folder and find MPRAGE subfolder if there is no subfolder, produce error
    catch cur_error
        error_msg = cur_error;
        disp(error_msg)
        return
    end
    
    disp('Specify coregistration (source image)');
    % filter here is start with "MPRAGE" and end with ".nii"
    % fill batch
    matlabbatch{1,4}.spm.spatial.coreg.estimate.source = {spm_select('FPList',pwd,'^MPRAGE.*\.nii')};
    disp('Specify coregistered t1 image for SPM segmentation');
    % fill batch
    matlabbatch{1,5}.spm.spatial.preproc.channel.vols = matlabbatch{1,4}.spm.spatial.coreg.estimate.source;
    
    % specify images for normalise:write
    disp('Specify images for normalise:write');
    cd(cur_home)
    cd(folders_to_pp{jj});
    %% CHANGE: if statement -> fnames_rel_1 must not change with every iteration!!
    if (jj==1)
    fnames_rel_1 = agk_add_prefix(fnames_rel_1,[fnames_rel_1_pwd filesep matlabbatch{1,3}.spm.spatial.realignunwarp.uwroptions.prefix matlabbatch{1,1}.spm.temporal.st.prefix]);
    end
    if  exist('fnames_rel_2')
        fnames_rel_2 = agk_add_prefix(fnames_rel_2,[fnames_rel_2_pwd filesep matlabbatch{1,3}.spm.spatial.realignunwarp.uwroptions.prefix matlabbatch{1,1}.spm.temporal.st.prefix]);
    end
    if exist('fnames_rel_2')
        matlabbatch{1,6}.spm.spatial.normalise.write.subj.resample = [fnames_rel_1; fnames_rel_2];
    else
        matlabbatch{1,6}.spm.spatial.normalise.write.subj.resample = fnames_rel_1;
    end
    
    % t1 normalise_write
    matlabbatch{1,8}.spm.spatial.normalise.write.subj.resample = matlabbatch{1,4}.spm.spatial.coreg.estimate.source;
    
    save('Preprocessing_pipeline_SPM12.mat','matlabbatch');
    
    % after second loop iteration: go back into first session's folder and
    % store the final batch there
    if jj == 2
        cd(cur_home)
        cd(folders_to_pp{1})
        save('Preprocessing_pipeline_SPM12.mat','matlabbatch');
    end
end
    
%% Add case 'subject has already been preprocessed' - otherwise there will be an error
    if run_pp
        try
            spm_jobman('initcfg');
            spm_jobman('run',matlabbatch);
            error_msg = ['The subject ' cur_sub ' series ' folders_to_pp{jj} '... was successful'];
            disp(error_msg)
        catch cur_error
            disp(['The subject ' cur_sub ' series ' folders_to_pp{jj} '.. produced an ERROR']);
            error_msg = cur_error;
            disp(error_msg)
        end
           
         
    else
        error_msg = ['The subject ' cur_sub ' series ' folders_to_pp{jj} '... did not run by user''s wish.'];
    end

end    




