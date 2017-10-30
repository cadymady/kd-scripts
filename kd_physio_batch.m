%-----------------------------------------------------------------------
% VPPG Studie, Charité Berlin 2017, last edit: 18.09.2017, 16:30
% Author: templates: physio_job_alex, KD, FC
%
%-----------------------------------------------------------------------

% additional comments
% creates a matlabbatch job to be run in SPM which will generate a
% design matrix with physio regressors
% only for first session of the data -> preparation for

% where u wanna get all results and output files saved?
function [error_msg, error_sub] = kd_physio_batch( base_path, cur_sub, session, ow )

% kd_physio_batch fills in batch and runs physio correction of each subject
% check if a results folder already exist, overwrite or leave

sub_path = fullfile(base_path, cur_sub); % full subject path
cd(sub_path)

disp(['Starting physio correction for subject ', cur_sub])

cd('Physio')
% fill first part of the bitch: where to create the results folder?
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = {pwd};
save_dir = fullfile(pwd,'results');  % where do you want to save the results/the batch?

%% Overwrite existing Physio correction?
% create saving directory
if exist(save_dir, 'file')
    if ow
        cmd_rmdir(save_dir)
        disp(['Existing Physio results folder removed.'])
        % create results directory
        matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'results';
        mkdir(save_dir)
    else
        error_msg = ['Physio results already exist. I will continue to the next subject.'];
        return
    end
else
    % create results directory
    matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'results';
    mkdir(save_dir)
end

% dependency: where is the save directory?
matlabbatch{2}.spm.tools.physio.save_dir(1) = cfg_dep('Make Directory: Make Directory ''results''', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));


%% define all folders necessary (dicoms, niftis etc)
% works for both epi and MoCo files
cd([sub_path, '\MRT\'])
dcm = ls('1*');             % dcm folder

cd(['NIFTI\', session])    % go into nifti folder (niftis are organized by session)

% list MoCo or epi folders
which_folders = {'*MoCo*', '*epi*'};
for jj = 1:length(which_folders)
    tmp = cellstr(ls(which_folders{jj}));
    for t = 1:length(tmp)
        if exist(tmp{t}, 'file')       % if tmp is a file in the folder
            nif_sess{t} = tmp{t};           % store it as 'nif_sess'
        end
    end
end
cd(sub_path)         % go back home

 cd('Physio')
    % list all physio files
    physio_files = cellstr(ls(['*',session,'*']));
    
% loop over all sessions: one batch for each session
for nn = 1:length(nif_sess) 
    %% fill in batch
    % what scanner did you use?
    matlabbatch{2}.spm.tools.physio.log_files.vendor = 'Siemens';
    
    % load in cardiac files
    if exist([sub_path, '\Physio\', physio_files{nn}])
        % store path of physio file in cardiac
        matlabbatch{2}.spm.tools.physio.log_files.cardiac = {[sub_path, '\Physio\', physio_files{nn}]};
    else
        error_msg = ['Physio file does not exist. Something wrong with path? \r\n'...
            [sub_path, '\Physio\', cur_sub,'_',session,'_',num2str(nn),'.puls']];
        continue
    end
    
    %% where is the .puls file?
    disp(['For subject ', cur_sub, ' I found ' num2str(numel(nif_sess)), ...
        ' runs. Filling batch for session ', session,' run ', num2str(nn)]);
    
    % store .puls file in batch (if exists)
   
    % we did not not do respiration so leave blank with {''}
    matlabbatch{2}.spm.tools.physio.log_files.respiration = {''};
    
    %% scan timing: use DICOMS
    % indidicate the path to first (or last) .dcm
    % watch out; should order "increasing" (happens automtically in ls() cmd)
    cd(sub_path)
    cd('MRT')
    cd(dcm) % go into dcm folder
    cd(nif_sess{nn})
    % list content of first moco/epi folder(same name as nifti folders!)
    all_dcms = cellstr(ls('1*'));
    try
        % try to use first dicom for the scan timing
        % 'first' if first dicom provided,
        matlabbatch{2}.spm.tools.physio.log_files.scan_timing = ...
            {fullfile(sub_path, '\MRT\', dcm, nif_sess{nn}, all_dcms{1})};
        % document which file you used
        matlabbatch{2}.spm.tools.physio.log_files.align_scan = 'first';
    catch
        % if it doesn't work, use the last dicom
        disp('Last DICOM has been used for the scan timing.')
        matlabbatch{2}.spm.tools.physio.log_files.scan_timing = ...
            {fullfile(sub_path, '\MRT\', dcm, nif_sess{nn}, all_dcms{length(all_dcms)})};
        % 'last' if last dicom provided
        matlabbatch{2}.spm.tools.physio.log_files.align_scan = 'last';
    end
    
  
    %% fill in the rest according to 'physio_job_alex' 
    % ignore for now
    matlabbatch{2}.spm.tools.physio.log_files.sampling_interval = [];
    % ignore for now
    matlabbatch{2}.spm.tools.physio.log_files.relative_start_acquisition = 0;
    
    % you will have to provide this info (check phoenix or...)
    matlabbatch{2}.spm.tools.physio.scan_timing.sqpar.Nslices = 33;
    % ignore
    matlabbatch{2}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
    % TR time (find in phoenix or...)
    matlabbatch{2}.spm.tools.physio.scan_timing.sqpar.TR = 2;
    % do you have dummy scans? I think not
    matlabbatch{2}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
    
    % number of scans = number of dicoms
    matlabbatch{2}.spm.tools.physio.scan_timing.sqpar.Nscans = length(all_dcms);
    % interleaved descending... how did you measure what is first slice?
    % should be your reference slice from preprocessing
    matlabbatch{2}.spm.tools.physio.scan_timing.sqpar.onset_slice = 17;
    % detail: ignore
    matlabbatch{2}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
    % ignore
    matlabbatch{2}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
    % ignore
    matlabbatch{2}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
    % 'PPU' or 'ECG' (I think, check that!)
    matlabbatch{2}.spm.tools.physio.preproc.cardiac.modality = 'PPU';
    % ignore it
    matlabbatch{2}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
    % no idea; an outputfile?
    matlabbatch{2}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
    % no idea
    matlabbatch{2}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
    % output file name
    matlabbatch{2}.spm.tools.physio.model.output_multiple_regressors = ['multiple_regressors_run_',num2str(nn),'.txt'];
    % output file name
    matlabbatch{2}.spm.tools.physio.model.output_physio = ['physio_run_', num2str(nn),'.mat'];
    % ignorieren
    matlabbatch{2}.spm.tools.physio.model.orthogonalise = 'none';
    % retroicore only useful if you have heart rate and breathing
    matlabbatch{2}.spm.tools.physio.model.retroicor.yes.order.c = 3;
    matlabbatch{2}.spm.tools.physio.model.retroicor.yes.order.r = 4;
    matlabbatch{2}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
    matlabbatch{2}.spm.tools.physio.model.rvt.no = struct([]);
    matlabbatch{2}.spm.tools.physio.model.hrv.yes.delays = 0;
    matlabbatch{2}.spm.tools.physio.model.noise_rois.no = struct([]);
  
    % movement params: empty for now
    
    %matlabbatch{2}.spm.tools.physio.model.movement.yes.file_realignment_parameters = {''};
    
    %% using mov params from pp not working
%    if pp done already you can provide the movement params here
    cd(fullfile(sub_path,'MRT\NIFTI', session, nif_sess{nn}))
    rp_file = ls('rp*');
    
    if exist(rp_file, 'file')
        rp_file_path = fullfile(pwd, rp_file);
        matlabbatch{2}.spm.tools.physio.model.movement.yes.file_realignment_parameters = {rp_file_path};
        
        % or leave empty if not yet
    else
        disp(strcat('no movemement parameters available for subject ', cur_sub))
       matlabbatch{2}.spm.tools.physio.model.movement.yes.file_realignment_parameters = {''};
    end
   %%
    % how many movement regressors (default 6, but higher order possible)
    % VOLTERRA movement regressors
    matlabbatch{2}.spm.tools.physio.model.movement.yes.order = 6;
    % maybe automatic censoring also; if too much movement will be excluded
    % (the volume)
    matlabbatch{2}.spm.tools.physio.model.movement.yes.outlier_translation_mm = 1;
    matlabbatch{2}.spm.tools.physio.model.movement.yes.outlier_rotation_deg = 1;
    matlabbatch{2}.spm.tools.physio.model.other.no = struct([]);
    matlabbatch{2}.spm.tools.physio.verbose.level = 2;
    matlabbatch{2}.spm.tools.physio.verbose.fig_output_file = '';
    matlabbatch{2}.spm.tools.physio.verbose.use_tabs = false;

    cd(save_dir)
    
    save(['PhysiO_batch_',session,'_', num2str(nn),'_','.mat'], 'matlabbatch')
    disp(['PhysiO batch filled for subject ', cur_sub, ' session ', session,' ', num2str(nn)]);
    
    spm_jobman('initcfg');  % necessary for dependencies?

    try
        spm_jobman('run',matlabbatch);
        error_msg = ['Physio correction succesful.'];
        error_sub = {''};       
    catch ME
        % document error and which subject produced the error
        error_msg = getReport(ME)
        error_sub = cur_sub;
    end
    
end


end

