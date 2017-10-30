% PIPELINE FOR RUNNING PREPROCESSING AND 1st LEVEL

% Author: Alexander Genauck
% Work address:
% email: alexander.genauck@charite.de
% Website:
% Feb 2016; Last revision: 17-01-2017

%------------- BEGIN CODE --------------
clear all
%% ############### Generell Settings ###############
% User-name
% comp_name = getenv('C:\Program Files\MATLAB');
% Script-Libary-Path: where did you save this script?
base_dir_lib  = 'E:\Master\Praktikum Charité\VPPG\scripts\PDT\ss_models';
base_dir_lib2 = 'T:\Library\MATLAB\PDT\MRI\ss_models\classical';
% where are the subject folders with MRI data?
base_dir_pl = 'E:\Master\Praktikum Charité\VPPG\data';
% Set spm Mask Threshold
sd_level.spm.spmMaskThresh = 0.2; % spm_default: 0.8 on 2nd Lvl
addpath(genpath(base_dir_lib))
addpath(genpath(base_dir_lib2))
%addpath(genpath(base_dir_lib2))
%addpath('T:\Library\MATLAB\PDT\dicm2nii');
addpath('C:\Program Files\spm12');

% which_spm(12,comp_name,1)
% preproc job templates
preproc_tmpl{1}  = fullfile('T:\Library\MATLAB\PDT\pp_subf', 'PDT_preprocessing_1sess_tmpl.mat');
preproc_tmpl{2}  = fullfile('T:\Library\MATLAB\PDT\pp_subf', 'PDT_preprocessing_2sess_tmpl.mat');

% what to run?
d2n         = 0; % dicom 2 nifti
pp          = 0; % preprocessing
ph          = 0; % make physio regressors
ss          = 1; % ss-level
ss_addcons  = 0; % ss-level addcons
ss_newcons  = 0; % ss-level make cons
ss_rm       = 0; % ss-level review design matrices
ss_extr     = 0; % ss-level extract eigvar from ROI

%% ############### DICOM2NIFTI #####################
% for loop for getting into folders
if d2n
    cd(base_dir_pl);
    all_subf = cellstr(ls('VPPG*'));  %list folders starting w vppg
    which_folders = {'*epi*','*MoCoSeries', '*Fieldmap*','*MPRAGE*'}; % patterns to know which folders
    for ii = 1:length(all_subf)
%         cd([base_dir_pl,all_subf{ii}]);
%         cd('MRT');
%         %delete all existing nifti files first
%         system('rmdir /s /q "NIFTI"')

        %go back home
       % cd([base_dir_pl,all_subf{ii}])
        try
            ii
            cv_feedback{ii} = agk_PDT_wrapper_dcm2nifti(all_subf{ii},which_folders);
        catch ME
            disp(ME.identifier);
            disp(['There was an error in d2n. ' all_subf{ii}])
            rethrow(ME);
        end
    end
end

%% ############### PREPROCESS #####################
if pp
    cd(base_dir_pl);
    all_subf         = cellstr(ls('VPPG*'));
    which_folders    = {'*epi*','*MoCoSeries'}; % patterns to know which folders to pp
    pm_defaults_file = fullfile('T:\Library\MATLAB\PDT\pp_subf', 'pm_defaults_AGK.m');
    ow_pp            = 1; % overwrite existing preprocessing?
    run_pp           = 0 % run the PP or just save batch?
    
    % add warning message if ow_pp = 1 -> press key to continue
    
    for ii = 36:size(all_subf,1)
        cd(base_dir_pl);
        try
            pp_feedback{ii} = PDT_preprocess([base_dir_pl,all_subf{ii}],which_folders,pm_defaults_file,preproc_tmpl,ow_pp,run_pp);
        catch MExc
            rethrow(MExc)
            disp(['There was an error in pp. ' all_subf{ii}])
        end
    end
end

%% ############### PHYSIO #########################
if ph
    % only creates physio regressors and saves all the peak detection
    % output; does not run any analysis
    % this must be done in an ss model
    cd(base_dir_pl);
    all_subf         = cellstr(ls('VPPG*'));
    cur_tmpl         = 'C:\Users\genaucka\Google Drive\Library\MATLAB\PDT\ss_models\physio\physio_job.m';
    ow_ph            = 1; % overwrite existing physio regressors?
    run_it           = 1; % run it or just save matlabbatch?
    
        
    for ii = [11,34]
        cd(base_dir_pl);
        try
            ph_feedback{ii} = agk_make_PDT_physio_design_00(cur_sub,cur_tmpl,run_it,ow_ph);
        catch MExc
            rethrow(MExc)
            disp(['There was an error in ph. ' all_subf{ii}])
        end
    end
end


%% ############### ss-LEVEL #######################
if ss
    cd(base_dir_pl);
    all_subf         = cellstr(ls('VPPG*'));
    %preproc_job      = load(['C:\Users\' comp_name '\Google Drive\Library\MATLAB\LA\Preprocessing\SPM12\LA_preprocessing_template_AGK_SPM12.mat']);
    %expl_mask = fullfile('T:\Library\MATLAB\PDT', 'gm_mask.nii'); % grey matter mask
    %expl_mask        = ['C:\Users\' comp_name '\Google Drive\Library\MATLAB\LA\2L_Analysis\ROI\gm_mask.nii'];
    expl_mask = '';
    cur_tmpl{1}      = ['T:\Library\MATLAB\PDT\MRI\ss_models', '\LA_Ss_template_without_acc_reject.mat']; % ss batch templates (1: without acc_reject 2: with accept_reject)
    cur_tmpl{2}      = ['T:\Library\MATLAB\PDT\MRI\ss_models','\LA_Ss_template.mat'];
    cur_tmpl{3}      = ['T:\Library\MATLAB\PDT\MRI\ss_models','\LA_Ss_template_without_acc_reject_2runs.mat'];
   %LA_charpentier   = ['C:\Users\' comp_name '\Google Drive\Promotion\VPPG\VPPG_Exchange\Library\MLE+Model\LA_Charpentier\MLE_results_Charpentier.txt']; % behav model LA Charpentier results
   % LA_classic       = ['C:\Users\' comp_name '\Google Drive\Promotion\VPPG\VPPG_Exchange\Library\MLE+Model\LA_cl\MLE_results_AG.txt'];          % behav model LA classical results
    run_it           = 0;
    aggr             = 3;
    acc_rec          = 1;
    ow               = 1; % overwrite existing ss results?
    fixed_dur        = [];
    do_center        = 1;
    physio_inc       = 1;
    classical        = 1; % estimate classical model?
    
    % THE PDT_SS_DESIGN
    % RERUN THE agk_make_PDT_ss_design_01; few subjects do not have it
    % yet;
    
    for ii = 1%length(all_subf)
        cd(base_dir_pl);
        if strcmp(all_subf{ii}, 'VPPG0115')
            ss_feedback = agk_make_PDT_ss_pspml_pic_0115(all_subf{ii},cur_tmpl,aggr,run_it,acc_rec,expl_mask,ow,do_center,fixed_dur, physio_inc);
        else
            if acc_rec
                if classical
                    ss_feedback = agk_make_PDT_ss_design_edfix_03(all_subf{ii},cur_tmpl,aggr, ...
                        run_it,acc_rec,expl_mask,ow, physio_inc);
                else
                    ss_feedback = agk_make_PDT_ss_pspml_only_pic_tworuns_accint(all_subf{ii},cur_tmpl,aggr, ...
                        run_it,acc_rec,expl_mask,ow,do_center,fixed_dur, physio_inc);
                end
            else
                ss_feedback = agk_make_PDT_ss_pspml_only_pic_tworuns_new(all_subf{ii},cur_tmpl,aggr, ...
                    run_it,acc_rec,expl_mask,ow,do_center,fixed_dur, physio_inc);
            end
        end
        save(['report_sslevel_PDT_ss_design_ed_03_' date '.mat'],'ss_feedback')
    end
        
end

%% ############### Review design matrices #########
if ss_rm
    cd(base_dir_pl);
    all_subf         = cellstr(ls('VPPG*'));
    
    ssdir = 'PDT_ss_design_LAcl_01';
    
    for ii = 1:length(all_subf)
        try
            %for ii = 4:length(all_subf)
            cd(base_dir_pl);
            cd(all_subf{ii})
            cd('MRT\NIFTI\results')
            agk_review_design_matrices(cd,ssdir)
        catch
            disp(['Failed to display ' ssdir ' of ' all_subf{ii}])
        end
    end
    
end

%% ############### ss-LEVEL add cons ##############
if ss_addcons
    cd(base_dir_pl);
    all_subf         = cellstr(ls('VPPG*'));
    ss_name          = 'PDT_ss_design_01'; % name of ss design to add cons to
    cur_templ        = 'C:\Users\genaucka\Google Drive\Library\MATLAB\PDT\ss_models\con_man_template.mat';
    tcon.names       = {'pic.on.gam_gr_pic.on.pos','pic.on.gam_kl_pic.on.pos','pic.on.gam_gr_pic.on.neg','pic.on.gam_kl_pic.on.neg'};
    tcon.codes       = {[0 1 0 -1], [0 -1 0 1], [0 1 -1 0], [0 -1 1 0]};
    del              = 0
    
    %FOR PDT_ss_design_00
    %tcon.names       = {'PICGAM>PIC','PICGAM.loss>PICGAMOPT.loss'};
    %tcon.codes       = {[-1 1], [0 0 0 1 0 0 0 -1]};
    
    % ADDING CONS TO THE SPM.mat of SS ANALYSIS
    for ii = 1:length(all_subf)
        cd(base_dir_pl);
        ss_addcon_feedback{ii} = agk_add_cons(base_dir_pl,all_subf{ii},ss_name, cur_templ,tcon,del);
    end
    
end

%% ############### ss-LEVEL new cons ##############
if ss_newcons
    cd(base_dir_pl);
    all_subf         = cellstr(ls('VPPG*'));
    ss_name          = 'PDT_ss_design_LAcl_01'; % name of ss design to add cons to
    cur_templ        = 'C:\Users\genaucka\Google Drive\Library\MATLAB\PDT\ss_models\con_man_template.mat';
    del              = 1;
    
    %FOR PDT_ss_design_00
    %tcon.names       = {'PICGAM>PIC','PICGAM.loss>PICGAMOPT.loss'};
    %tcon.codes       = {[-1 1], [0 0 0 1 0 0 0 -1]};
    
    % ADDING CONS TO THE SPM.mat of SS ANALYSIS
    for ii = 2:length(all_subf)      
        cd(base_dir_pl);
        ss_addcon_feedback{ii} = agk_make_cons_ss(base_dir_pl,all_subf{ii},ss_name, cur_templ,del);
    end
end

%% ############### extract ss-eigvar ##############
if ss_extr
    % make ROIs
    % names table 2:
    names_ROIs       = {'r_caudate_head','l_medial_frontal_gyrus','l_posterior_cingulate', 'l_cingulate_gyrus', 'r_claustrum', 'l_middle_temporal_gyrus',  'l_precuneus', 'l_thalamus' , 'l_fusiform_gyrus', 'l_anterior_cingulate', 'r_inferior_occipital_gyrus', 'l_thalamus', 'l_parahippocampal_gyrus', 'l_middle_frontal_gyrus'};
    % coord_ROI table2:
    coord_ROIs       = {[8 4 -2], [-2 52 6],[-4 -50 18],[-4 -38 30], [32 10 -6], [-58 -30 -4], [-2 -62 60], [-20 -22 12], [-50 -42 -12], [-2 24 16], [28 -88 -10], [-2 -6 10], [-20 -12 -16], [-20 22 58]};
    
    % names table 3:
    % names_ROIs  = {'t3_l_posterior_cingulate', 't3_l_superior_temporal_gyrus', 't3_r_precuneus', 't3_l_precuneus'};
    % coord_ROI table3:
    % coord_ROIs = {[-6 -36 26], [-52 10 -14], [34 -74 36], [-12 -58 26]};
    
    cur_SPM          = 'E:\Daten\VPPG\MRT\MRT\VPPG0104\MRT\NIFTI\results\PDT_ss_design_01\SPM.mat';
    for kk = 1:length(names_ROIs)
        cur_coord        = coord_ROIs{kk};
        cur_name         = cellstr(names_ROIs{kk});
        cur_dest         = 'C:\Users\genaucka\Google Drive\Library\MATLAB\PDT\sl\ROIs';
        cd (cur_dest)
        create_sphere_image(cur_SPM,cur_coord,cur_name,repmat(20,3,1));
    end
    
    % calculate a combined cue reactivity map
    name_out = 'table2.nii';
    names_in = {};
    for gg = 1:length(coord_ROIs) 
        names_in{gg} = [names_ROIs{gg} '_mask.nii'];
    end
    f = 'sum(X)';
    spm_imcalc(names_in,name_out,f,{1})

    
    for kk = 1:1
        % prep extr.
        cd(base_dir_pl);
        %cur_ROI          = 'C:\Users\genaucka\Google Drive\Library\MATLAB\LA\2L_Analysis\ROI\accumbens\accumbens.nii';
        cur_ROI          = ['C:\Users\genaucka\Google Drive\Library\MATLAB\PDT\sl\ROIs\table3.nii']; % path to ROI where to extract
        all_subf         = cellstr(ls('VPPG*'));
        ss_name          = 'PDT_ss_design_01'; % name of ss design for extr.
        
        % EXTRACTING BOLD RESPONSE PER TRIAL DURING PIC-PHASE AT ROI
        for ii = 1:1
            cd(base_dir_pl);
            cur_sub = all_subf{ii};
            %ss_extr_feedback{ii} = agk_get_pic_mri_scores(cur_sub,ss_name,cur_ROI);
            ss_extr_feedback{ii} = agk_get_pic_mri_timeseries(cur_sub,ss_name,cur_ROI);
        end
    end
    
    % deconvolving
    cur_ts  = cell2mat(ss_extr_feedback{ii}{5});
    cur_mn_tsg = mean(cur_ts(ss_extr_feedback{1}{4} == 1,:));
    plot(cur_mn_tsg)
    cur_hrf = spm_hrf(1);
    
    v          = cur_hrf;
    x          = (1:1:length(v)); % original coordinates of sampling
    xq         = linspace(x(1), x(end), (length(x))*1000); % TR of 2s to 1000Hz, interpolated sampling points;
    Y_1        = interp1(x,v,xq);
    cur_hrf    = zscore(Y_1(1:16000));
    
    [q,r] = deconv(zscore(cur_mn_tsg(2:16000)),cur_hrf(2:16000));
    
    
    % DONT FORGET TO MOVE EXTRACTS AFTER EXTRACTING! (USING BEHAV MOVE)
    
end