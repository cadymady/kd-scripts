% agk_make_PDT_ss_pspml_only_pic_tworuns
% for one subject
% makes the ss design of the PDT task based on the Tom et al. (2007) paper
% pictures here are put in as categorical vectors
% to allow a cue reactivity contrast (between categories modeling)

% INPUT
% cur_sub       : name of current subject
% which_folders : name of folder with niftis of task (pciks first) ???
% aggr          : aggregation level gambling matrix (we usually take 3)
% run_it        : should the job run? (will be only save if not)
% acc_rec       : should accept be included as a factor in model (def: no)
% expl_mask     : usually use a gray matter mask here
% cur_tmpl      : the ss job template to be used here, cell, first field is
%                 the one without accept reject as factor and second with
% ow            : overwrite exisiting results files? (0,1)
% do_center     : should be 1
% fixed_dur     : e.g. 0 if stick functions wanted; or another scalar,
% fixed_dur     : if empty, then duration until next onset will be used

% OUTPUT
% error_message: success message or explanation what hasn't worked

% IMPLICIT OUTPUT
% writes the ss model into the results_ssdir in the subject folder

%THIS MODEL ('pspml')
% 4 onsets pic: PIC_NEU, PIC_GAM, PIC_NEG, PIC_POS
% 1 onset gamble with 3 modulators: gain, abs(loss), ed

% DETAILS
% preproc_job is needed to get the slice timing parameters
% (to create time bins), will be looked for in the nifti folder provided
% will look for "preprocessing..." pattern
% abs values of losses are used
% nifti data need to be SPM12 preprocessed
% takes in an aggregation factor: e.g. 1 means no change; 2 means, the 12
% by 12 matrix will be aggregated to a 6-by-6 matrix, and so on
% the picture category is modeled by three dummy variables (i.e. all 0 is
% neutral)
% PIC+GAM+OPT and FEEDBACK are not parametrically modulated (leads to too
% much intercorrelation in design matrix

% AUTHORSHIP
% author: Alexander Genauck
% revision: Katharina Dücker
% date  : 26.10.2017
% email : alexander.genauck@charite.de

function [error_message] = agk_make_PDT_ss_pspml_only_pic_tworuns_new(cur_sub,cur_tmpl,aggr, ...
    run_it,acc_rec,expl_mask,ow,do_center,fixed_dur, physio_inc)
try
    %% PREPARATIONS
    warning('Still have to include a censoring regressor for missing trials!')
    % change into subs directory
    root_dir = pwd;
    cd(cur_sub)
    
    % name of this analysis determines the results directory
    if physio_inc
        results_ssdir = 'agk_make_PDT_ss_pspml_only_pic_tworuns_new_missing_physio';
    else
        results_ssdir = 'agk_make_PDT_ss_pspml_only_pic_tworuns_new_missing';
    end
    
    % messaging
    disp(['Estimating subject ' cur_sub ' using ' results_ssdir])
    %     disp(cur_sub.name);
    %     disp(['First level model' results_ssdir]);
    
    % gain and loss ranges
    gain_min = 14;
    gain_max = 36;
    loss_max = 18;
    loss_min = 7;
    
    % for euclidean distance
    vec = [2;1;0];   % slope vector [gain; loss; 0]
    sp  = [26;13;0]; % sp support point; point on the diagonal [gain; loss; 0]
    
    % calculate the aggregated possible values
    cur_gains = (gain_min:2:gain_max);
    cur_losss = loss_min:loss_max    ;
    [vg,osg,nsg]  = agk_downsample_steps(cur_gains,aggr);
    gain_min  = min(vg);
    gain_max  = max(vg);
    [vl,osl,nsl]  = agk_downsample_steps(cur_losss,aggr);
    loss_min  = min(vl);
    loss_max  = max(vl);
    
    %%%%%%%% get behav data  %FC: moved up because we need some parameters of
    % P structure before filling the batch
    behav_data = 0;
    try
        cd('Behav\PDT')
        % special case: VPPG0115 -> take second behavioral file
        load(ls('P*'))
        behav_data = 1;
    catch
        behav_data = 0;
        error_message = [results_ssdir ' ' cur_sub ' No P struct found.'];
        cd(root_dir)
        return
    end
    tworuns   = isfield(P.t,'triggerscannertpostpause');
    %% PREPARE FILLING BATCH
    % preparing the design information variables
    
    % tworun corrections
    %check if this is two runs and get onset correction
    
    if tworuns
        % new starting point in time for second run
        time_corr_second_run_fixed = (P.t.triggerscannert0 - ...
            P.t.triggerscannertpostpause);
        % length of each run
        length_run = length(P.cur.choice)/2;     
        num_run = 2;
       % special case: VPPG0115
      
    else
        time_corr_second_run_fixed = 0;
        length_run = length(P.cur.choice);
        num_run = 1;
    end
    

        names       = {'Pic.neu','Pic.gam','Pic.neg','Pic.pos', 'Pic.gam.on', ...           
            'Pic.gam.opt.on'};
        if acc_rec
            names{length(names)+1} = 'acceptance';
        end
        onsets          = cell(1,length(names));
        durations       = cell(1,length(names));
        
        pmod            = struct('name',{},'poly',{},'param',{});
        pmod(5).name{1} = 'gain';
        % trial: add pmod.poly and name
        pmod(5).poly{1} = 1;
        pmod(5).name{2} = 'loss';
        pmod(5).poly{2} = 1;
        pmod(5).name{3} = 'ed';
        pmod(5).poly{3} = 1;
        % prep params
        for kk = 1:length(pmod)
            params{kk,1} = cell(1,length(pmod(kk).name));
        end
        % creating upper level onsets, durations, pmod cells such that
        % it's easier to save them for each run (see below!)
        names_cell       = cell(1,num_run);
        onsets_cell      = cell(1,num_run);
        durations_cell   = cell(1,num_run);
        pmod_cell        = cell(1,num_run);
        params_cell      = cell(1,num_run);
    
    for n = 1:num_run
        names_cell{n}     = names;
        onsets_cell{n}    = onsets;
        durations_cell{n} = durations;
        pmod_cell{n}      = pmod;
        params_cell{n}    = params; 
    end
    
    % getting the gain mean and the loss mean
    
    for r1 = 1:length(P.cur.choice)
        all_gains(r1) = str2double(cell2mat(P.gain.strings(P.cur.gamble{r1}(1))));
        all_losss(r1) = str2double(cell2mat(P.loss.strings(P.cur.gamble{r1}(2))));
        
    end

    
    % calculate mean gain and loss
    mean_gain = mean(all_gains);
    mean_loss = mean(all_losss);
    
    % mean centering?
    if do_center == 0
        mean_gain = 0;
        mean_loss = 0;
    end
    
    %% fill onsets, pmod p's
    %     try
    % special case: VPPG0115
    
    for ii = 1 : length(P.cur.choice)
        % two run correction
        
        
        if tworuns
            if ii <= length_run
                % first run
                nn = 1;
                time_corr_second_run = 0;
                % run length correction
                rlc = 0;
            else
                % second run
                nn = 2;
                time_corr_second_run = time_corr_second_run_fixed;
                % run length correction
                rlc = length_run;
                if strcmp(cur_sub, 'VPPG0115')
                    cd(root_dir)
                    cd(cur_sub)
                    cd('Behav\PDT')
                    cd(p_structs{2})
                    load(ls('P*'))
                end
            end
        else
            % single session (like a long first run)
            time_corr_second_run = 0;
            nn = 1;
            % run length correction
            rlc = 0;
        end
        
        % MISSING AS MULTIPLE REGRESSOR
        if ~(P.cur.choice(ii) > 0 && P.cur.choice(ii) < 5) && floor((P.t.cur_trial.stimLA_off(ii)+time_corr_second_run)/2) < size(tmp,1)
           % find limits of the trial -> don't cut out too much.
           % (max(round..)
           trial{nn}{ii-rlc} = [ceil((P.t.cur_trial.iti_on(ii)+time_corr_second_run)/2), floor((P.t.cur_trial.stimLA_off(ii)+time_corr_second_run)/2)];
           missing{nn}(trial{nn}{ii-rlc}(1):trial{nn}{ii-rlc}(2)) = 1;

            continue
        end
        
        % PIC ON
        switch P.cur.cat(ii)
            case 1
                ons_code = 2;
            case 2
                ons_code = 3;
            case 3
                ons_code = 4;
            case 6
                ons_code = 1;
        end
        
        onsets_cell{nn}(ons_code) = {[cell2mat(onsets_cell{nn}(ons_code)) P.t.cur_trial.stim_on(ii) + time_corr_second_run]}; %changed nn to num_run
        if ~isempty(fixed_dur)
            
            durations_cell{nn}(ons_code) = {[cell2mat(durations_cell{nn}(ons_code)) fixed_dur]};
        else
            
            durations_cell{nn}(ons_code) = {[cell2mat(durations_cell{nn}(ons_code)) (P.t.cur_trial.stimLA_watch(ii) - P.t.cur_trial.stim_on(ii))]};
        end
        
        % PIC PLUS GAMBLE (WITHOUT CHOICE OPTIONS)
       
        onsets_cell{nn}(5) = {[cell2mat(onsets_cell{nn}(5)) P.t.cur_trial.stimLA_watch(ii) + time_corr_second_run]}; %{[cell2mat(onsets_cell{nn}(5)) P.t.cur_trial.stimLA_on(ii) + time_corr_second_run]};
        
        if ~isempty (fixed_dur)
            durations_cell{nn}(5) = {[cell2mat(durations_cell{nn}(5)) fixed_dur]};
            
        else
            durations_cell{nn}(5) = {[cell2mat(durations_cell{nn}(5)) P.t.cur_trial.stimLA_on(ii) - P.t.cur_trial.stimLA_watch(ii)]};
        end
        
        % gain
        cur_gain     = agk_recode(str2double(cell2mat(P.gain.strings(P.cur.gamble{ii}(1)))),osg,nsg) - mean_gain;
        gain_orig    = cur_gain + mean_gain;
        % loss (here changed loss to abs. loss)
        cur_loss     = agk_recode(abs(str2double(cell2mat(P.loss.strings(P.cur.gamble{ii}(2))))),osl,nsl) - abs(mean_loss);
        loss_orig    = cur_loss+abs(mean_loss);
        % ed
        cur_point    = [gain_orig;loss_orig;0];
        ed           = agk_get_ed(cur_point,sp,vec);
        
        % fill into params cell
        params_cell{nn}{5}(1) = {[cell2mat(params_cell{nn}{5}(1)) cur_gain]};
        params_cell{nn}{5}(2) = {[cell2mat(params_cell{nn}{5}(2)) cur_loss]};
        params_cell{nn}{5}(3) = {[cell2mat(params_cell{nn}{5}(3)) ed]};
        
        % PIC PLUS GAMBLE ON (WITH CHOICE OPTIONS)
        
        onsets_cell{nn}(6) = {[cell2mat(onsets_cell{nn}(6)) P.t.cur_trial.stimLA_on(ii) + time_corr_second_run]}; %{[cell2mat(onsets_cell{nn}(5)) P.t.cur_trial.stimLA_on(ii) + time_corr_second_run]};
        
        if ~isempty (fixed_dur)
            durations_cell{nn}(6) = {[cell2mat(durations_cell{nn}(6)) fixed_dur]};
            
        else
            durations_cell{nn}(6) = {[cell2mat(durations_cell{nn}(6)) P.t.cur_trial.stimLA_off(ii) - P.t.cur_trial.stimLA_on(ii)]};
        end
        
         % fill the pmods with params
         
         if ~isempty(pmod_cell{nn}(kk).name)
             pmod_cell{nn}(5).param = params_cell{nn}{5};
         end
         
         pmod_cell{nn}(6).name = [];
         pmod_cell{nn}(6).poly = [];
         pmod_cell{nn}(6).param = [];
         
    end

     % Concatenate Missing and multiple regressor .txt file
    if physio_inc && exist(fullfile(root_dir, cur_sub, 'Physio'))
        cd(fullfile(root_dir, cur_sub, 'Physio\results'))
        corr_m_r = cellstr(ls('corr*'));
        if isempty(corr_m_r{1})
            for mm = 1:num_run
                file_mult = load(mult_reg{mm});
                file_mult_corr = horzcat(file_mult, missing{mm});
                save(['corr_multiple_regressors_run_',num2str(mm),'.txt'],'file_mult_corr', '-ascii')
            end
        end
        
    elseif physio_inc == 0 || ~exist(fullfile(root_dir, cur_sub, 'Physio'))
        for mm = 1:num_run
            disp(['No multiple regressor file for subject ', cur_sub,'. Missing trials are added to rp file'])
            cd(fullfile(root_dir, cur_sub, 'MRT\NIFTI\PDT',cur_fold{mm}))
            corr_rp = cellstr(ls('corr*'));
            if isempty(corr_rp{1})
                rp = load(ls('rp*'));
                rp_corr = horzcat(rp, missing{mm});
                save(['corr_rp_aepi_run',num2str(mm),'_MoCoAP_MoCo_00001.txt'], 'rp_corr', '-ascii')
            end
        end
    end
  
    
    %%%%%% % get the EPIs of this subject
    cd(root_dir)
    cd(cur_sub)
    try
        cd('MRT\NIFTI\PDT');
    catch
        disp([results_ssdir ' ' cur_sub ' no NIFTI dir. Skipping.'])
        error_message =  [results_ssdir ' ' cur_sub '  no NIFTI dir. Skipped.']; % ??? !!!
        return
    end
    
    
  
    found_epi = 0;
    cur_epi_dirs  = cellstr(ls('*_epi*'));
    cur_MoCo_dirs = cellstr(ls('*_MoCo*'));
    
    if length(cur_epi_dirs{1}) > 0 % check if we're working with epis
        found_epi = 1;
        cur_fold = cur_epi_dirs;
    elseif length(cur_MoCo_dirs{1}) > 0 % check if we're working with MoCo
        
        found_epi = 1;
        cur_fold = cur_MoCo_dirs;       % store current folder name
    end
    
    cd(cur_fold{1})
    load(ls('Preprocessing*')); % for microtime resolution
    preproc_job = matlabbatch;
    cur_nslices   = preproc_job{1}.spm.temporal.st.nslices;  % store slices etc in batch
    cur_refslice  = preproc_job{1}.spm.temporal.st.refslice;
    
    % load multiple regressor (physio corrected) .txt file
    
    cd(root_dir);  % add path
    cd(cur_sub);
    try
        cd('Physio\results')
        mult_reg = cellstr(ls('multiple_regressors*'));
    catch
        disp('No Physio folder for this subject')
    end
   
    % load template for ss model
    if acc_rec == 1
        load(cur_tmpl{2})
    elseif acc_rec == 0
        if tworuns
            load(cur_tmpl{3})
        else
            load(cur_tmpl{1})
        end
    end
    
    %% fill in number of slices and multiple regressors/rp file
    % CHANGE kd: in loop over number of runs
    
        
    if ~tworuns
        cd(root_dir);  % add path
        cd(cur_sub);
        cd('MRT/NIFTI/PDT')
        cd(cur_fold{1})
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = [];
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('FPList',pwd,'swuaepi'));
        
        % if physio should be included, go into Physio results folder
        % and load .txt files from physio correction
        if physio_inc && exist(fullfile(root_dir, cur_sub, 'Physio'))
            % create multi_reg subfield
           matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(root_dir, cur_sub, 'Physio/results',mult_reg{1})};
            
        else
            % else load rp file from preprocession
           matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(spm_select('FPList',pwd,'^rp'));
           results_ssdir = 'agk_make_PDT_ss_pspml_only_pic_tworuns_new_missing';
        end
        
        if found_epi == 0 % incase no scans were found: throw away
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {''};
        end
        
        
        cd(root_dir)
        cd(cur_sub)
        cd('MRT/NIFTI/PDT')
    else
        for nr = 1:num_run
            cd(root_dir);  % add path
            cd(cur_sub);
            cd('MRT/NIFTI/PDT')
            cd(cur_fold{nr})
            matlabbatch{1}.spm.stats.fmri_spec.sess(nr).scans = [];
            matlabbatch{1}.spm.stats.fmri_spec.sess(nr).scans = cellstr(spm_select('FPList',pwd,'swuaepi'));
            
            % if physio should be included, go into Physio results folder
            % and load .txt files from physio correction
            if physio_inc && exist(fullfile(root_dir, cur_sub, 'Physio'))
                % create multi_reg subfield
                matlabbatch{1}.spm.stats.fmri_spec.sess(nr).multi_reg = [];
                matlabbatch{1}.spm.stats.fmri_spec.sess(nr).multi_reg = {fullfile(root_dir, cur_sub, 'Physio/results',mult_reg{nr})};
                
            else
                % else load rp file from preprocession
            
                matlabbatch{1}.spm.stats.fmri_spec.sess(nr).multi_reg = [];
                matlabbatch{1}.spm.stats.fmri_spec.sess(nr).multi_reg = cellstr(spm_select('FPList',pwd,'^rp'));
            end
            
            
            if found_epi == 0 % incase no scans were found: throw away
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {''};
            end
            
                      
        end
    end
    
   
    %%%%%%%% create the results dir
    mkdir(pwd,'results')
    cd('results')
    
    if exist([pwd filesep results_ssdir])
        if ow == 1
            cmd_rmdir([pwd filesep results_ssdir])
        elseif ow == 0
            disp('Results dir already present. Overwrite not allowed, I will skip this subject.')
            error_message =  [results_ssdir ' ' cur_sub ' results dir already present. Skipped.'];
            return
        end
    end
    
    agk_mkdir_ex(pwd,results_ssdir)
    cd(results_ssdir)
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(pwd);
    
    if ~tworuns
        if behav_data == 1
            names = names_cell{1};
            onsets = onsets_cell{1};
            durations = durations_cell{1};
            
            pmod = pmod_cell{1};
            save('mult_cond.mat','names','onsets','durations','pmod');
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = cellstr(spm_select('FPList',pwd,'mult_cond.mat'));
        end
    else % in case of two sessions
        if behav_data == 1
            for n = 1:num_run
                names     = names_cell{n};
                onsets    = onsets_cell{n};
                
                durations = durations_cell{n};
                
                pmod      = pmod_cell{n};
                save(['mult_cond_',num2str(n),'.mat'],'names','onsets','durations','pmod');
                matlabbatch{1}.spm.stats.fmri_spec.sess(n).multi = cellstr(spm_select('FPList',pwd,['mult_cond_',num2str(n),'.mat']));
            end
        end
    end
    
    % fill in microtime resolution
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = cur_nslices;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = cur_refslice;
    % specifiy an explicit mask, which voxels shall only be analysed
  %  matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.mask = cellstr(expl_mask);
    % specify TR-time
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.0; %% CHANGE!
    
    %%%%%%%%%%%TO DO: CHANGE THIS adjusting to new contrasts    % contrast manager
    con.contrastNames = {'pic.on','pic.on.gam','pic.on.neg','pic.on.pos','picgam.on.gain','picgam.on.loss','picgam.on.ed', 'picgam.opt.on'};  %changed FC
    con.contrastType = cellstr(repmat('t',length(con.contrastNames),1))';
    con.contrastWeights = agk_basic_t_cons(length(con.contrastType));
    if tworuns
        con.contrastWeights = agk_basic_t_cons_2sess(length(con.contrastType),length(con.contrastType));
    end
    con.contrastRep = cellstr(repmat('none',length(con.contrastNames),1))';
    
    % dependencies: estimation to model specification
    %               contrast manager to model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

    for jj = 1:length(con.contrastType)
        matlabbatch{3}.spm.stats.con.consess{jj}.tcon.name    = con.contrastNames{jj};
        matlabbatch{3}.spm.stats.con.consess{jj}.tcon.convec  = con.contrastWeights{jj};
        matlabbatch{3}.spm.stats.con.consess{jj}.tcon.sessrep = con.contrastRep{jj};
    end
    
    if behav_data == 1
        save('design.mat','matlabbatch');
    end
    
    try
        if run_it == 1
            spm_jobman('run',matlabbatch);
        end
        error_message = [results_ssdir ' ' cur_sub ' Estimation successfull.'];
    catch MExc
        error_message = {MExc,[results_ssdir ' ' cur_sub ' Estimation not successfull.']};
        cd(root_dir)
        return
    end
    cd(root_dir);

catch MExc
    disp('Something went wrong');
    error_message = MExc;
end
