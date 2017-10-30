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

function [error_message] = agk_make_PDT_ss_pspml_only_pic_tworuns_accint(cur_sub,cur_tmpl,aggr, ...
    run_it,acc_rec,expl_mask,ow,do_center,fixed_dur, physio_inc)
% try
    %% PREPARATIONS
    %warning('Still have to include a censoring regressor for missing trials!')
    % change into subs directory
    root_dir = pwd;
    cd(cur_sub)
    
    % name of this analysis determines the results directory
    if physio_inc
        results_ssdir = 'agk_make_PDT_ss_pspml_only_pic_tworuns_acc_missing_physio';
    else
        results_ssdir = 'agk_make_PDT_ss_pspml_only_pic_tworuns_acc_missing';
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
    % get the ed as well
   
    
    %%%%%%%% get behav data  %FC: moved up because we need some parameters of
    % P structure before filling the batch
    behav_data = 0;
    try
        cd('Behav\PDT')
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
    % check if this is two runs and get onset correction
    
    if tworuns
        % new starting point in time for second run
        time_corr_second_run_fixed = (P.t.triggerscannert0 - ...
            P.t.triggerscannertpostpause);
        % length of each run
        length_run = length(P.cur.choice)/2;
        num_run = 2;
        
    else
        time_corr_second_run_fixed = 0;
        length_run = length(P.cur.choice);
        num_run = 1;
    end
    
    edabs = [];
    for ff = 1:length(P.cur.choice)
        cur_gain  = str2double(cell2mat(P.gain.strings(P.cur.gamble{ff}(1))));
        cur_loss  = str2double(cell2mat(P.loss.strings(P.cur.gamble{ff}(2))));
        cur_point = [cur_gain;cur_loss;0];
        edabs(ff) = agk_get_ed(cur_point,sp,vec);
    end
    
    % downsample also edabs
    edabs           = round(edabs);
    [edabs,ose,nse] = agk_downsample_steps(edabs,10);
    
    % FIRST PHASE: picture on
    names       = {'Pic.gam','Pic.neg','Pic.pos'};
    
    % SECOND PHASE: Pic + gamble
    
        % add the gain,loss,edabs Pic.gam onsets
        % which gain should be dummy coded? 1: lower, 2: upper
        gain_steps = unique(nsg);
        gn_stps_on = [];
        for kk = 2:length(gain_steps)
            cur_pg_name                      = ['Pic.Gam.gain' ...
                num2str(gain_steps(kk))];
            names{length(names)+1}           = cur_pg_name;
            gn_stps_on(length(gn_stps_on)+1) = length(names);
        end
        gn_stps_on = [0, gn_stps_on]; % dummy code
        
        loss_steps = unique(nsl);
        ls_stps_on = [];
        
        for kk = 2:length(loss_steps)
            cur_pg_name                      = ['Pic.Gam.loss' ...
                num2str(loss_steps(kk))];
            names{length(names)+1}           = cur_pg_name;
            ls_stps_on(length(ls_stps_on)+1) = length(names);
        end
        ls_stps_on = [0, ls_stps_on];
        
        
        edabs_steps = unique(nse);
        ed_stps_on = [];
        
        for kk = 2:length(edabs_steps)
            cur_pg_name                      = ['Pic.Gam.edabs' ...
                num2str(edabs_steps(kk))]; 
            names{length(names)+1}           = cur_pg_name;
            ed_stps_on(length(ed_stps_on)+1) = length(names);
        end
        ed_stps_on = [0, ed_stps_on];
        
        
        % add acceptance 
        names{length(names)+1}                  = 'acceptance';
        % THIRD PHASE: Pic + gamble + opt
        names{length(names)+1} = 'Pic.Gam.Opt';
        
        % prepare cells for onsets and durations
        onsets          = cell(1,length(names));
        durations       = cell(1,length(names));
    

    % pmod now only used for acceptance rate for pic onsets
    pmod            = struct('name',{},'poly',{},'param',{});
    for ii = 1:3                    % 3 because we have 3 pic cats
        pmod(ii).name{1} = 'acceptance';
        pmod(ii).poly{1} = 1;
    end

    for pp = 5:length(names)
        pmod(pp).name = {};
    end
    
    % pmod for acceptance by picture category
    pmod(length(names)-1).name{1} = 'gam';
    pmod(length(names)-1).name{2} = 'neg';
    pmod(length(names)-1).name{3} = 'pos';
    
    for pl = 1:length(pmod(length(names)-1).name)
        pmod(length(names)-1).poly{pl} = 1;
    end
    
    % collecting the parametric modulators (for accept/reject)
    ps = cell(1,num_run);
    for pp = 1:num_run
        ps{pp} = cell(1,length(pmod));
        ps{pp}{length(pmod)-1} = cell(1,length(pmod(length(names)-1).name));
    end
    
    % creating upper level onsets, durations, pmod cells such that
    % it's easier to save them for each run (see below!)
    names_cell       = cell(1,num_run);
    onsets_cell      = cell(1,num_run);
    durations_cell   = cell(1,num_run);
    pmod_cell        = cell(1,num_run);
    
    
    for n = 1:num_run
        names_cell{n}     = names;
        onsets_cell{n}    = onsets;
        durations_cell{n} = durations;
        pmod_cell{n}      = pmod;
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
    
    % load multiple regressor (physio corrected) .txt file
    % load epis
    cd(fullfile(root_dir, cur_sub, 'MRT\NIFTI\PDT'))
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
    
    cd(root_dir);  % add path
    cd(cur_sub);
    try
        cd('Physio\results')
        mult_reg = cellstr(ls('multiple_regressors*'));
        for ml = 1:num_run
            tmp = load(mult_reg{ml});
            missing{ml} = zeros(size(tmp,1),1);
        end
    catch
        disp('No Physio folder for this subject')
        for ml = 1:num_run
            cd(fullfile(root_dir, cur_sub, 'MRT\NIFTI\PDT', cur_fold{ml}))
            tmp = (load(ls('rp*')));
            missing{ml} = zeros(size(tmp,1),1);
        end
    end
    
    
    %% fill onsets, pmod p's
    
    for ii = 1: length(P.cur.choice)
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
                ons_code = 1;
            case 2
                ons_code = 2;
            case 3
                ons_code = 3;
            case 6
                ons_code = nan;
        end
        
         
        
        if ~isnan(ons_code)
            onsets_cell{nn}(ons_code) = {[cell2mat(onsets_cell{nn}(ons_code)) P.t.cur_trial.stim_on(ii) + time_corr_second_run]}; %changed nn to num_run
            if ~isempty(fixed_dur)
                
                durations_cell{nn}(ons_code) = {[cell2mat(durations_cell{nn}(ons_code)) fixed_dur]};
            else
                
                durations_cell{nn}(ons_code) = {[cell2mat(durations_cell{nn}(ons_code)) (P.t.cur_trial.stimLA_watch(ii) - P.t.cur_trial.stim_on(ii))]};
            end
            % pmod for PIC ON (acceptance) = acceptance induced by pic category
            
            cur_acc = agk_recode(P.cur.choice(ii),[1,2,3,4,5],[1,1,0,0,5]);
            ps{nn}{ons_code} = [ps{nn}{ons_code} cur_acc];
            pmod_cell{nn}(ons_code).param = ps{nn}(ons_code);
        end
        
        % PIC PLUS GAMBLE (WITHOUT CHOICE OPTIONS)
        
        % recode current gain/loss to steps
        onsets_to_be_filled = [];
        cur_gain  = str2double(cell2mat(P.gain.strings(P.cur.gamble{ii}(1))));
        cur_gain  = agk_recode(cur_gain,osg,nsg);
        cur_gn_on = agk_recode(cur_gain,gain_steps,gn_stps_on);
        cur_loss  = abs(str2double(cell2mat(P.loss.strings(P.cur.gamble{ii}(2)))));
        cur_loss  = agk_recode(cur_loss,osl,nsl);
        cur_ls_on = agk_recode(cur_loss,loss_steps,ls_stps_on);
        cur_ed_on = agk_recode(edabs(ii),edabs_steps,ed_stps_on);
        
        % onsets for gain, loss, ed
        if cur_gn_on
            onsets_to_be_filled(1) = cur_gn_on;
        else
            onsets_to_be_filled(1) = nan;
        end
        
        if cur_ls_on
            onsets_to_be_filled(2) = cur_ls_on;
        else
            onsets_to_be_filled(2) = nan;
        end
        
        if cur_ed_on
            onsets_to_be_filled(3) = cur_ed_on;
        else
            onsets_to_be_filled(3) = nan;
        end

        
        for jj = 1:length(onsets_to_be_filled)
            if ~isnan(onsets_to_be_filled(jj))
                onsets_cell{nn}(onsets_to_be_filled(jj)) = {[cell2mat(onsets_cell{nn}(onsets_to_be_filled(jj))) P.t.cur_trial.stimLA_watch(ii)+ time_corr_second_run]};
                if ~isempty(fixed_dur)
                    durations_cell{nn}(onsets_to_be_filled(jj)) = {[cell2mat(durations_cell{nn}(onsets_to_be_filled(jj))) fixed_dur]};
                else
                    durations_cell{nn}(onsets_to_be_filled(jj)) = {[cell2mat(durations_cell{nn}(onsets_to_be_filled(jj))) (P.t.cur_trial.stimLA_on(ii) - P.t.cur_trial.stimLA_watch(ii))]};
                end
            end
        end
        
         % ACCEPTANCE (ONSET WITH OPT ON) will be set if accept == 1
        if cur_acc
            cur_cat = agk_recode(P.cur.cat(ii),[1,2,3,6],[1,2,3,0]);
%             disp('cur_cat is...')
%             disp(num2str(cur_cat))
            if cur_cat
                onsets_cell{nn}(length(onsets)-1) = {[cell2mat(onsets_cell{nn}(length(onsets)-1)) P.t.cur_trial.stimLA_watch(ii)+ time_corr_second_run]};
                if ~isempty (fixed_dur)
                    durations_cell{nn}(length(onsets)-1) = {[cell2mat(durations_cell{nn}(length(onsets)-1)) fixed_dur]};
                else
                    durations_cell{nn}(length(onsets)-1) = {[cell2mat(durations_cell{nn}(length(onsets)-1)) (P.t.cur_trial.stimLA_on(ii) - P.t.cur_trial.stimLA_watch(ii))]};
                end
                
                ps{nn}{length(onsets)-1}{1} = [ps{nn}{length(onsets)-1}{1} 0];
                ps{nn}{length(onsets)-1}{2} = [ps{nn}{length(onsets)-1}{2} 0];
                ps{nn}{length(onsets)-1}{3} = [ps{nn}{length(onsets)-1}{3} 0];
                
                ps{nn}{length(onsets)-1}{cur_cat}(length(ps{nn}{length(onsets)-1}{cur_cat})) = 1;
                pmod_cell{nn}(length(onsets)-1).param = ps{nn}{length(onsets)-1};
%                 
%                 disp('pmod params are ...')
                pmod_cell{nn}(length(onsets)-1).param{1};
                pmod_cell{nn}(length(onsets)-1).param{2};
                pmod_cell{nn}(length(onsets)-1).param{3}; 
            end
        end
        
        % PIC PLUS GAMBLE WITH OPTIONS 
        onsets_cell{nn}(length(onsets)) = {[cell2mat(onsets_cell{nn}(length(onsets))) P.t.cur_trial.stimLA_on(ii)+ time_corr_second_run]};
        if ~isempty (fixed_dur)
            durations_cell{nn}(length(onsets)) = {[cell2mat(durations_cell{nn}(length(onsets))) fixed_dur]};
        else
            durations_cell{nn}(length(onsets)) = {[cell2mat(durations_cell{nn}(length(onsets))) (P.t.cur_trial.stimLA_off(ii) - P.t.cur_trial.stimLA_on(ii))]};
        end

        
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
    
    cd(cur_fold{1})
    load(ls('Preprocessing*')); % for microtime resolution
    preproc_job = matlabbatch;
    cur_nslices   = preproc_job{1}.spm.temporal.st.nslices;  % store slices etc in batch
    cur_refslice  = preproc_job{1}.spm.temporal.st.refslice;
    

   
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
        if physio_inc && exist(fullfile(root_dir, cur_sub, 'Physio\results'))
            % create multi_reg subfield
           mult_reg = cellstr(ls('corr*'));
           matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(root_dir, cur_sub, 'Physio/results',mult_reg{1})};
            
        else
            % else load rp file from preprocession
           matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(spm_select('FPList',pwd,'^corr'));
           results_ssdir = 'agk_make_PDT_ss_pspml_only_pic_tworuns_acc_missing';
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
                cd(fullfile(root_dir, cur_sub, 'Physio\results'))
                % create multi_reg subfield
                mult_reg = cellstr(ls('corr*'));
                matlabbatch{1}.spm.stats.fmri_spec.sess(nr).multi_reg = [];
                matlabbatch{1}.spm.stats.fmri_spec.sess(nr).multi_reg = {fullfile(root_dir, cur_sub, 'Physio/results',mult_reg{nr})};
                
            else
                % else load rp file from preprocession
            
                matlabbatch{1}.spm.stats.fmri_spec.sess(nr).multi_reg = [];
                matlabbatch{1}.spm.stats.fmri_spec.sess(nr).multi_reg = cellstr(spm_select('FPList',pwd,'^corr'));
            end
            
            
            if found_epi == 0 % incase no scans were found: throw away
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {''};
            end
            
            cd(root_dir);  % add path
            cd(cur_sub);
            cd('MRT/NIFTI/PDT')        
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
    
    %keyboard
    
    % fill in microtime resolution
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = cur_nslices;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = cur_refslice;
    % specifiy an explicit mask, which voxels shall only be analysed
  %  matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.mask = cellstr(expl_mask);
    % specify TR-time
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.0; %% CHANGE!
    
%     %%%%%%%%%%%TO DO: CHANGE THIS adjusting to new contrasts    % contrast manager
    con.contrastNames = {'pic.on.gam','pic.on.neg','pic.on.pos','picgam.on.gain','picgam.on.loss','picgam.on.ed', 'picgam.opt.on', 'acceptance'};  %changed FC
    con.contrastType = cellstr(repmat('t',length(con.contrastNames),1))';
    con.contrastWeights = agk_basic_t_cons(length(con.contrastType));
    %con.contrastWeights = agk_basic_t_cons(length(con.contrastType));
    % repeat and reslice
    con.contrastRep = cellstr(repmat('replsc',length(con.contrastNames),1))'; 
    
    
    for jj = 1:length(con.contrastType)
        matlabbatch{3}.spm.stats.con.consess{jj}.tcon.name    = con.contrastNames{jj};

        matlabbatch{3}.spm.stats.con.consess{jj}.tcon.sessrep = con.contrastRep{jj};
    end
    
    for jj = 1:length(con.contrastWeights)
        matlabbatch{3}.spm.stats.con.consess{jj}.tcon.convec  = con.contrastWeights{jj};
    end
    
    % dependencies: estimation to model specification
    %               contrast manager to model estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

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

% catch MExc
%     disp('Something went wrong');
%     error_message = MExc;
% end
