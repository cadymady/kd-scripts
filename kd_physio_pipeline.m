% Physio correction
clear all

base_dir_lib  = 'T:\Library\MATLAB\PDT';
base_dir_lib2 = 'T:\Library\MATLAB\Physio\scripts';
% where are the subject folders with MRI data?
base_dir_pl = 'E:\Master\Praktikum Charité\VPPG\data\';
% Set spm Mask Threshold

 cd(base_dir_pl);
 all_subf = cellstr(ls('VPPG*'));  %list folders starting w vppg
 session = 'PDT';
 no_physio = {};
 pp = 0;
 ow = 1;
 
 for ii = 62
     cd(base_dir_pl)
     cd(all_subf{ii})
     if exist(fullfile(pwd,'Physio'), 'file')
        [error_msg{ii}, error_sub{ii}] = kd_physio_batch( base_dir_pl,all_subf{ii}, session, ow );
         
     else
         % if there is no physio file
         %disp(['No Physio folder found for subject ', all_subf{ii}, ' I will continue to the next one.']);
         pp = pp +1;
         no_physio{pp} = all_subf{ii};
         continue
     end
 end
 
 save('Error_report.mat', 'error_msg');