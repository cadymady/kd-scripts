% Charité Berlin, 12.09.2017 - VPPG Study, PDT data

% compare P structure for every subject and data_pdt matrix (created in R)
clear all

cd('T:\Library\R')
T = readtable('data_pdt_behav.csv', 'Delimiter', ','); % 
T = sortrows(T, {'subject', 'trial'}, {'ascend'});          % sort subject and trial in ascending order
C = table2cell(T);                                     % convert to cell to access

% columns of interest
subject = 2;
gain = 5;
loss = 6;
choice = 7;
rt = 8;
cat = 9;
st_dur = 11;
ed_abs = 18;

base_dir = 'E:\Master\Praktikum Charité\VPPG\data\';        % base directory
cd(base_dir)    
all_subj = cellstr(ls('VPPG*'));                            % list all subject folders
a = 0;              % counter for cell array with subjects that are not included yet (see below)

for ii = 1:size(all_subj,1)
    cd(base_dir)
    cd(all_subj{ii})            % go into each subject's folder
    cd('Behav\PDT')
    
    
    load(ls('P_*'))             % load file with BEHAV data
    COPYP = P;                  % copy to not mess with structure
    save('COPYP.mat', 'COPYP')

    C_subj = strcmp(C, all_subj{ii});       % look for current subject in C
    
    
    s = 1;          % find starting point of subject in data_pdt cell
    while C_subj(s,2) == 0 && s < size(C,1)
        s = s+1;
    end
    
    st = s;        % find end point
    while C_subj(st,2) == 1 && st < size(C,1)
        st = st+1;
    end
    
   % try
        if s < size(C,1)            % for those subjects that have been found, compare ed (and choice)
            disp('exists')
        C_choice = transpose(C(s:st,choice));
        C_ed_abs = transpose(C(s:st,ed_abs));
        C_cat = transpose(C(s:st,cat));
        c = isequal(P.cur.choice,cell2mat(C_choice)); % equal
        e = isequal(P.gmat.ed_abs,cell2mat(C_ed_abs));
       
        if e == 0 % if eds differ, change in MATLAB structure
            P.gmat.ed_abs = cell2mat(C_ed_abs);
        end
        save(ls('P_*'), 'P')      % save new MATLAB structure in subject's folder
        
    %catch
        elseif s == size(C,1)           % if subject hasn't been found, continue
            disp(['Ed_abs for Subject ', all_subj{ii}, ' missing . I will continue to the next one.'])
            a = a +1;
            no_ed_abs{a} = all_subj{ii};    
            
        end
    %end   
   
    
end
cd('E:\Master\Praktikum Charité\VPPG\data')
save('subj_missing.mat', 'no_ed_abs')
