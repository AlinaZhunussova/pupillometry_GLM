% Analysis developed in collaboration with my supervisor, Dorothea Hämmerer
% First-level GLM (eye_glm.m) by Tobias Hauser (2015).
% Second-level OLS (ols.m) by Tim Behrens (2004).

%clear
%close all

projroot = pwd;   % or e.g. '/path/to/pupil_project_2'

addpath(fullfile(projroot, 'raacampbell-shadedErrorBar-19cf3fe'));
addpath(fullfile(projroot, 'fieldtrip'));

pathdat   = [fullfile(projroot, 'EM_task', 'EM_cleaneddata_test') filesep];
pathsave  = [fullfile(projroot, 'EM_task', 'GLM_outputs') filesep];
pathbehav = [fullfile(projroot, 'EM_task', 'behav_em') filesep];
pathglm   = [fullfile(projroot, 'EM_task', 'files_for_GLM') filesep];

missingst1 = [106 107 111 115 117 202 212 216 222 302 311 319 332];
missingst2 = [105 110 202 211 214 302 315 317 319 329 330 332];

rangeids = [105 106 107 108 109 110 111 113 114 115 116 117 118 119 121 122 ...
   123 124 125 127 128 129 130 131 132 133 134 135 202 203 204 205 206 207 ...
   208 209 211 212 213 214 215 216 217 220 221 222 223 225 226 227 228 ...
   229 230 301 302 304 305 307 308 309 310 311 312 313 314 315 316 317 318 319 320 ...
   321 323 324 325 326 327 328 329 330 331 332 333];

    has_st1 = ~ismember(rangeids, missingst1);
    has_st2 = ~ismember(rangeids, missingst2);
    both_sessions = has_st1 & has_st2;
    n_both = sum(both_sessions);
    ids_with_both_sessions = rangeids(both_sessions)';


    yas = find(ids_with_both_sessions < 200);
    oas = find(ids_with_both_sessions > 200 & ids_with_both_sessions < 300);
    mci = find(ids_with_both_sessions > 300);


% permutation settings
do_randpermute = 1;
length_randpermute_indiv = 100; % how many shuffles at individual level
length_randpermute = 1000; % how many shuffles at the group level
tw = 701:2501; % 701:2501 0.5 to 2.3s - time-window for multiple comparison on timeseries 

clear cib vib tib store_R storedat storebslnoise

for i = 1:length(ids_with_both_sessions)
    test_id = ids_with_both_sessions(i);
    clear beh  

    % test session 1
    load([pathdat 'r' num2str(test_id) '_eye_rv200_1.mat'], 'eye_rv');
    eye_rv1 = eye_rv;

    % test session 2
    load([pathdat 'r' num2str(test_id) '_eye_rv200_2.mat'], 'eye_rv');
    eye_rv2 = eye_rv;

    % encoding session (scanner data)
    load([pathbehav num2str(test_id) '/behav/EM/' num2str(test_id) '_emomem.mat']); % contains variable called exp
    task = exp; clear exp

    if test_id == 325
        eye_rv1.fsample = 1000;   % weirdly labeled as 1000.001
        eye_rv2.fsample = 1000;
    end

    % resample each session to 1000 Hz if needed
    for s = 1:2
        if s == 1, tmp = eye_rv1; else, tmp = eye_rv2; end
        if tmp.fsample ~= 1000
            cfg = [];
            cfg.resamplefs = 1000;
            tmp = ft_resampledata(cfg, tmp);
            for tr = 1:length(tmp.trial)
                tmp.trial{tr} = tmp.trial{tr}(1:2701);
                tmp.time{tr}  = tmp.time{tr}(1:2701);
            end
        end
        if s == 1, eye_rv1 = tmp; else, eye_rv2 = tmp; end
    end

    if test_id == 325
        eye_rv.fsample = 1000;
        eye_rv2.trlraus(60:88) = 1;   
    end
        
    trlnum_s1 = find(eye_rv1.trlraus == 0);
    trlnum_s2 = find(eye_rv2.trlraus == 0);

    zero_exclude = struct();
    zero_exclude(1).id    = 119;
    zero_exclude(1).trials = 14;
    zero_exclude(2).id    = 122;
    zero_exclude(2).trials = 124;

    eye_rv = eye_rv1;
    eye_rv.trial = [eye_rv1.trial       eye_rv2.trial];
    eye_rv.time = [eye_rv1.time        eye_rv2.time];
    eye_rv.stimname       = [eye_rv1.stimname(trlnum_s1)       eye_rv2.stimname(trlnum_s2)];
    eye_rv.stimtypeneuemo = [eye_rv1.stimtypeneuemo(trlnum_s1) eye_rv2.stimtypeneuemo(trlnum_s2)];
    eye_rv.respcert       = [eye_rv1.respcert(trlnum_s1)       eye_rv2.respcert(trlnum_s2)];
    
    fprintf('trial=%d time=%d stimname=%d stimtype=%d respcert=%d\n', ...
    length(eye_rv.trial), length(eye_rv.time), ...
    length(eye_rv.stimname), length(eye_rv.stimtypeneuemo), length(eye_rv.respcert));

    good_trials = true(1, length(eye_rv.trial)); 

    excl_idx = find([zero_exclude.id] == test_id);
    if ~isempty(excl_idx)
        bad_trials = zero_exclude(excl_idx).trials;
        good_trials(bad_trials) = false;
        
        eye_rv.trial          = eye_rv.trial(good_trials);
        eye_rv.time           = eye_rv.time(good_trials);
        eye_rv.stimname       = eye_rv.stimname(good_trials);
        eye_rv.stimtypeneuemo = eye_rv.stimtypeneuemo(good_trials);
        eye_rv.respcert       = eye_rv.respcert(good_trials);

        fprintf('Subject %d: excluded trial(s) %s (zero trial)\n', test_id, num2str(bad_trials));
    end

    length(eye_rv1.trial) 
    length(eye_rv2.trial)
    length(eye_rv.trial)

    % determine which stimuli are old
    eye_rv.oldnew = zeros(1,length(eye_rv.stimname));
    
    for ste1 = 1:length(eye_rv.stimname)
        for sta = 1:length(task.dat.stimname)
            if strmatch(eye_rv.stimname{ste1},task.dat.stimname{sta})
                eye_rv.oldnew(ste1) = 1;
            end
        end
    end

    % determine certain/uncertain responses
    if ~ismember(28,eye_rv.respcert)
        rGUESS = 26; rSURE = 13;
    elseif ismember(28,eye_rv.respcert)
        rGUESS = 28; rSURE = 31;
    end
    eye_rv.respcert(find(eye_rv.respcert == rGUESS)) = 1;
    eye_rv.respcert(find(eye_rv.respcert == rSURE)) = 2;
        
    %n_trials_s1 = length(trlnum_s1);
    %n_trials_s2 = length(trlnum_s2);
    %n_trials    = n_trials_s1 + n_trials_s2;
    
    n_trials_s1 = length(trlnum_s1);
    n_trials_s2 = length(trlnum_s2);
    retention_raw_full = [zeros(1, n_trials_s1), ones(1, n_trials_s2)]; % create a vector for indicating whether it is test 1 or 2
    
    
    emotionality_all = eye_rv.stimtypeneuemo; 
    recognition_all  = eye_rv.oldnew;         
    certainty_all    = eye_rv.respcert;       
    retention_raw = retention_raw_full(good_trials);  
    n_trials = sum(good_trials); 

    % z-score and baseline correct
    dat = cell2mat(eye_rv.trial'); % trials x timepoints
    meandat = nanmean(dat, 2); % each value is an average pupil size across the entire trial duration # mean(dat(1:2701,1))
    stddat = nanstd(dat'); % timepoints x trials - each value is std of pupil size across time for that trial
    dat = dat-repmat(meandat, 1, size(dat, 2)); % subtract mean from raw pupil data for each trial
    dat = dat ./ repmat(stddat', 1, size(dat, 2)); % zscore = (X-MEAN(X)) ./ STD(X)

    storebslnoise(i) = nanmean(nanstd(dat(:,1:200)'));% mean of std within 200ms baseline per person, assess whether baseline noise is different
        
    % baseline correct (0 to 200 ms time window)
     for trl = 1:size(dat, 1)
         dat(trl,:) = dat(trl, :)- nanmean(dat(trl,1:200)); % count the mean of first 200 data points and substract from each timepoint
     end
 
    % recode as binary
    beh.emotionality = double(emotionality_all == 2); % 1 - emotional, 0 - neutral
    beh.recognition = double(recognition_all == 1); % old (shown in th scanner) - 1, new - 0           
    beh.certainty = double(certainty_all == 2); % 1 - sure, 0 - guess
    beh.retention = retention_raw; % test 1 - 0, test 2 - 1

    % z-score everything
    beh.emotionality = zscore(beh.emotionality);
    beh.recognition = zscore(beh.recognition);
    beh.certainty = zscore(beh.certainty);
    beh.retention = zscore(beh.retention);
    % interactions
    beh.recognition_x_emotionality = beh.recognition .* beh.emotionality;
    beh.recognition_x_certainty = beh.recognition .* beh.certainty;
    
    beh.recognition_x_emotionality = zscore(beh.recognition .* beh.emotionality);
    beh.recognition_x_certainty    = zscore(beh.recognition .* beh.certainty);
    
    % catch zero-variance regressors that zscore turns into NaN 
    raw_regs = struct('emotionality', double(emotionality_all==2), ...
                      'recognition',  double(recognition_all==1), ...
                      'certainty',    double(certainty_all==2), ...
                      'retention',    retention_raw);
  
    fn = fieldnames(raw_regs);
    for r = 1:numel(fn)
        if nanstd(raw_regs.(fn{r})) == 0
            fprintf('  ZERO-VARIANCE: subject %d, regressor "%s" (n=%d) -> zscore = NaN\n', ...
                    test_id, fn{r}, numel(raw_regs.(fn{r})));
        end
    end

    % check collinearity
    reg_matrix = [beh.emotionality', beh.recognition', beh.certainty', beh.retention', beh.recognition_x_emotionality', beh.recognition_x_certainty'];

    store_R(i, :, :) = corrcoef(reg_matrix);

    pupil_custom = struct();
    pupil_custom.trial = cell(n_trials, 1);
    pupil_custom.time  = repmat({eye_rv.time{1}}, n_trials, 1);
   % pupil_custom.sampleinfo = eye_rv.sampleinfo; % not used in eye_glm.m
   % for this analysis
   % pupil_custom.fsample = eye_rv.fsample; % not used in eye_glm.m
   % pupil_custom.label = eye_rv.label; % not used in eye_glm.m

    for trl = 1:n_trials
        pupil_custom.trial{trl} = dat(trl, :);
    end

    reg_names = {'emotionality', 'recognition', 'certainty', 'retention', ...
                  'recognition_x_emotionality', 'recognition_x_certainty'};

    contrast = [1	0	0	0 	0   0   0   % emotionality
                0	1	0	0	0   0   0   % recognition
                0	0	1	0	0   0   0   % certainty
                0	0	0	1	0	0   0   % testing period
                0   0   0   0   1   0   0   % recognition_x_emotionality
                0   0   0   0   0   1   0   % recognition_x_certainty
                0	0	0	0	0   0   1]; % intercept

    con_name = [reg_names, {'intercept'}];
    
    options.rejBad = 0;
    options.bl_corr = 0;
    options.verbose = 1;

    options.detrend = 0;
    glm_result = eye_glm(pupil_custom, beh, reg_names, contrast, con_name, options);

    storedat(i).glm = glm_result;
    storedat(i).id  = test_id;


    if do_randpermute
        for ri = 1:length_randpermute_indiv
            ind = randperm(n_trials); % makes a random reordering of the trial numbers
            
            % shuffle labels in all regressors  
            beh_perm = beh;
            for bl = 1:length(reg_names)
                beh_perm.(reg_names{bl}) = beh.(reg_names{bl})(ind); % the same ind is applied to every regressor
            end
            
            dum = eye_glm(pupil_custom, beh_perm, reg_names, contrast, con_name, options);
            
            for st = 1:length(reg_names)+1
                cib(i, st, ri, :) = dum.con{1, st}; % i - subject, st - contrast, ri - permutation number, timepoints
                vib(i,st,ri,:) = dum.res{1,st};
                tib(i,st,ri,:) = dum.tstat{1,st};
            end
            clear ind dum beh_perm
        end
    end
end    

% correlation matrix plot
reg_labels = {'emotionality', 'recognition', 'certainty', 'retention', 'recognition_x_emotionality', 'recognition_x_certainty'};
mean_R = squeeze(nanmean(store_R, 1));
std_R  = squeeze(nanstd(store_R, 1));
imagesc(mean_R);
colorbar;
colormap(jet);
clim([-1 1]);
set(gca, 'XTick', 1:length(reg_labels), 'XTickLabel', reg_labels, ...
         'YTick', 1:length(reg_labels), 'YTickLabel', reg_labels);
title('Regressor correlation matrix');

% within group colors
col_ya  = [0.106 0.624 0.467];  % #1B9E77 teal-green
col_oa  = [0.216 0.494 0.722];
col_mci = [0.835 0.243 0.310];  % muted crimson red

% significance colors
col_between_10 = [0.6 0.6 0.6];  % grey for p<.10
col_between_05 = [0.2 0.2 0.2];  % dark grey/black for p<.05

col_ya_05  = [0.050 0.350 0.250];  % darker teal-green
col_oa_05  = [0.100 0.280 0.450];  % darker cornflower blue
col_mci_05 = [0.500 0.100 0.150]; 

% significance bar y-positions per contrast 
sign_ypos(1).name    = 'emotionality';
sign_ypos(1).within  = [-0.11, -0.12, -0.10];
sign_ypos(1).between = [-0.14, -0.15, -0.16];

sign_ypos(2).name    = 'recognition';
sign_ypos(2).within  = [-0.10, -0.12, -0.14];
sign_ypos(2).between = [-0.28, -0.18, -0.23];

sign_ypos(3).name    = 'certainty';
sign_ypos(3).within  = [-0.06, -0.15, -0.17];
sign_ypos(3).between = [-0.20, -0.24, -0.28];

sign_ypos(4).name    = 'retention';
sign_ypos(4).within  = [-0.31, -0.33, -0.35];
sign_ypos(4).between = [-0.38, -0.40, -0.42];

sign_ypos(5).name    = 'recognition_x_emotionality';
sign_ypos(5).within  = [-0.06, -0.08, -0.04];  
sign_ypos(5).between = [-0.09, -0.11, -0.13];

sign_ypos(6).name    = 'recognition_x_certainty';
sign_ypos(6).within  = [-0.06, -0.07, -0.08];
sign_ypos(6).between = [-0.09, -0.10, -0.13];


sign_ylim = [-0.22  0.20;   % emotionality
             -0.45  0.30;   % recognition
             -0.35  0.25;   % certainty
             -0.45  0.20;   % retention
             -0.15  0.20;   % recognition_x_emotionality 
             -0.15  0.20];  % recognition_x_certainty

% plot across people
time = eye_rv.time{1};  
range = 1:length(time);
n_subs    = length(ids_with_both_sessions);
n_reg     = length(reg_names);
n_time    = length(time);

group_dm = zeros(n_subs, 3);
group_dm(yas, 1) = 1;
group_dm(oas, 2) = 1;
group_dm(mci, 3) = 1;

load(fullfile(pathglm, 'RTs_Emo_minus_Neu_T1_for_GLM.mat'));
load(fullfile(pathglm, 'RTs_Emo_minus_Neu_T2_for_GLM.mat'));

behav_group = nan(length(ids_with_both_sessions), 1);
for i = 1:length(ids_with_both_sessions)
    idx_T1 = find(rt_T1.id == ids_with_both_sessions(i));
    idx_T2 = find(rt_T2.id == ids_with_both_sessions(i));
    val_T1 = NaN; val_T2 = NaN;
    if ~isempty(idx_T1), val_T1 = rt_T1.behav_group(idx_T1); end
    if ~isempty(idx_T2), val_T2 = rt_T2.behav_group(idx_T2); end
    behav_group(i) = nanmean([val_T1, val_T2]); % if a subject has RT for only one session, they still get a valid value
end

disp(['behav_group NaNs: ' num2str(sum(isnan(behav_group)))]);
disp(['behav_group complete: ' num2str(sum(~isnan(behav_group)))]);

load(fullfile(pathglm, 'dprime_Emo_minus_Neu_T1_for_GLM.mat'));
load(fullfile(pathglm, 'dprime_Emo_minus_Neu_T2_for_GLM.mat'));

behav_rec = nan(length(ids_with_both_sessions), 1);
for i = 1:length(ids_with_both_sessions)
    idx_T1 = find(dp_T1.id == ids_with_both_sessions(i));
    idx_T2 = find(dp_T2.id == ids_with_both_sessions(i));
    val_T1 = NaN; val_T2 = NaN;
    if ~isempty(idx_T1), val_T1 = dp_T1.behav_group(idx_T1); end
    if ~isempty(idx_T2), val_T2 = dp_T2.behav_group(idx_T2); end
    behav_rec(i) = nanmean([val_T1, val_T2]);
end

disp(['behav_rec NaNs: ' num2str(sum(isnan(behav_group)))]);
disp(['behav_rec complete: ' num2str(sum(~isnan(behav_group)))]);

group_con = [1 -1  0;   % YA > OA
             1  0 -1;   % YA > MCI
             0  1 -1];  % OA > MCI

n_group_con = size(group_con, 1);  % = 3

clear cgb vgb tgb
cgb = nan(n_reg, length_randpermute, 4, n_time); % 4 - 1. YA > OA, 2. YA > MCI, 3. OA > MCI, 4. the behavioural covariate (RT for contrIND 1, d′ for contrIND 2)
vgb = nan(n_reg, length_randpermute, 4, n_time);
tgb = nan(n_reg, length_randpermute, 4, n_time);

store_group_glm = cell(1, n_reg);

display_names = {'Emotionality', ...
                 'Recognition', ...
                 'Certainty', ...
                 'Testing period', ...
                 'Recognition x Emotionality', ...   
                 'Recognition x Certainty', ...
                 'Intercept'};

fprintf('\n=== NaN diagnostics - behavioural covariates ===\n\n');
fprintf('behav_group NaNs: %d subjects\n', sum(isnan(behav_group)));
nan_idx = find(isnan(behav_group));
for k = 1:length(nan_idx)
    fprintf('  Subject %d (index %d)\n', ids_with_both_sessions(nan_idx(k)), nan_idx(k));
end

fprintf('behav_rec NaNs: %d subjects\n', sum(isnan(behav_rec)));
nan_idx = find(isnan(behav_rec));
for k = 1:length(nan_idx)
    fprintf('  Subject %d (index %d)\n', ids_with_both_sessions(nan_idx(k)), nan_idx(k));
end

% text label y-positions 
% format: [YA>OA, YA>MCI, OA>MCI]
text_ypos = containers.Map(...
    {'emotionality', 'recognition', 'certainty', 'retention', ...
     'recognition_x_emotionality', 'recognition_x_certainty'}, ...
    {[-0.14, -0.15, -0.16], ...   % emotionality
     [-0.28, -0.18, -0.23], ...   % recognition
     [-0.20, -0.24, -0.28], ...   % certainty
     [-0.38, -0.40, -0.42], ...   % retention
     [-0.09, -0.11, -0.13], ...   % recognition_x_emotionality
     [-0.09, -0.10, -0.13]});     % recognition_x_certainty

fprintf('\n=== peak latencies ===\n\n');
contrast_names = {'emotionality', 'recognition', 'certainty', 'retention', ...
                  'recognition_x_emotionality', 'recognition_x_certainty'};

group_names = {'YA', 'OA', 'MCI'};
group_idx_list = {yas, oas, mci};

for contrIND = 1:length(con_name)

    dat1 = zeros(length(ids_with_both_sessions), length(time));
    for i = 1:length(ids_with_both_sessions)
        dat1(i,:) = storedat(i).glm.con{1, contrIND};
    end
    c = dat1;

    mean_yas = nanmean(dat1(yas, :)); % mean across all YAs subjects at that timepoint
    se_yas   = nanstd(dat1(yas, :)) ./ sqrt(length(yas));
    mean_oas = nanmean(dat1(oas, :));
    se_oas   = nanstd(dat1(oas, :)) ./ sqrt(length(oas));
    mean_mci = nanmean(dat1(mci, :));
    se_mci   = nanstd(dat1(mci, :)) ./ sqrt(length(mci));

    if contrIND <= n_reg
       fprintf('--- %s ---\n', contrast_names{contrIND});
    
        for g = 1:3
            grp = group_idx_list{g};
            mean_grp = nanmean(dat1(grp, :));
            
            % find both the peak (max) and trough (min) within the window
            [peak_val,   peak_idx]   = max(mean_grp(tw));
            [trough_val, trough_idx] = min(mean_grp(tw));
            
            if abs(peak_val) >= abs(trough_val)
                chosen_idx = peak_idx;
                direction  = 'positive';
            else
                chosen_idx = trough_idx;
                direction  = 'negative';
            end
            
            peak_time = time(tw(chosen_idx));
            fprintf('%s peak: %.2f s (%s)\n', group_names{g}, peak_time, direction);
        end
        fprintf('\n');
    end

   if contrIND <= n_reg
        bad_subs  = any(isnan(dat1), 2);
        good_subs = ~bad_subs;
            
        c_clean        = dat1(good_subs, :);
        group_dm_clean = group_dm(good_subs, :);

        bg = behav_group(good_subs);
        bg_mean = nanmean(bg);
        bg_std  = nanstd(bg);
        behav_group_z = (bg - bg_mean) ./ bg_std;

        % set YAs and OAs to 0 (mean after zscoring) - only MCI contributes to RT effect
        behav_group_z_mci = behav_group_z;
        
        % find which indices in good_subs correspond to YAs and OAs
        good_subs_idx = find(good_subs);
        yas_in_good = ismember(good_subs_idx, yas);
        oas_in_good = ismember(good_subs_idx, oas);

        if contrIND == 1
        % z-score RT within MCI only first
        mci_in_good = ismember(good_subs_idx, mci);
        rt_mci_vals = behav_group_z(mci_in_good);
        behav_group_z_mci(mci_in_good) = (rt_mci_vals - nanmean(rt_mci_vals)) ./ nanstd(rt_mci_vals);
    
        % then zero out YAs and OAs
        behav_group_z_mci(yas_in_good) = 0;
        behav_group_z_mci(oas_in_good) = 0;
            
        % additionally exclude subjects with NaN in behav_group_z_mci
        nan_covar = isnan(behav_group_z_mci);
        if any(nan_covar)
            fprintf('contrIND 1: excluding %d subject(s) with NaN RT covariate\n', sum(nan_covar));
            c_clean        = c_clean(~nan_covar, :);
            group_dm_behav = [group_dm_clean(~nan_covar, :), behav_group_z_mci(~nan_covar)];
        else
            group_dm_behav = [group_dm_clean, behav_group_z_mci];
        end

        elseif contrIND == 2
            br = behav_rec(good_subs);
            behav_rec_z = nan(size(br));
            
            good_subs_idx = find(good_subs);
            yas_in_good = ismember(good_subs_idx, yas);
            oas_in_good = ismember(good_subs_idx, oas);
            mci_in_good = ismember(good_subs_idx, mci);
    
            % z-score within each group separately
            for grp_mask = {yas_in_good, oas_in_good, mci_in_good}
                m = grp_mask{1};
                vals = br(m);
                grp_mean = nanmean(vals);
                grp_std  = nanstd(vals);
                if grp_std > 0
                    behav_rec_z(m) = (vals - grp_mean) ./ grp_std;
                else
                    behav_rec_z(m) = 0;  % edge case: no variance in group
                    fprintf('Warning: zero std in dprime for a group - setting covariate to 0\n');
                end
            end
        
            fprintf('Within-group z-scored dprime stats:\n');
            fprintf('  YAs:  mean=%.4f, std=%.4f\n', nanmean(behav_rec_z(yas_in_good)), nanstd(behav_rec_z(yas_in_good)));
            fprintf('  OAs:  mean=%.4f, std=%.4f\n', nanmean(behav_rec_z(oas_in_good)), nanstd(behav_rec_z(oas_in_good)));
            fprintf('  MCI:  mean=%.4f, std=%.4f\n', nanmean(behav_rec_z(mci_in_good)), nanstd(behav_rec_z(mci_in_good)));
    
            % additionally exclude subjects with NaN in behav_rec_z
            nan_covar = isnan(behav_rec_z);
            if any(nan_covar)
                fprintf('contrIND 2: excluding %d subject(s) with NaN dprime covariate\n', sum(nan_covar));
                c_clean        = c_clean(~nan_covar, :);
                group_dm_clean = group_dm_clean(~nan_covar, :);
                behav_rec_z    = behav_rec_z(~nan_covar);
                group_dm_behav = [group_dm_clean, behav_rec_z];
            else
                group_dm_behav = [group_dm_clean, behav_rec_z];
            end
        else 
            group_dm_behav = [group_dm_clean, zeros(size(group_dm_clean, 1), 1)];
        end  


        group_con_behav = [group_con, zeros(3,1);  
                            0  0  0  1];  

        [cg, vg, tg] = ols(c_clean, group_dm_behav, group_con_behav);

        store_group_glm{contrIND}.cg = cg; % contrast estimates (differences) of how much bigger is the #YA effect than the #OA effect at each timepoint
        store_group_glm{contrIND}.vg = vg;
        store_group_glm{contrIND}.tg = tg;

        disp(['starting permutations for contrIND=' num2str(contrIND)]) 

        % update n_good after covariate NaN exclusion
        if contrIND == 1 || contrIND == 2
            n_good = size(c_clean, 1);
        else
            n_good = sum(good_subs);
        end

        for r = 1:length_randpermute
            ind    = randperm(n_good);
            c_perm = c_clean(ind, :);
            [tmp_cg, tmp_vg, tmp_tg] = ols(c_perm, group_dm_behav, group_con_behav);
            cgb(contrIND, r, :, :) = reshape(tmp_cg, [1 1 4 n_time]);
            vgb(contrIND, r, :, :) = reshape(tmp_vg, [1 1 4 n_time]);
            tgb(contrIND, r, :, :) = reshape(tmp_tg, [1 1 4 n_time]);
        end


        store_group_glm{contrIND}.cgb = squeeze(cgb(contrIND,:,:,:)); % cgb - null (fake) contrast estimates — the difference (YA>OA, YA>MCI, OA>MCI, RT) for every one of the 1000 shuffles
        store_group_glm{contrIND}.vgb = squeeze(vgb(contrIND,:,:,:)); % cg - real contrast estimates
        store_group_glm{contrIND}.tgb = squeeze(tgb(contrIND,:,:,:));

   end

    plotsign_yas_10 = nan(1, n_time);
    plotsign_yas_05 = nan(1, n_time);
    plotsign_oas_10 = nan(1, n_time);
    plotsign_oas_05 = nan(1, n_time);
    plotsign_mci_10 = nan(1, n_time);
    plotsign_mci_05 = nan(1, n_time);
    plotsign_10     = nan(n_group_con, n_time);
    plotsign_05     = nan(n_group_con, n_time);

    if do_randpermute && contrIND <= n_reg
        yp = sign_ypos(contrIND);

        % YAs
        datdum = squeeze(nanmean(cib(yas, contrIND, :, :), 1)); % cib = fake contrast estimates at the within-subject level, 100 fake group-level contrast estimates for YAs
      
        for rep = 1:length_randpermute_indiv
            datmin(rep) = min(datdum(rep, min(tw):max(tw))); % a 100-long list of minimums
            datmax(rep) = max(datdum(rep, min(tw):max(tw))); % a 100-long list of maximums
        end

        plotsign_yas_10(find(mean_yas(min(tw):max(tw)) < prctile(datmin, 10) | ...
                             mean_yas(min(tw):max(tw)) > prctile(datmax, 90)) + min(tw)) = yp.within(1); % Here mean_yas (the real one) is compared at every timepoint against the one threshold number from min(datdum) or max(datdum)
        plotsign_yas_05(find(mean_yas(min(tw):max(tw)) < prctile(datmin,  5) | ...
                             mean_yas(min(tw):max(tw)) > prctile(datmax, 95)) + min(tw)) = yp.within(1); % you compare each point literally to that one prctile(dat...) number, and only points above it get marked
                                                                                                         % 90% of the fake peaks fall below. So if the real curve is above it, the real effect is bigger than 90% of what chance produced, in the top 10% of the null
                                                                                                         % So if the real curve drops below this threshold, the real dip is deeper than 90% of chance dips, in the bottom 10% of the null

        clear datmin datmax rep datdum

        % OAs
        datdum = squeeze(nanmean(cib(oas, contrIND, :, :), 1));
        for rep = 1:length_randpermute_indiv
            datmin(rep) = min(datdum(rep, min(tw):max(tw)));
            datmax(rep) = max(datdum(rep, min(tw):max(tw)));
        end
        plotsign_oas_10(find(mean_oas(min(tw):max(tw)) < prctile(datmin, 10) | ...
                             mean_oas(min(tw):max(tw)) > prctile(datmax, 90)) + min(tw)) = yp.within(2);
        plotsign_oas_05(find(mean_oas(min(tw):max(tw)) < prctile(datmin,  5) | ...
                             mean_oas(min(tw):max(tw)) > prctile(datmax, 95)) + min(tw)) = yp.within(2);

        clear datmin datmax rep datdum

        % MCI
        datdum = squeeze(nanmean(cib(mci, contrIND, :, :), 1));
        for rep = 1:length_randpermute_indiv
            datmin(rep) = min(datdum(rep, min(tw):max(tw)));
            datmax(rep) = max(datdum(rep, min(tw):max(tw)));
        end
        plotsign_mci_10(find(mean_mci(min(tw):max(tw)) < prctile(datmin, 10) | ...
                             mean_mci(min(tw):max(tw)) > prctile(datmax, 90)) + min(tw)) = yp.within(3);
        plotsign_mci_05(find(mean_mci(min(tw):max(tw)) < prctile(datmin,  5) | ...
                             mean_mci(min(tw):max(tw)) > prctile(datmax, 95)) + min(tw)) = yp.within(3);

        clear datmin datmax rep datdum

        store_sigwindow(contrIND).yas_05 = plotsign_yas_05;
        store_sigwindow(contrIND).oas_05 = plotsign_oas_05;
        store_sigwindow(contrIND).mci_05 = plotsign_mci_05;
        store_sigwindow(contrIND).yas_10 = plotsign_yas_10;
        store_sigwindow(contrIND).oas_10 = plotsign_oas_10;
        store_sigwindow(contrIND).mci_10 = plotsign_mci_10;

        if contrIND == 1  % only for emotionality
            for rep = 1:size(cgb, 2)
                datmin(rep) = min(squeeze(cgb(contrIND, rep, 4, min(tw):max(tw))));
                datmax(rep) = max(squeeze(cgb(contrIND, rep, 4, min(tw):max(tw))));
            end

        plotsign_rt_10 = nan(1, n_time);
        plotsign_rt_05 = nan(1, n_time);
        plotsign_rt_10(find(cg(4, min(tw):max(tw)) < prctile(datmin, 10) | ...
                            cg(4, min(tw):max(tw)) > prctile(datmax, 90)) + min(tw)) = yp.between(3) - 0.02;
        plotsign_rt_05(find(cg(4, min(tw):max(tw)) < prctile(datmin,  5) | ...
                        cg(4, min(tw):max(tw)) > prctile(datmax, 95)) + min(tw)) = yp.between(3) - 0.02;
        clear datmin datmax
        end

        if contrIND == 2  % dprime significance for recognition
            for rep = 1:size(cgb, 2)
                datmin(rep) = min(squeeze(cgb(contrIND, rep, 4, min(tw):max(tw))));
                datmax(rep) = max(squeeze(cgb(contrIND, rep, 4, min(tw):max(tw))));
            end
            plotsign_dp_10 = nan(1, n_time);
            plotsign_dp_05 = nan(1, n_time);
            plotsign_dp_10(find(cg(4, min(tw):max(tw)) < prctile(datmin, 10) | ...
                                cg(4, min(tw):max(tw)) > prctile(datmax, 90)) + min(tw)) = yp.between(1) - 0.1;
            plotsign_dp_05(find(cg(4, min(tw):max(tw)) < prctile(datmin,  5) | ...
                                cg(4, min(tw):max(tw)) > prctile(datmax, 95)) + min(tw)) = yp.between(1) -  0.1;
            clear datmin datmax
        end
    end

    if do_randpermute && contrIND <= n_reg
        for grcon = 1:size(group_con,1) % group contrasts
            for rep = 1:length_randpermute
                datmin(rep) = min(squeeze(cgb(contrIND,rep,grcon,min(tw):max(tw))));
                datmax(rep) = max(squeeze(cgb(contrIND,rep,grcon,min(tw):max(tw))));
            end
    
            plotsign_10(grcon, find(cg(grcon,min(tw):max(tw)) < prctile(datmin, 10)) + min(tw)) = yp.between(grcon); % we take cg which are the real contrast estimates for a specific group comparison (e.g. YAs vs OAs) and compare each timepoint to a single min or max threshold from prctile(datmin/datmax) built from the 1000 shuffles
            plotsign_10(grcon, find(cg(grcon,min(tw):max(tw)) > prctile(datmax, 90)) + min(tw)) = yp.between(grcon);
            plotsign_05(grcon, find(cg(grcon,min(tw):max(tw)) < prctile(datmin,  5)) + min(tw)) = yp.between(grcon);

            plotsign_05(grcon, find(cg(grcon,min(tw):max(tw)) > prctile(datmax, 95)) + min(tw)) = yp.between(grcon);
            clear datmin datmax
        end
    end   

    figure('Name', con_name{contrIND});
    set(gcf, 'Color', 'white');
    set(gcf, 'Position', [100 100 1800 1200]); 

    shadedErrorBar(time(range), mean_yas(range), se_yas(range), ...
        'lineProps', {'Color', col_ya, 'LineWidth',2}); hold on
    shadedErrorBar(time(range), mean_oas(range), se_oas(range), ...
        'lineProps', {'Color', col_oa, 'LineWidth',2});
    shadedErrorBar(time(range), mean_mci(range), se_mci(range), ...
        'lineProps', {'Color', col_mci, 'LineWidth',2});
    
    if contrIND <= n_reg
        % WITHIN-GROUP significance bars (grey group color p<.10, darker p<.05 on top)
        plot(time(range), plotsign_yas_10(range), '-', 'Color', col_ya,     'LineWidth', 8);
        plot(time(range), plotsign_yas_05(range), '-', 'Color', col_ya_05,  'LineWidth', 6);
        plot(time(range), plotsign_oas_10(range), '-', 'Color', col_oa,     'LineWidth', 8);
        plot(time(range), plotsign_oas_05(range), '-', 'Color', col_oa_05,  'LineWidth', 6);
        plot(time(range), plotsign_mci_10(range), '-', 'Color', col_mci,    'LineWidth', 8);
        plot(time(range), plotsign_mci_05(range), '-', 'Color', col_mci_05, 'LineWidth', 6);
        

        % BETWEEN-GROUP significance bars (grey p<.10, dark grey p<.05)
        % YA > OA
        plot(time(range), plotsign_10(1,range), '-', 'Color', col_between_10, 'LineWidth', 8);
        plot(time(range), plotsign_05(1,range), '-', 'Color', col_between_05, 'LineWidth', 6);
        
        % YA > MCI
        plot(time(range), plotsign_10(2,range), '-', 'Color', col_between_10, 'LineWidth', 8);
        plot(time(range), plotsign_05(2,range), '-', 'Color', col_between_05, 'LineWidth', 6);
        
        % OA > MCI
        plot(time(range), plotsign_10(3,range), '-', 'Color', col_between_10, 'LineWidth', 8);
        plot(time(range), plotsign_05(3,range), '-', 'Color', col_between_05, 'LineWidth', 6);
        
        % add text labels at the left edge of the plot for each comparison
        % get text y positions for this contrast
        ty = text_ypos(con_name{contrIND});
        
   
        cg_real = store_group_glm{contrIND}.cg;
        
        % YA vs OA (grcon = 1)
        if any(~isnan(plotsign_05(1,range))) || any(~isnan(plotsign_10(1,range)))
            effect_val = nanmean(cg_real(1, min(tw):max(tw)));
            if effect_val > 0
                label1 = 'YAs > OAs';
            else
                label1 = 'OAs > YAs';
            end
            text(0.25, ty(1), label1, 'FontName', 'Serif', 'FontSize', 45, 'Color', col_between_05);
        end
        
        % YA vs MCI (grcon = 2)
        if any(~isnan(plotsign_05(2,range))) || any(~isnan(plotsign_10(2,range)))
            effect_val = nanmean(cg_real(2, min(tw):max(tw)));
            if effect_val > 0
                label2 = 'YAs > MCI';
            else
                label2 = 'MCI > YAs';
            end
            text(0.25, ty(2), label2, 'FontName', 'Serif', 'FontSize', 45, 'Color', col_between_05);
        end
        
        % OA vs MCI (grcon = 3)
        if any(~isnan(plotsign_05(3,range))) || any(~isnan(plotsign_10(3,range)))
            effect_val = nanmean(cg_real(3, min(tw):max(tw)));
            if effect_val > 0
                label3 = 'OAs > MCI';
            else
                label3 = 'MCI > OAs';
            end
            text(0.25, ty(3), label3, 'FontName', 'Serif', 'FontSize', 45, 'Color', col_between_05);
        end
        ylim(sign_ylim(contrIND, :));
        xlim([-0.5 2.5]);

        if contrIND == 1  % only for emotionality
            plot(time(range), plotsign_rt_10(range), '-', 'Color', [0.6 0.9 0.4], 'LineWidth', 8); % p<.10 light green
            plot(time(range), plotsign_rt_05(range), '-', 'Color', [0.1 0.6 0.1], 'LineWidth', 6); % p<.05 dark green
            if any(~isnan(plotsign_rt_05(range))) || any(~isnan(plotsign_rt_10(range)))
                text(0.25, yp.between(3) - 0.02, 'RT (MCI)', 'FontName', 'Serif', ...
                     'FontSize', 45, 'Color', [0.1 0.6 0.1]);
            end

        end

        if contrIND == 2
            plot(time(range), plotsign_dp_10(range), '-', 'Color', [0.6 0.9 0.4], 'LineWidth', 8);
            plot(time(range), plotsign_dp_05(range), '-', 'Color', [0.1 0.6 0.1], 'LineWidth', 6);
            if any(~isnan(plotsign_dp_05(range))) || any(~isnan(plotsign_dp_10(range)))
                text(0.25, yp.between(1) - 0.1, 'Recog. memory', 'FontName', 'Serif', ...
                'FontSize', 45, 'Color', [0.1 0.6 0.1]);
            end
        end
    end    

    hold off

    set(gca, 'FontName', 'Serif', 'FontSize', 50);
    ax = gca;
    ax.XAxis.LineWidth = 2;
    ax.YAxis.LineWidth = 2;
    ax.XAxis.Color = 'k';
    ax.YAxis.Color = 'k';
    xlabel('Time (s)', 'FontName', 'Serif', 'FontSize', 50);
    ylabel('Beta weights (PD, z-scored)', 'FontName', 'Serif', 'FontSize', 50);
    title(display_names{contrIND}, 'FontName', 'Serif', 'FontSize', 50, 'FontWeight', 'normal');
    lgd = legend({'YAs', 'OAs', 'MCI'}, ...
        'FontName', 'Serif', 'FontSize', 50, ...
        'Box', 'off', ...
        'Location', 'eastoutside');
    title(lgd, 'Group', 'FontSize', 50, 'FontName', 'Serif', 'FontWeight', 'normal');
    lgd.PlotChildren(1).LineWidth = 4;  % YAs line
    lgd.PlotChildren(2).LineWidth = 4;  % OAs line  
    lgd.PlotChildren(3).LineWidth = 4;  

    xlim([-0.2 2.5]);

    if contrIND <= n_reg
        ylim(sign_ylim(contrIND, :));
    else
        ylim([-0.50 0.90]);  % intercept range
    end

    exportgraphics(gcf, ...
        fullfile(pathsave, ['GLM_' con_name{contrIND} '.png']), ...
        'Resolution', 300);
end

fprintf('\n=== NaN diagnostics - cgb ===\n\n');
for contrIND = 1:n_reg
    n_nan = sum(isnan(squeeze(cgb(contrIND, 1, 1, :))));
    fprintf('  contrIND %d (%s): %d NaN timepoints in cgb\n', contrIND, reg_names{contrIND}, n_nan);
end

fprintf('\ndat1 NaN subjects per contrast:\n');
for contrIND = 1:2
    dat1 = zeros(length(ids_with_both_sessions), length(time));
    for i = 1:length(ids_with_both_sessions)
        dat1(i,:) = storedat(i).glm.con{1, contrIND};
    end
    bad_subs = find(any(isnan(dat1), 2));
    fprintf('  contrIND %d (%s): %d subjects with NaN betas\n', contrIND, reg_names{contrIND}, length(bad_subs));
    for k = 1:length(bad_subs)
        fprintf('    Subject %d\n', ids_with_both_sessions(bad_subs(k)));
    end
end

fprintf('\n=== within-group significance diagnostics ===\n\n');

groups = {'YAs', 'OAs', 'MCI'};
group_idx = {yas, oas, mci};

for contrIND = 1:length(reg_names)
    fprintf('--- Regressor: %s ---\n', reg_names{contrIND});
    
    dat1 = zeros(length(ids_with_both_sessions), length(time));
    for i = 1:length(ids_with_both_sessions)
        dat1(i,:) = storedat(i).glm.con{1, contrIND};
    end
    
    for g = 1:3
        grp = group_idx{g};
        grp_name = groups{g};
        mean_grp = nanmean(dat1(grp, :));
        real_max = max(mean_grp(min(tw):max(tw)));
        real_min = min(mean_grp(min(tw):max(tw)));
        
        datdum = squeeze(nanmean(cib(grp, contrIND, :, :), 1));
        for rep = 1:length_randpermute_indiv
            datmax(rep) = max(datdum(rep, min(tw):max(tw)));
            datmin_perm(rep) = min(datdum(rep, min(tw):max(tw)));
        end
        
        p90 = prctile(datmax, 90);  p95 = prctile(datmax, 95);
        p10 = prctile(datmin_perm, 10); p05 = prctile(datmin_perm, 5);
        
        % check both positive and negative effects
        if real_max > p95 || real_min < p05
            sig = 'SIGNIFICANT p<0.05';
        elseif real_max > p90 || real_min < p10
            sig = 'TREND p<0.10';
        else
            sig = 'not significant';
        end
        
        fprintf('%s: max=%.4f (90th=%.4f, 95th=%.4f) min=%.4f (10th=%.4f, 5th=%.4f) → %s\n', ...
            grp_name, real_max, p90, p95, real_min, p10, p05, sig);
        clear datmax datmin_perm datdum

    end
    fprintf('\n');
end

fprintf('\n=== between-group significance diagnostics ===\n\n');
group_con_names = {'YA > OA', 'YA > MCI', 'OA > MCI'};

for contrIND = 1:length(reg_names)
    fprintf('--- Regressor: %s ---\n', reg_names{contrIND});
    cg_real = store_group_glm{contrIND}.cg;
    
    for grcon = 1:3
        real_max = max(cg_real(grcon, min(tw):max(tw)));
        real_min = min(cg_real(grcon, min(tw):max(tw)));
        
        for rep = 1:length_randpermute
            datmax(rep) = max(squeeze(cgb(contrIND, rep, grcon, min(tw):max(tw))));
            datmin_perm(rep) = min(squeeze(cgb(contrIND, rep, grcon, min(tw):max(tw))));
        end
        
        p90 = prctile(datmax, 90); p95 = prctile(datmax, 95);
        p10 = prctile(datmin_perm, 10); p05 = prctile(datmin_perm, 5);
        
        if real_max > p95 || real_min < p05
            sig = 'SIGNIFICANT p<0.05';
        elseif real_max > p90 || real_min < p10
            sig = 'TREND p<0.10';
        else
            sig = 'not significant';
        end
        
        fprintf('%s: max=%.4f (90th=%.4f, 95th=%.4f) min=%.4f (10th=%.4f, 5th=%.4f) → %s\n', ...
            group_con_names{grcon}, real_max, p90, p95, real_min, p10, p05, sig);
        clear datmax datmin_perm
    end
    fprintf('\n');
end

fprintf('\n=== RT effect significance (emotionality only) ===\n\n');
cg_real = store_group_glm{1}.cg;
real_max = max(cg_real(4, min(tw):max(tw)));
real_min = min(cg_real(4, min(tw):max(tw)));

for rep = 1:length_randpermute
    datmax(rep) = max(squeeze(cgb(1, rep, 4, min(tw):max(tw))));
    datmin_perm(rep) = min(squeeze(cgb(1, rep, 4, min(tw):max(tw))));
end

p90 = prctile(datmax, 90); p95 = prctile(datmax, 95);
p10 = prctile(datmin_perm, 10); p05 = prctile(datmin_perm, 5);

if real_max > p95 || real_min < p05
    sig = 'SIGNIFICANT p<0.05';
elseif real_max > p90 || real_min < p10
    sig = 'TREND p<0.10';
else
    sig = 'not significant';
end

fprintf('RT effect: max=%.4f (90th=%.4f, 95th=%.4f) min=%.4f (10th=%.4f, 5th=%.4f) → %s\n', ...
    real_max, p90, p95, real_min, p10, p05, sig);

fprintf('\n=== dprime effect significance (recognition only) ===\n\n');
cg_real = store_group_glm{2}.cg;
real_max = max(cg_real(4, min(tw):max(tw)));
real_min = min(cg_real(4, min(tw):max(tw)));

for rep = 1:length_randpermute
    datmax(rep) = max(squeeze(cgb(2, rep, 4, min(tw):max(tw))));
    datmin_perm(rep) = min(squeeze(cgb(2, rep, 4, min(tw):max(tw))));
end

p90 = prctile(datmax, 90); p95 = prctile(datmax, 95);
p10 = prctile(datmin_perm, 10); p05 = prctile(datmin_perm, 5);

if real_max > p95 || real_min < p05
    sig = 'SIGNIFICANT p<0.05';
elseif real_max > p90 || real_min < p10
    sig = 'TREND p<0.10';
else
    sig = 'not significant';
end

fprintf('Dprime effect: max=%.4f (90th=%.4f, 95th=%.4f) min=%.4f (10th=%.4f, 5th=%.4f) → %s\n', ...
    real_max, p90, p95, real_min, p10, p05, sig);

fprintf('\n=== significant time windows ===\n\n');
group_names_sw = {'YAs', 'OAs', 'MCI'};
field_names_05 = {'yas_05', 'oas_05', 'mci_05'};

for contrIND = 1:n_reg
    fprintf('--- %s ---\n', reg_names{contrIND});
    for g = 1:3
        sig = store_sigwindow(contrIND).(field_names_05{g});
        sig_times = time(~isnan(sig));
        if ~isempty(sig_times)
            fprintf('%s: %.2f to %.2f s (duration: %.2f s)\n', ...
                group_names_sw{g}, min(sig_times), max(sig_times), max(sig_times)-min(sig_times));
        else
            fprintf('%s: not significant\n', group_names_sw{g});
        end
    end
    fprintf('\n');
end