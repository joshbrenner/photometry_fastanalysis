function mice_arr = photometry_allmice_analyze_trials_Hemo_pretrialBL

%this version is used to generate mouse arrays with the 405 channel used as
%a baseline for the general mismatch charts.

path = 'C:\Users\Tuica\Documents\MATLAB\photometry processing\Single_Fiber_Josh-221006-121345\';
% names = {'m3_day1','m3_day2a','m3_day2b','m4_day1','m4_day2','m735_day1','m735_day2','m772_day2a','m772_day2b'};
% onset_cell = [1,2537,24235,102519,2073,22786,1288,487,7284]; 
% offset_cell = [3973520,485159,1801190,3649790,2792470,5874560,1466900,2792470,1115310];

names = {'m4_day1','m4_day2','m735_day1','m735_day2','m772_day2a','m772_day2b'}; %different sessions and mice
onset_cell = [102519,2073,22786,1288,487,7284];%This device has limited digital inputs so we have to mark our onset and offset manually
offset_cell = [3649790,2792470,5874560,1466900,2792470,1115310];
mice_array = cell([1,length(names)]);

for exp = 1:length(names)
   mice = analyzethosemice
mice_arr{1,exp} = mice;

end


    function mice = analyzethosemice
    
    %% Load and preprocess photometry data - we will have to adjust onset/offset times for each experiment since we don't have a shutter.
    photometry_block_name = strcat(path, names{exp});
    
    photometry_data = TDTbin2mat(photometry_block_name);
    raw_diode = photometry_data.streams.Wav2.data;
    Fs = photometry_data.streams.x470A.fs; % Sampling rate
    onset = onset_cell(exp); %we don't have a clear on signal from the recording rig - just toss the first few thousand frames and look for the first up in the diff diode signal.
    offset = offset_cell(exp); %and the last off..
    photLim = 2575;
    
    ds_rate =60;
    stim_len_lim = (Fs*15)/ds_rate; %Any stims longer than 30 seconds don't have a replay, so we discard them 
    
    
    sig470 = photometry_data.streams.x470A.data(onset:offset); %from calcium signal
    sig470  = nanmean(reshape([sig470(:); nan(mod(-numel(sig470),ds_rate),1)],ds_rate,[])); %downsample by a factor of ds_rate to smooth some of that noise.
    sig405 = photometry_data.streams.x405A.data(onset:offset); %from hemodynamic signal
    sig405  = nanmean(reshape([sig405(:); nan(mod(-numel(sig405),ds_rate),1)],ds_rate,[])); %downsample by a factor of ds_rate to smooth some of that noise.
    
    sig_run = photometry_data.streams.Wav1.data(onset:offset); %from linear encoder
    sig_run = nanmean(reshape([sig_run(:); nan(mod(-numel(sig_run),ds_rate),1)],ds_rate,[])); %downsample by a factor of ds_rate to smooth some of that noise.
    sig_diode = photometry_data.streams.Wav2.data(onset:offset); %from photodiode (ie trial onset/offset)
    sig_diode = nanmean(reshape([sig_diode(:); nan(mod(-numel(sig_diode),ds_rate),1)],ds_rate,[])); %downsample by a factor of ds_rate to smooth some of that noise.
    
    %% Subtract hemodynamic signal 
    controlFit = [];
    binind = ceil(linspace(1,length(sig470),10));
    binind(end) = binind(end)+1;
    for i = 1:numel(binind)-1

    reg = polyfit(sig405(binind(i):binind(i+1)-1), sig470(binind(i):binind(i+1)-1), 1);
    a = reg(1);
    b = reg(2);
    controlFit = [controlFit (a.*sig405(binind(i):binind(i+1)-1) + b)]; %linear fit hemo signal to phot signal
    end
    controlFit = movmean(controlFit,20000,'omitnan','Endpoints','shrink');
    reg = polyfit(sig405,sig470,1);
    a = reg(1);
    b = reg(2);
    controlFit2= a.*sig405 + b;
    
    controlFit3 = (controlFit.*controlFit2)./nanmean(controlFit2);

    df_photometry = (sig470-controlFit3)./controlFit3; %f
    
    
    

    %% Find trial times
    
    binary_trial_times = sig_diode > photLim; %mV signal that indicates diode is on
    diff_bin_tt = diff(binary_trial_times(1:end));
    stimLimits = [find(diff_bin_tt>.5);find(diff_bin_tt<-.5)]; %first row is onset, second is offset
    stim_len = stimLimits(2,:)-stimLimits(1,:)+1; %add one since this otherwise doesn't count the first frame.

    %Get rid of the first pair - signal tends to be contaminated by changing
    %LED strength etc.
    
    stimLimits(:,1:2) = [];
    stim_len(1:2) = [];
    stimLimits(:,end-1:end)= [];
    stim_len(end-1:end)= [];


    if exp == 10 %This experiment has some issues detecting trials without a replay. We go in and identify them manually.
    rem_mat = [37 38 43 44 47 48 51 54 55 68 73 84 85 86 111 116 119 120 173 194 195 196 201 204 239 270 277 294 313 324 325 344 345 356];
    stim_len(rem_mat) = [];
    stimLimits(:,rem_mat) = [];


    else
      
% flag = zeros([1 length(stim_len)]);
% for i = 1:length(stim_len)
%     if i == 1
%         if stim_len(i) > stim_len_lim*1.3 && (abs(stim_len(i)-stim_len(i+1))/stim_len(i)) > .1
%             flag(i) = 1;
%         end
%     elseif i == length(stim_len)
%         if stim_len(i) > stim_len_lim*1.3 && (abs(stim_len(i)-stim_len(i-1))/stim_len(i)) > .1
%             flag(i) = 1;
%         end
%     elseif stim_len(i) > stim_len_lim*1.3 && (abs(stim_len(i)-stim_len(i+1))/stim_len(i)) > .1 && (abs(stim_len(i)-stim_len(i-1))/stim_len(i)) > .1
%         flag(i) = 1;
%     end
% end
% %    
%     stim_len(flag) = [];
%     stimLimits(flag) = [];

    stimLimits(:,stim_len>=stim_len_lim*1.3) = []; %gets rid of most too-long trials
    stim_len(stim_len>=stim_len_lim*1.3) = []; 
    

    stim_removal = stim_len >= stim_len_lim;
    
    
    for i = 1:numel(stim_len) 
        if stim_removal(i) && mod(i,2)==0
            stim_removal(i-1) = 1;
        elseif stim_removal(i) && stim_len(i) < stim_len_lim*1.3 && i ~= length(stim_len)%Make sure we don't remove only one of the two paired trials by accident -- 
            %of course, if you have a very long self-gen trial it will not generate a paired exo trial, so check for that as well
            stim_removal(i+1) = 1;
        end
    end
    
    stimLimits(:,stim_removal) = []; 
    stim_len(:,stim_removal) = []; 
    end
    stim_len_pairs = [];
    stim_len_diff = [];
    %Get rid of pairs with more than +/-5% difference in trial length
    stim_len_pairs = [stim_len(1:2:end);stim_len(2:2:end)];
    stim_len_diff = (stim_len_pairs(2,:) - stim_len_pairs(1,:))./stim_len_pairs(1,:);
    stim_len_diff_size = (stim_len_pairs(2,:) - stim_len_pairs(1,:));
    disp("maximum replay/self-gen length difference = ")
    disp(max(abs(stim_len_diff))) %check that there's nothing too big - could indicate a mispairing.
    
    err_removal_mat = zeros([1,numel(stim_len)]);
    for i = 1:numel(stim_len_diff)
        if abs(stim_len_diff(i)) >= .05 && abs(stim_len_diff_size(i)) > 2
            err_removal_mat((2*i)-1:2*i) = 1;
        end
    end
    
    stimLimits(:,logical(err_removal_mat)) = [];
    stim_len(logical(err_removal_mat)) = [];
    

    %% Find mean responses and baselines for each epoch and calculate dff
% BLtime = round(Fs/(3*ds_rate));
% 
% phot_resp(1:BLtime) = NaN;
% for i = BLtime+1:length(sig470)
%     phot_resp(i) = (sig470(i)-nanmean(sig470(i-BLtime:i)))./nanmean(sig470(i-BLtime:i));
% end

    %% Find mean responses and baselines for each epoch and calculate dff
%     BLtime = round((2*Fs)/(ds_rate)); %take a 2 second baseline
%     resp_delay = ceil((Fs*.3)/ds_rate);
%     for i = 1:size(stimLimits,2)
%     phot_resp(i) = nanmean(sig470(stimLimits(1,i)+resp_delay:stimLimits(2,i)+resp_delay)); %iteratively take responses between the stim limits defined by the diode
%     phot_BL(i) = nanmean(sig470(stimLimits(1,i)-BLtime:stimLimits(1,i)-1)); %take a baseline of 1 second before the stimulus onset
%     phot_hemo(i) = nanmean(controlFit(stimLimits(1,i)+resp_delay:stimLimits(2,i)+resp_delay)); 
%     end
% 
%     %smooth the baseline
%     phot_BL_smooth = movmean(phot_BL,10,'omitnan','Endpoints','shrink');
%     phot_dff = (phot_resp-phot_BL_smooth)./abs(phot_BL_smooth); %calculate dff
%     phot_dff2 = phot_dff - controlFit;
%     
    % phot_dff_pair = [phot_dff(1:2:end);phot_dff(2:2:end)];
    % mean(phot_dff_pair,2)
    % phot_lim = mean(phot_BL)*.25; %Remove any trials with unusually large baselines
    % phot_removal = phot_BL > phot_lim+mean(phot_BL) | phot_BL < mean(phot_BL) - phot_lim;
%     phot_resp(1:BLtime) = NaN;
%     for i = BLtime+1:length(sig470)
%         phot_resp(i) = (sig470(i)-nanmean(sig470(i-BLtime:i)))./nanmean(sig470(i-BLtime:i));
%     end
    

    %% Find running speed during each epoch
    
    dif_sig_run  = [0 diff(sig_run)]; %convert our running signal into a change in voltage
    dif_sig_run((dif_sig_run)>300|dif_sig_run<-20) = NaN; %If we have a large jump this is either noise or the encoder switching 0 <-> +5 Volts as the mouse completes a circuit. This has to be tuned for each mouse because there
    %are inconsistent sized artifacts in the encoder.
    
    for i = 1:2:size(stimLimits,2) %find the running speed across each epoch
        rs_arr{1,ceil(i/2)} = (dif_sig_run(stimLimits(1,i):stimLimits(2,i)));
        rs_arr{2,ceil(i/2)} = (dif_sig_run(stimLimits(1,i+1):stimLimits(2,i+1)));
    end
    
    for i = 1:size(rs_arr,2) %Find the mismatch between rs in self-gen and exogenous epochs
        stretchSize = size(rs_arr{2,i}); %resize the self-generated to the exogenous for easier comparison
%         mm_arr{1,i} = (imresize(rs_arr{1,i},stretchSize,'nearest') - rs_arr{2,i})./(imresize(rs_arr{1,i},stretchSize,'nearest'));
        NaN_arr = (abs(imresize(rs_arr{1,i},stretchSize,'nearest')) < 2.5) & (abs(rs_arr{2,i}) < 2.5); %get rid of low mm frames where the mouse isn't running at all in either condition...
%         mm_arr{1,i}(NaN_arr) = NaN; %This creates a weird spike around 0
        %negative mm indicates a visual stimulus moving slower than expected; positive means faster than expected
    end
    
    %% Sort frames by mismatch
    all_mm = [];
    exo_sig = [];
    self_sig = [];
    rs_sig_exo = [];
    rs_sig_self = [];
    resp_delay = ceil((Fs*.15)/ds_rate);
    resp_delay_trial = ceil((Fs*1.1)/ds_rate);
    binSelf = [];
    binExo = [];
%     
%     for i = 1:size(rs_arr,2)  %for 405 bl
%         all_mm = [all_mm mm_arr{i}]; %take all mismatch values
%         rs_sig_exo = [rs_sig_exo ((dif_sig_run(stimLimits(1,2*i):stimLimits(2,2*i))))];
%         rs_sig_self = [rs_sig_self ((imresize(dif_sig_run(stimLimits(1,(2*i)-1):stimLimits(2,(2*i)-1)), [1 stim_len(2*i)], 'nearest')))];
%         exo_sig = [exo_sig ((df_photometry(stimLimits(1,2*i)+resp_delay:stimLimits(2,2*i)+resp_delay)))];
%         self_sig = [self_sig ((imresize(df_photometry(stimLimits(1,(2*i)-1)+resp_delay:stimLimits(2,(2*i)-1)+resp_delay), [1 stim_len(2*i)], 'nearest')))];
%       
%     end
% %     
%     
%     for i = 1:size(rs_arr,2) %for a moving window bl
%         all_mm = [all_mm mm_arr{i}]; %take all mismatch values
%         rs_sig_exo = [rs_sig_exo ((dif_sig_run(stimLimits(1,2*i):stimLimits(2,2*i))))];
%         rs_sig_self = [rs_sig_self ((imresize(dif_sig_run(stimLimits(1,(2*i)-1):stimLimits(2,(2*i)-1)), [1 stim_len(2*i)], 'nearest')))];
%         exo_sig = [exo_sig ((phot_resp(stimLimits(1,2*i)+resp_delay:stimLimits(2,2*i)+resp_delay)))];
%         self_sig = [self_sig ((imresize(phot_resp(stimLimits(1,(2*i)-1)+resp_delay:stimLimits(2,(2*i)-1)+resp_delay), [1 stim_len(2*i)], 'nearest')))];
%       
%     end
%     
%     
% 
    for i = 1:size(rs_arr,2) %for a trial based 405 bl 
    exo_sig(i) = nanmean(df_photometry(stimLimits(1,2*i):stimLimits(2,2*i)+resp_delay_trial));  
    self_sig(i) = nanmean(df_photometry(stimLimits(1,(2*i)-1):stimLimits(2,(2*i)-1)+resp_delay_trial)); 
    rs_sig_exo(i) = nanmean(dif_sig_run(stimLimits(1,2*i):stimLimits(2,2*i)));
    rs_sig_self(i) = nanmean(dif_sig_run(stimLimits(1,(2*i)-1):stimLimits(2,(2*i)-1)));

    exo_sigf{i} = (df_photometry(stimLimits(1,2*i)-5:stimLimits(2,2*i)+resp_delay_trial));  
    self_sigf{i} = (df_photometry(stimLimits(1,(2*i)-1)-5:stimLimits(2,(2*i)-1)+resp_delay_trial)); 
    
    if rs_sig_exo(i) < 0
        rs_sig_exo(i) = 0;
    end
    if rs_sig_self(i) < 0
        rs_sig_self(i) = 0;
    end
    all_mm(i) = (rs_sig_self(i) - rs_sig_exo(i))./(rs_sig_self(i)+rs_sig_exo(i));
    end
    


    [mm_sort, sort_ind] = sort(all_mm);
    sort_self = self_sig(sort_ind);
    sort_exo = exo_sig(sort_ind);
    rs_sig_exo = rs_sig_exo(sort_ind);
    rs_sig_self = rs_sig_self(sort_ind);
    sort_selff = self_sigf(sort_ind);
    sort_exof = exo_sigf(sort_ind);

    
binSelf = movmean(sort_self, 1,'omitnan');
binExo = movmean(sort_exo,1,'omitnan');
x = movmean(mm_sort,10,'omitnan');

    
    figure
    plot(x, binSelf)
    hold on
    plot(x, binExo)
    
    mice.name = names{exp};
    mice.mmsort = mm_sort;
    mice.sortself = sort_self;
    mice.sortexo = sort_exo;
    mice.sortind = sort_ind;
    mice.rs_sig_exo = rs_sig_exo;
    mice.rs_sig_self = rs_sig_self;
    mice.sortselff = sort_selff;
    mice.sortexof = sort_exof;
    
    
    % 
    % figure
    % scatter(1:numel(sort_exo),sort_exo)
    % hold on
    % x = -.2:.1:.3;
    % plot(x,x)
    end

end



