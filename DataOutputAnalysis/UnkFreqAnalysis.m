%% ------ NOTES ------ %%
% Nicole Ortega @ 7/28/2017
%
% Revision of Tapping Data that allows for trials testing any combination of perturbation
% frequencies. Conditions tested (combination of perturbation frequencies)
% are determined by readDevFreq_rev1().
% 
% This version makes sure to check input to microprocessor with buzzer
% output to prove accuracy of input. This version organizes data by subject
% and differentiates between rhythm and nonrhythm trials. This
% version breaks down Tapping() into three seperate functions.
%
% Individual trials have their own output file. 
% Allows for multiple subjects and trials to be read. Threshold is constant. 
%
% Fix: Tapping function is broken up into three new functions --> 
% StoreLogData(), AvgMetDev(), and CalcDeviation()
% Fix: Each trial has its own output file. Previously, one file contained
% all trials. 
% Fix: Calculates the FFT based on new discrete perturbations 
% Fix: Threshold value is constant
% Fix: Added for loops  to enable as many subjects and trials as necessary
% NOTE: when Nans are still present in signal, the FFT does not work
% Haptic environment: pushing down against button
%
%% ------ Clear ------ %%
clear all; 
close all;

%% ------ Constants ------ %%
sample_freq  = 1000;                % 1k Hz
met_period = .333;
met_freq = 1/met_period;                  % freq of steady metronome in Hz
fundamental_period = 80;            %events per cycle, AKA perturbations per period
fundamental_freq = 1/fundamental_period;    %cycles per event
repeats = 2;                        %number of times a complete trial will repeat
thresh_val = -105;                  % Estimating threshold value

% these prime multiples are approx log spaced. dropping 2 to keep it to 7 
dev_harmonics = [3 5 7 11 17 23 31]; % 23 29 ** 37 max for period = 80
% phase shifts were randomly chosen using rand(7,1)
dev_phase_shift = [0.5291 0.0409 0.2208 0.8799 0.0536 0.7651 0.1022];

%% ------ Variables ------ %%
startSubj = 1;
num_subjects = 8;                   % number of subjects
date = '7.12.2017';                 % date for experiment
rest =  5;                          % in seconds

% Tested frequencies
conditions = dev_harmonics * fundamental_freq;
condition_set = [conditions(1:3); conditions(4:5) 0; conditions(6) 0 0];

fixed_trials = [8 12 13 0 0;... 
    0 0 0 0 0;...
    4 5 7 0 0;...
    4 0 0 0 0;...
    6 0 0 0 0;...
    6 0 0 0 0;...
    0 0 0 0 0;...
    0 0 0 0 0;...
    0 0 0 0 0;...
    5 0 15 0 0;...
    14 18 0 0 0;...
    9 12 0 0 0;...
    0 0 0 0 0];
unsuccessful_trials = [0 0 0 0 0;... 
    8 11 0 0 0;...
    0 0 0 0 0;...
    0 0 0 0 0;...
    0 0 0 0 0;...
    0 0 0 0 0;...
    16 17 18 0 0;...
    0 0 0 0 0;...
    4 5 0 0 0;...
    8 0 0 0 0;...
    0 0 0 0 0;...
    0 0 0 0 0;...
    18 0 0 0 0];

signal_len = 260;

% initializing variables 
reps = zeros(size(condition_set,1), 1);
reps_old = zeros(size(condition_set,1), 1);
reps_r = zeros(size(condition_set,1), 1);
reps_nr = zeros(size(condition_set,1), 1);
perc_tapsvsbuzz = nan(size(condition_set,1), num_subjects);
perc_tapsthrown = nan(size(condition_set,1), num_subjects);

subject_data = struct;
rhythm_data = struct;

%%
% ---- Open Files and Store Data ---- % 
% Store raw data (metronome times, buzzer times, tapping times) and find
% average deviation from steady metronome. Data will later be used to
% determine the buzzer/tapping deviation from steady metronome
for nSubj = startSubj:num_subjects
%     if nSubj == 6
%         continue;
%     end
    rhythm_tap_file = csvread(strcat('/',date, '/', date,'_Subject_',num2str(nSubj),'_Rhythm.txt'));
    num_trials = length(rhythm_tap_file);
    % Frequencies tested per subject in the order they appear
    freq_devs_by_subject = readDevFreq_rev1(nSubj, date, num_trials);
    avg_met_dev = nan(num_trials, 1);
    rhythm_reps = zeros(length(condition_set), num_trials, 2);
   
    for nTr = 1:num_trials
        if find(nTr == nonzeros(unsuccessful_trials(nSubj,:)))
            continue;
        end
        % initialize variables
        n_met_all = nan(signal_len,1); n_buzz_all = nan(signal_len,1); n_tap_all = nan(signal_len,1);
        buzz_dev_all = nan(signal_len,1); tap_dev_all = nan(signal_len,1);
        
        % extract deviation freq being tested during trial
        freq_devs_by_trial = nonzeros(freq_devs_by_subject(nTr,:));
        % determine condition being tested
        n_cond = 4 - length(freq_devs_by_trial);
       
        
        % read input and output file
        log_file = csvread(strcat('/',date, '/', date,'_Subject_',num2str(nSubj),'.',num2str(nTr),'_Res.txt'));
        input_file = csvread(strcat('/',date, '/', date,'_Subject_',num2str(nSubj),'.',num2str(nTr),'_DecPerts.txt'));
        % skip unsuccessful trials
        if find(nTr == nonzeros(fixed_trials(nSubj,:)))
            log_file = fixTrial(log_file);
        end
        log_file = log_file(1:end-2,:); %cutting out extra data points at the end of file
        
        % Extract steady metronome times, buzzer times, and tapping times
        % Store in vectors of the same length
        [ n_met, n_buzz, n_tap ] = ExtractLogData_rev0(met_period, log_file, freq_devs_by_trial, thresh_val);
        n_met_all(1:length(n_met)) = n_met;  n_buzz_all(1:length(n_buzz)) = n_buzz; n_tap_all(1:length(n_tap)) = n_tap;
        
        % Calculate the average delay during a steady metronome 
        [ subject_data(nSubj).avg_met_dev(nTr) ] = AvgMetDev_rev0(met_period, n_buzz, n_tap);

        rhythm_cond = rhythm_tap_file(nTr) + 1;
        rhythm_reps(n_cond, nTr, rhythm_cond) = 1;
        rep = sum(rhythm_reps(n_cond, :, rhythm_cond));
            
        subject_data(nSubj).rhythm_condition(rhythm_cond).trials_by_condition(n_cond).n_met(:,rep) = n_met_all;
        subject_data(nSubj).rhythm_condition(rhythm_cond).trials_by_condition(n_cond).n_buzz(:,rep) = n_buzz_all;
        subject_data(nSubj).rhythm_condition(rhythm_cond).trials_by_condition(n_cond).n_tap(:,rep) = n_tap_all;
        
        close all;
        
    end
    
end
close all;
%% 

for nSubj = startSubj:size(subject_data, 2)
    rhythm_tap_file = csvread(strcat('/',date, '/', date,'_Subject_',num2str(nSubj),'_Rhythm.txt'));
    num_trials = length(rhythm_tap_file);
    % Frequencies tested per subject in the order they appear
    freq_devs_by_subject = readDevFreq_rev1(nSubj, date, num_trials);
    
    for rhCond = 1:size(subject_data(nSubj).rhythm_condition, 2)
        
        for trCond = 1:size(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition, 2)
            
            for rep = 1: size(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met, 2)
                %initializing variables
                buzz_dev_all = nan(signal_len,1); tap_dev_all = nan(signal_len,1);
                
                % extract deviation freq being tested during trial
                freq_devs_by_trial = nonzeros(freq_devs_by_subject(nTr,:));
                avg_met_dev = mean(nonzeros(subject_data(nSubj).avg_met_dev(:)));

                n_met = subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met(:,rep);
                n_buzz = subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_buzz(:,rep);
                n_tap = subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_tap(:,rep);

                % take out nan values
                n_met=(n_met(~isnan(n_met)));
                n_buzz=(n_buzz(~isnan(n_buzz)));
                n_tap=(n_tap(~isnan(n_tap)));


                [ buzz_dev, tap_dev, percent_taps ] = CalcDeviation_rev0(n_met, n_buzz, n_tap, avg_met_dev);
                buzz_dev_all(1:length(n_met)) = buzz_dev; tap_dev_all(1:length(n_met)) = tap_dev;
                subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev(:,rep) = buzz_dev_all;
                subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev(:,rep) = tap_dev_all;
                
                close all;
                %statistic of tapping/buzzer percentage (raw data)
                p_tb(1:2) = [perc_tapsvsbuzz(trCond, nSubj)  length(n_tap)/length(n_buzz)];
                perc_tapsvsbuzz(trCond, nSubj) = nanmean(p_tb);

                %statistics of percentage taps kept to find deviation
                p_tt(1:2) = [perc_tapsthrown(trCond, nSubj) percent_taps]; 
                perc_tapsthrown(trCond, nSubj) = nanmean(p_tt);

            end
        end
    end
end

close all;
%%
% Calculate nanmean per subject 
for nSubj = startSubj:size(subject_data, 2)
    for rhCond = 1:size(subject_data(nSubj).rhythm_condition, 2)
        for trCond = 1:size(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition, 2)
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met_avg = nanmean(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met,2);
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_buzz_avg = nanmean(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_buzz,2);
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_tap_avg = nanmean(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_tap,2);
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev_avg = nanmean(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev,2);
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev_avg = nanmean(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev,2);
        end
    end

end

for nSubj = startSubj:size(subject_data, 2)
    for rhCond = 1:size(subject_data(nSubj).rhythm_condition, 2)
        for trCond = 1:size(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition, 2)
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met_avg = nanmean([subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met_avg(42:121) subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met_avg(122:201)],2);
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_buzz_avg = nanmean([subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_buzz_avg(42:121) subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_buzz_avg(122:201)],2);
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_tap_avg = nanmean([subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_tap_avg(42:121) subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).n_tap_avg(122:201)],2);
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev_avg = nanmean([subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev_avg(42:121) subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev_avg(122:201)],2);
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev_avg = nanmean([subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev_avg(42:121) subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev_avg(122:201)],2);
        end
    end
end
%%
% Calculate FFT for individual subjects 
for nSubj = startSubj:size(subject_data, 2)
    for rhCond = 1:size(subject_data(nSubj).rhythm_condition, 2)
        for trCond = 1:size(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition, 2)
            buzz_dev = subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev_avg;
            tap_dev = subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev_avg;
            [FFT_buzz, FFT_tap, mag_buzz, mag_tap, angle_buzz, angle_tap, f_axis] = CalcFFT_rev1(fundamental_period, 1, buzz_dev, tap_dev);
            str = sprintf('FFT of Frequency Deviations: %.2f %.2f %.2f', freq_devs_by_trial);
            title(str);
            xlabel('Frequency (cycles per event)');
            ylabel('DFT Values');
            legend('Buzzer', 'Tapping');
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).FFT_buzz = FFT_buzz;
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).FFT_tap = FFT_tap;
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).mag_buzz = mag_buzz;
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).mag_tap = mag_tap;
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).angle_buzz = angle_buzz;
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).angle_tap = angle_tap;
            subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).f_axis = f_axis;
        end
    end
end
close all;
%%
% Rearrange data by storing by rhythm/no_rhythm  
for rhCond = 1:size(subject_data(nSubj).rhythm_condition, 2)
    for trCond = 1:size(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition, 2)
        for nSubj = startSubj:size(subject_data, 2)
            if nSubj == 6 || nSubj == 5 || nSubj == 7 || nSubj ==4  
                continue;
            end
            mag_buzz_by_subj = subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).mag_buzz;
            mag_tap_by_subj = subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).mag_tap;
            angle_buzz_by_subj = subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).angle_buzz;
            angle_tap_by_subj = subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition(trCond).angle_tap; 
            
            rhythm_data(rhCond).condition(trCond).mag_buzz(:,nSubj) = mag_buzz_by_subj;
            rhythm_data(rhCond).condition(trCond).mag_tap(:,nSubj) = mag_tap_by_subj;
            rhythm_data(rhCond).condition(trCond).angle_buzz(:,nSubj) = angle_buzz_by_subj;
            rhythm_data(rhCond).condition(trCond).angle_tap(:,nSubj) = angle_tap_by_subj;         
            
        end
    end
end
%%
% Calculate bode plot values for each subject 
freq_indx = [4 6 8 12 18 24];
cond  = 1;
for rhCond = 1:size(subject_data(nSubj).rhythm_condition, 2)
    cond  = 1;
    for trCond = 1:size(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition, 2)
        for j = 1:length(nonzeros(condition_set(trCond,:)))
            indx = freq_indx(cond);
            mag_buzz_by_subj = rhythm_data(rhCond).condition(trCond).mag_buzz(indx,:);
            mag_tap_by_subj = rhythm_data(rhCond).condition(trCond).mag_tap(indx,:);
            angle_buzz_by_subj = rhythm_data(rhCond).condition(trCond).angle_buzz(indx,:);
            angle_tap_by_subj = rhythm_data(rhCond).condition(trCond).angle_tap(indx,:);

            mag_bode_by_subj = 20*log(mag_tap_by_subj./mag_buzz_by_subj);
            angle_bode_by_subj = (angle_tap_by_subj-angle_buzz_by_subj)/pi*180;
            angle_bode_by_subj(angle_bode_by_subj == 0) = NaN;
            rhythm_data(rhCond).freq_condition(cond).mag_bode_by_subj = mag_bode_by_subj;
            rhythm_data(rhCond).freq_condition(cond).angle_bode_by_subj = angle_bode_by_subj;
            cond = cond + 1;
        end
    end
end
%%
% Take average mag/angle for condition
for rhCond = 1:size(subject_data(nSubj).rhythm_condition, 2)
    cond  = 1;
    for trCond = 1:size(subject_data(nSubj).rhythm_condition(rhCond).trials_by_condition, 2)
        for j = 1:length(nonzeros(condition_set(trCond,:)))
            rhythm_data(rhCond).mag_bode_avg(cond) = nanmean(rhythm_data(rhCond).freq_condition(cond).mag_bode_by_subj);
            rhythm_data(rhCond).angle_bode_avg(cond) = nanmean(rhythm_data(rhCond).freq_condition(cond).angle_bode_by_subj);
            rhythm_data(rhCond).std_mag_bode_avg(cond) = nanstd(rhythm_data(rhCond).freq_condition(cond).mag_bode_by_subj);
            rhythm_data(rhCond).std_angle_bode_avg(cond) = nanstd(rhythm_data(rhCond).freq_condition(cond).angle_bode_by_subj);
            cond = cond + 1;
        end
    end
end

%% 
% Calculate bode plot values for each condition using avg of subjects
figure;
colors = ['r', 'b'];
for rhCond = 1:2
    hold on; 
    subplot(2,1,1);
    x = conditions(1:6)';
    y = rhythm_data(rhCond).mag_bode_avg';
    dy = rhythm_data(rhCond).std_mag_bode_avg';
    %errorbar(x,y,err)
    h1(rhCond) = fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(rhCond),'linestyle','none');
    line(x,y,'Color',colors(rhCond))
    set(h1(rhCond),  'facealpha',0.3)
    %lgd1(1) = semilogx(conditions(1:6), rhythm_data(rhCond).mag_bode_avg);
    
    hold on;
    subplot(2,1,2);
    x = conditions(1:6)';
    y = rhythm_data(rhCond).angle_bode_avg';
    dy = rhythm_data(rhCond).std_angle_bode_avg';
    %errorbar(x,y,err)
    h2(rhCond) = fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(rhCond),'linestyle','none');
    line(x,y, 'Color',colors(rhCond))
    set(h2(rhCond),  'facealpha',0.3)
%     lgd2(1) = semilogx(conditions(1:6), rhythm_data(rhCond).angle_bode_avg);
end
hold on;
subplot(2,1,1);
set(gca,'xscale','log')
title('Magnitude Plot');
xlabel('Frequency (cycles/event)');
ylabel('Magnitude (dB)');
axis([.05 .2875 -50 10]);

hold on;
subplot(2,1,2);
set(gca,'xscale','log')
title('Phase Plot');
xlabel('Frequency (cycles/event)');
ylabel('Phase (degrees)');
axis([.05 .2875 -100 10]);

legend(h1, 'No Rhythm', 'Rhythm', 'Location', 'southwest');
legend(h2, 'No Rhythm', 'Rhythm', 'Location', 'southwest');

