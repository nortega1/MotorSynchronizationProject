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

%% ------ Clear ------ %%
clear all; 
close all;

%% ------ Variables ------ %%
% Edit as needed depending on experiment 
% Variable E holds experimental variables/constants
E = struct;                         
E.start_subj = 1;
E.total_subj = 1;                   % total number of subjects
E.date = '9.1.2017';                % E.date for experiment (type of5experiment)
E.trials_per_subj = 33;

%% ------ Constants ------ %%
E.sample_freq  = 1000;                % 1k Hz
E.met_period = .333;
E.met_freq = 1/E.met_period;            % freq of steady metronome in Hz
E.fundamental_period = 80;            % events per cycle, AKA perturbations per period
E.fundamental_freq = 1/E.fundamental_period;    %cycles per event
E.repeat = 2;                        %number of times a complete trial will repeat
E.thresh_val = -105;                  % Estimating threshold value
E.rest =  5;                          % time in seconds before start of trial


devfreq_info = dlmread(strcat(E.date, '/', E.date,'_DevFreq_Info.txt'));
E.dev_harmonics = devfreq_info(1,:);
E.dev_phase_shift = devfreq_info(2,:);
E.dev_amplitude = devfreq_info(3,:);

% Tested frequencies
E.conditions = E.dev_harmonics * E.fundamental_freq;
E.condition_set = [E.conditions(1:3); E.conditions(4) 0 0; ...
    E.conditions(5) 0 0; E.conditions(6) 0 0; E.conditions(7) 0 0];

E.signal_len = 260;     %length of trial

%% ------ Initializing Variables ------ %%
reps = zeros(size(E.condition_set,1), 1);
reps_old = zeros(size(E.condition_set,1), 1);
reps_r = zeros(size(E.condition_set,1), 1);
reps_nr = zeros(size(E.condition_set,1), 1);
perc_tapsvsbuzz = nan(size(E.condition_set,1), E.total_subj);
perc_tapsthrown = nan(size(E.condition_set,1), E.total_subj);

subject_data = struct;
rhythm_data = struct;

%% ---- Open Files and Store Data ---- % 
% Store raw data (metronome times, buzzer times, tapping times) and find
% average deviation from steady metronome. Data will later be used to
% determine the buzzer/tapping deviation from steady metronome

for curr_subj = E.start_subj:E.total_subj
    
    if  curr_subj == 7 
        continue;
    end
    curr = struct;
    
    % Rhythms tested for subject in the order they appear
    % Total number of trials for subject 
    % Files with fixed and unfixed corruptions
    rhythm_file = dlmread(strcat(E.date, '/', E.date,'_Subject_',num2str(curr_subj),'_Rhythm.txt'));
    fixed_trials = dlmread(strcat(E.date, '/', E.date,'_Subject_',num2str(curr_subj),'_Corr.txt'));
    unfixed_trials = dlmread(strcat(E.date, '/', E.date,'_Subject_',num2str(curr_subj),'_Skip.txt'));

    % Frequencies tested for subject in the order they appear
    freq_devs_by_subject = readBinDevFreq(curr_subj, E);
    
    avg_met_dev = nan(E.trials_per_subj, 1);
    rhythm_reps = zeros(length(E.condition_set), E.trials_per_subj, 2);
    
    plot_on = 0;    %produce plots from this section
    for curr_trial = 1:E.trials_per_subj
        if find(curr_trial == unfixed_trials)
            continue;
        end
        % initialize variables
        plot_on = 0;    %produce plots from this section
        n_met_all = nan(E.signal_len,1); n_buzz_all = nan(E.signal_len,1); n_tap_all = nan(E.signal_len,1);
        buzz_dev_all = nan(E.signal_len,1); tap_dev_all = nan(E.signal_len,1);
        
        % extract deviation freq being tested during trial
        freq_devs_by_trial = nonzeros(freq_devs_by_subject(curr_trial,:));
        % determine condition being tested by comparing trial condition to
        % condition set
        [r,c] = find(freq_devs_by_trial(1) == E.condition_set);
        n_cond = r;
       
        % read input and output file
        log_file = dlmread(strcat(E.date, '/', E.date,'_Subject_',num2str(curr_subj),'.',num2str(curr_trial),'_Res.txt'));
        input_file = dlmread(strcat(E.date, '/', E.date,'_Subject_',num2str(curr_subj),'.',num2str(curr_trial),'_DecPerts.txt'));
        % skip unsuccessful trials
        if find(curr_trial == fixed_trials)
            log_file = fixTrial(log_file);
        end
        log_file = log_file(1:end-2,:); %cutting out extra data points at the end of file
        
        % Extract steady metronome times, buzzer times, and tapping times
        % Store in vectors of the same length
        [ n_met, n_buzz, n_tap ] = extractLogData(E, log_file, freq_devs_by_trial, plot_on);
        n_met_all(1:length(n_met)) = n_met;  n_buzz_all(1:length(n_buzz)) = n_buzz; n_tap_all(1:length(n_tap)) = n_tap;
        
        % Calculate the average delay during a steady metronome 
        [ subject_data(curr_subj).avg_met_dev(curr_trial) ] = avgMetDev(E, n_buzz, n_tap, plot_on);

        rhythm_cond = rhythm_file(curr_trial) + 1;
        rhythm_reps(n_cond, curr_trial, rhythm_cond) = 1;
        rep = sum(rhythm_reps(n_cond, :, rhythm_cond));
        
        if rhythm_cond == 3
            n_cond = 1;
        end
            
        subject_data(curr_subj).rhythm_condition(rhythm_cond).trials_by_condition(n_cond).n_met(:,rep) = n_met_all;
        subject_data(curr_subj).rhythm_condition(rhythm_cond).trials_by_condition(n_cond).n_buzz(:,rep) = n_buzz_all;
        subject_data(curr_subj).rhythm_condition(rhythm_cond).trials_by_condition(n_cond).n_tap(:,rep) = n_tap_all;
        
        close all;
        
    end
    
end
close all;
%% 

for curr_subj = E.start_subj:size(subject_data, 2)
    rhythm_file = dlmread(strcat(E.date, '/', E.date,'_Subject_',num2str(curr_subj),'_Rhythm.txt'));
    num_trials = length(rhythm_file);
    % Frequencies tested per subject in the order they appear
    freq_devs_by_subject = readBinDevFreq(curr_subj, E.date, num_trials, dev_harmonics);
    
    for rhCond = 1:size(subject_data(curr_subj).rhythm_condition, 2)
        
        for trCond = 1:size(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition, 2)
            
            for rep = 1: size(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met, 2)
                %initializing variables
                buzz_dev_all = nan(signal_len,1); tap_dev_all = nan(signal_len,1);
                
                % extract deviation freq being tested during trial
                freq_devs_by_trial = nonzeros(freq_devs_by_subject(curr_trial,:));
                avg_met_dev = mean(nonzeros(subject_data(curr_subj).avg_met_dev(:)));

                n_met = subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met(:,rep);
                n_buzz = subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_buzz(:,rep);
                n_tap = subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_tap(:,rep);

                % take out nan values
                n_met=(n_met(~isnan(n_met)));
                n_buzz=(n_buzz(~isnan(n_buzz)));
                n_tap=(n_tap(~isnan(n_tap)));


                [ buzz_dev, tap_dev, percent_taps ] = calcDeviation(n_met, n_buzz, n_tap, avg_met_dev);
                buzz_dev_all(1:length(n_met)) = buzz_dev; tap_dev_all(1:length(n_met)) = tap_dev;
                subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev(:,rep) = buzz_dev_all;
                subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev(:,rep) = tap_dev_all;
                
                close all;
                %statistic of tapping/buzzer percentage (raw data)
                p_tb(1:2) = [perc_tapsvsbuzz(trCond, curr_subj)  length(n_tap)/length(n_buzz)];
                perc_tapsvsbuzz(trCond, curr_subj) = nanmean(p_tb);

                %statistics of percentage taps kept to find deviation
                p_tt(1:2) = [perc_tapsthrown(trCond, curr_subj) percent_taps]; 
                perc_tapsthrown(trCond, curr_subj) = nanmean(p_tt);

            end
        end
    end
end

close all;
%%
% Calculate nanmean per subject 
for curr_subj = E.start_subj:size(subject_data, 2)
    for rhCond = 1:size(subject_data(curr_subj).rhythm_condition, 2)
        for trCond = 1:size(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition, 2)
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met_avg = nanmean(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met,2);
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_buzz_avg = nanmean(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_buzz,2);
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_tap_avg = nanmean(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_tap,2);
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev_avg = nanmean(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev,2);
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev_avg = nanmean(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev,2);
        end
    end

end

for curr_subj = E.start_subj:size(subject_data, 2)
    for rhCond = 1:size(subject_data(curr_subj).rhythm_condition, 2)
        for trCond = 1:size(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition, 2)
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met_avg = nanmean([subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met_avg(42:121) subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_met_avg(122:201)],2);
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_buzz_avg = nanmean([subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_buzz_avg(42:121) subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_buzz_avg(122:201)],2);
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_tap_avg = nanmean([subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_tap_avg(42:121) subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).n_tap_avg(122:201)],2);
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev_avg = nanmean([subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev_avg(42:121) subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev_avg(122:201)],2);
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev_avg = nanmean([subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev_avg(42:121) subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev_avg(122:201)],2);
        end
    end
end
%%
% Calculate FFT for individual subjects 
for curr_subj = E.start_subj:size(subject_data, 2)
    for rhCond = 1:size(subject_data(curr_subj).rhythm_condition, 2)
        for trCond = 1:size(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition, 2)
            buzz_dev = subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).buzz_dev_avg;
            tap_dev = subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).tap_dev_avg;
            [FFT_buzz, FFT_tap, mag_buzz, mag_tap, angle_buzz, angle_tap, f_axis] = calcFFT(fundamental_period, 1, buzz_dev, tap_dev);
            str = sprintf('FFT of Frequency Deviations: %.2f %.2f %.2f', freq_devs_by_trial);
            title(str);
            xlabel('Frequency (cycles per event)');
            ylabel('DFT Values');
            legend('Buzzer', 'Tapping');
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).FFT_buzz = FFT_buzz;
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).FFT_tap = FFT_tap;
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).mag_buzz = mag_buzz;
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).mag_tap = mag_tap;
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).angle_buzz = angle_buzz;
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).angle_tap = angle_tap;
            subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).f_axis = f_axis;
        end
    end
end
close all;
%%
% Rearrange data
for rhCond = 1:size(subject_data(curr_subj).rhythm_condition, 2)
    for trCond = 1:size(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition, 2)
        for curr_subj = E.start_subj:size(subject_data, 2)
            if curr_subj == 7 
                continue;
            end
            mag_buzz_by_subj = subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).mag_buzz;
            mag_tap_by_subj = subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).mag_tap;
            angle_buzz_by_subj = subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).angle_buzz;
            angle_tap_by_subj = subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition(trCond).angle_tap; 
            
            rhythm_data(rhCond).condition(trCond).mag_buzz(:,curr_subj) = mag_buzz_by_subj;
            rhythm_data(rhCond).condition(trCond).mag_tap(:,curr_subj) = mag_tap_by_subj;
            rhythm_data(rhCond).condition(trCond).angle_buzz(:,curr_subj) = angle_buzz_by_subj;
            rhythm_data(rhCond).condition(trCond).angle_tap(:,curr_subj) = angle_tap_by_subj;         
            
        end
    end
end
%%
% Calculate bode plot values for each subject 
freq_indx = [4 6 8 12 18 21 24];
cond  = 1;
for rhCond = 1:size(subject_data(curr_subj).rhythm_condition, 2)
    cond  = 1;
    for trCond = 1:size(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition, 2)
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
for rhCond = 1:size(subject_data(curr_subj).rhythm_condition, 2)
    cond  = 1;
    for trCond = 1:size(subject_data(curr_subj).rhythm_condition(rhCond).trials_by_condition, 2)
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
x = struct;
y = struct;
for rhCond = 1:2
    hold on; 
    subplot(2,1,1);
    x.mag.rhy(:,rhCond) = conditions';
    y.mag.rhy(:,rhCond) = rhythm_data(rhCond).mag_bode_avg';
    dy = rhythm_data(rhCond).std_mag_bode_avg';
    %errorbar(x,y,err)
    h1(rhCond) = fill([x.mag.rhy(:,rhCond);flipud(x.mag.rhy(:,rhCond))],[y.mag.rhy(:,rhCond)-dy;flipud(y.mag.rhy(:,rhCond)+dy)],colors(rhCond),'linestyle','none');
    line(x.mag.rhy(:,rhCond),y.mag.rhy(:,rhCond),'Color',colors(rhCond))
    set(h1(rhCond),  'facealpha',0.3)
    %lgd1(1) = semilogx(conditions, rhythm_data(rhCond).mag_bode_avg);
    
    hold on;
    subplot(2,1,2);
    x.phase.rhy(:,rhCond) = conditions';
    y.phase.rhy(:,rhCond) = rhythm_data(rhCond).angle_bode_avg';
    dy = rhythm_data(rhCond).std_angle_bode_avg';
    %errorbar(x,y,err)
    h2(rhCond) = fill([x.phase.rhy(:,rhCond);flipud(x.phase.rhy(:,rhCond))],[y.phase.rhy(:,rhCond)-dy;flipud(y.phase.rhy(:,rhCond)+dy)],colors(rhCond),'linestyle','none');
    line(x.phase.rhy(:,rhCond), y.phase.rhy(:,rhCond), 'Color',colors(rhCond))
    set(h2(rhCond),  'facealpha',0.3)
    %lgd2(1) = semilogx(conditions, rhythm_data(rhCond).angle_bode_avg);
end
hold on;
subplot(2,1,1);
set(gca,'xscale','log')
title('Magnitude Plot');
xlabel('Frequency (cycles/event)');
ylabel('Magnitude (dB)');
axis([.05 .3 -60 20]);

hold on;
subplot(2,1,2);
set(gca,'xscale','log')
title('Phase Plot');
xlabel('Frequency (cycles/event)');
ylabel('Phase (degrees)');
axis([.05 .3 -200 50]);

legend(h1, 'No Rhythm', 'Rhythm', 'Location', 'southwest');
legend(h2, 'No Rhythm', 'Rhythm', 'Location', 'southwest');

%% for single subjects
% no error included, no .25 magnitude
figure;
x.mag.rhy(6,:) = [];
x.phase.rhy(6,:) = [];
y.phase.rhy(6,:) = [];
y.mag.rhy(6,:) = [];
subplot(2,1,1);
plot(x.mag.rhy, y.mag.rhy);
title('Magnitude Plot');
xlabel('Frequency (cycles/event)');
ylabel('Magnitude (dB)');
set(gca,'xscale','log')
axis([.05 .3 -60 20]);
legend('No Rhythm', 'Rhythm', 'Location', 'southwest');

subplot(2,1,2);
plot(x.phase.rhy, y.phase.rhy);
title('Phase Plot');
xlabel('Frequency (cycles/event)');
ylabel('Phase (degrees)');
set(gca,'xscale','log')
axis([.05 .3 -200 50]);
legend('No Rhythm', 'Rhythm', 'Location', 'southwest');





