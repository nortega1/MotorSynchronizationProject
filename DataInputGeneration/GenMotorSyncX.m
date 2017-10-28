%GENMOTORSYNCX
%Third version of genTappingX.m
%Generates frequencies and their corresponding sinusoidal 
%perturbations for the tapping experiment. Uses genPerts_DT() to 
%generate decimal perturbations.
%
%There are five different types of trials, with two/three conditions: 
%       1. sum of sines: 0.0375, 0.0625, 0.0875 (no rhythm, rhythm)
%       2. single sine: 0.1375 (no rhythm, rhythm)
%       3. single sine: 0.2125 (no rhythm, rhythm)
%       4. single sine: 0.25   (no rhythm, rhythm, enhanced rhythm)
%       5. single sine: 0.2875 (no rhythm, rhythm)
%
%genPertsDT is used to generate perturbations
%
%NOTE: adds 10 pert buffer to the beginning and end of the trial
%NOTE: used to generate experiment 9.1.2017

close all; clear all;
%inputs by user
date = '9.1.2017';
subj_num = '1';
num_repeats = 3;    %number of repeats for each trial type with condition
% 5 different trial types. nan in front represented unsigned trial
% condition
trial_type = [nan 0 0 0 0 1 1 1; nan 0 0 0 1 0 0 0; nan 0 0 1 0 0 0 0; 
    nan 0 1 0 0 0 0 0; nan 1 0 0 0 0 0 0];
% 3 different conditions. row- trial type, column- condition. nan
% represents unused third condition for trial type 
conditions = [0 1 nan; 0 1 nan; 0 1 nan; 0 1 2; 0 1 nan];
total_cond = length(conditions(~isnan(conditions)));

%Make trials with conditions, then randomize
trials_all = NaN(total_cond*num_repeats, 8);
i = 1;
for t = 1:size(trial_type,1)
    type = trial_type(t, :);
    for c = 1:size(conditions,2)
        cond = conditions(t,c);
        for rep = 1:num_repeats
            if isnan(cond)
                break;
            end
            type(1) = cond;
            trials_all(i,:) = type;
            i = i +1;
        end
    end
end

%randomize trials
trials_all_random = trials_all(randperm(size(trials_all,1)),:);


%Open files to store data fed to microprocessor
bin_dev_freq = strcat(date,'/',date,'_Subject_',subj_num,'_BinDevFreqs.txt');
dec_dev_freq = strcat(date,'/',date,'_Subject_',subj_num,'_DecDevFreqs.txt');
rhy = strcat(date,'/',date,'_Subject_',subj_num,'_Rhythm.txt');


mkdir(date);
bin_dev_freq_file = fopen(bin_dev_freq,'w');   
dec_dev_freq_file = fopen(dec_dev_freq,'w'); 
rhythm_file = fopen(rhy,'w'); 

%Write on files
for i = 1:size(trials_all_random,1)
    bin_freq_sum = dec2bin(bin2dec(num2str(trials_all_random(i,2:8))),7);
    bin_freq_sum = strcat(num2str(trials_all_random(i,1)), bin_freq_sum);
    [dec_freq_sum, rhythm, dec_sum_perts_no_buff] = genPertsDT(bin_freq_sum);
    close all;
    buffer = dec_sum_perts_no_buff(1:10);
    dec_sum_perts_w_buffer = [buffer dec_sum_perts_no_buff buffer];
    fprintf(bin_dev_freq_file,'%s\n', bin_freq_sum);
    for l = 1:length(dec_freq_sum)
        fprintf(dec_dev_freq_file,'%1.4f ', dec_freq_sum(l));
    end
    fprintf(dec_dev_freq_file,'\n', dec_freq_sum(l));
    fprintf(rhythm_file, '%d\n', rhythm);
    
    dec_perts = strcat(date,'/',date,'_Subject_',subj_num,'.',num2str(i),'_DecPerts.txt');
    dec_perts_file = fopen(dec_perts,'w'); 
    fprintf(dec_perts_file,'%1.3f\n', dec_sum_perts_w_buffer);    
    fclose(dec_perts_file);
end

fclose(bin_dev_freq_file);
fclose(dec_dev_freq_file);
fclose(rhythm_file);
close all;
