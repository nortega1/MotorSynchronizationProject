% Decodes data files and fixes trials with corruptions

%% Subject information 
% Edit as needed 
date = '7.12.2017';
subject_num = '1';
num_trials = '18';

%% Decode all subject files
% Edit number of trials of subject
path_subj_decode = '"/Users/Nicole/Documents/JHU/LIMBS/MotorSynchronizationProject/LogDecoder/LogDecoderSubject_MAC"';
format_spec = '%s %s %s %s';
input = sprintf(format_spec, path_subj_decode, date, subject_num, num_trials);
system(input);

%% Fix Corruption
% Check and fix corruptions
curr_trial = 1; %current trial
fix = 1;
path_trial_decode = '"/Users/Nicole/Documents/JHU/LIMBS/MotorSynchronizationProject/LogDecoder/LogDecoderTrial_MAC"';
format_spec = '%s %s %s %s %s';

while(curr_trial <= str2double(num_trials))
    log_file = csvread(strcat('/',date, '/', date,'_Subject_',subject_num,'.',num2str(curr_trial),'_Res.txt'));
    log_file = log_file(1:end-2,:); %cutting out extra data points at the end of file
    
    %check for corruption
    time = log_file(600:end-2,3);
    d_time = diff(time);
    disorder = d_time(d_time < 0); %check if time chronological
    if ~isempty(disorder)
        input = sprintf(format_spec, path_trial_decode, date, subject_num, num2str(curr_trial), num2str(fix));
        system(input);
        fix = fix + 1;
        if fix  == 8
            fix = 1;
            curr_trial = curr_trial + 1;
        end
    else
        curr_trial = curr_trial + 1;
        fix = 1;
    end
    
        
end



