% Nicole Ortega (c) 10/26/2017
% Decodes data files and fixes trials with corruptions

%% Subject information 
% Edit as needed 
date = '9.1.2017';
subject_num = '1';
num_trials = '33';

%% Decode All Subject Files
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
corr_files = NaN(str2num(num_trials),1);

while(curr_trial <= str2double(num_trials))
    log_file = csvread(strcat('/',date, '/', date,'_Subject_',subject_num,'.',num2str(curr_trial),'_Res.txt'));
    log_file = log_file(1:end-2,:); %cutting out extra data points at the end of file
    
    time = log_file(600:end,3);
    d_time = diff(time);
    disorder = d_time(d_time < 0); %chronological in time
    
    %check for corruption and attempt to fix it
    if length(disorder) > 400 || ~isempty(find(isnan(time)))  
        corr_files(curr_trial) = curr_trial;
        input = sprintf(format_spec, path_trial_decode, date, subject_num, num2str(curr_trial), num2str(fix));
        system(input);
        if fix  == 8    %unable to fix corruption
            fix = 1;
            fprintf('Unable to fix corrupted trial: %d\n', curr_trial);
            curr_trial = curr_trial + 1;
        end
        fix = fix + 1;
        fprintf('Fixing corrupted trial: %d\n', curr_trial);
    else
        curr_trial = curr_trial + 1;
        fix = 1;
    end
        
end

%% Record Corrupted Files
fileID = fopen(fullfile(date, strcat(date,'_Subject_',subject_num,'_Corr.txt')), 'w');
corr_files(isnan(corr_files)) = [];
if isempty(corr_files) 
    corr_files(1) = 0;
end
display(corr_files');
fprintf(fileID,'%d ',corr_files);
fclose(fileID);



