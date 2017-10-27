% Decodes data files and fixes trials with corruptions

%% Subject information 
% Edit as needed 
date = '9.1.2017';
subject_num = '6';
num_trials = '33';

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
corr_files = NaN(str2num(num_trials),1);

while(curr_trial <= str2double(num_trials))
    log_file = csvread(strcat('/',date, '/', date,'_Subject_',subject_num,'.',num2str(curr_trial),'_Res.txt'));
    log_file = log_file(1:end-2,:); %cutting out extra data points at the end of file
    
    %check for corruption
    l_file = fixTrial(log_file);
    time = log_file(600:end-2,3);
    d_time = diff(time);
    disorder = d_time(d_time < 0); %check if time chronological
    fixTrial(log_file);
    if length(disorder) > 400 || ~isempty(find(isnan(time)))
        %record corruption
        corr_files(curr_trial) = curr_trial;
        input = sprintf(format_spec, path_trial_decode, date, subject_num, num2str(curr_trial), num2str(fix));
        system(input);
        fix = fix + 1;
        if fix  == 8
            fix = 1;
            curr_trial = curr_trial + 1;
        end
        display(curr_trial);
    else
        curr_trial = curr_trial + 1;
        fix = 1;
    end
    
        
end

%record corrupted files
fileID = fopen(fullfile(date, strcat(date,'_Subject_',subject_num,'_Corr.txt')), 'w');
corr_files(isnan(corr_files)) = [];
fprintf(fileID,'%d ',corr_files);
fclose(fileID);



