% Plots data from subject. Can help visually inspect if corruptions have been
% resolved 
date = '9.1.2017';
subject_num = '1';
num_trials = '33';

corr_files = csvread(strcat('/',date, '/', date,'_Subject_',subject_num,'_Corr.txt'));
for curr_trial = 1:str2double(num_trials)
    log_file = csvread(strcat('/',date, '/', date,'_Subject_',subject_num,'.',num2str(curr_trial),'_Res.txt'));
    % corrupted files
    if find(curr_trial == corr_files)
        log_file = fixTrial(log_file);
    end
    
    % plot data 
    figure(); plot(log_file(600:end-2,3),log_file(600:end-2,4));
    xlabel('Time'); ylabel('Encoder Counts');
end
