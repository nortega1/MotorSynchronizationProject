function [ freq_trials ] = readBinDevFreq( subj_num, date, num_trials, dev_harmonics)
%READDEVFREQ reads the input txt file of a subject's deviation
%frequencies (represented as binary numbers).
%Outputs a structure with all the conditions (frequencies
%tested) on the subject, in the order that they were tested in. 
%
%This function is used for single or sum of sines testing.
%INPUT: subject number, date and total number of trials to open 
%text file containing a vertical list of dev frequencies (in a binary 
%representation) tested in the order they appeared to the participant.
%OUTPUT: dev frequencies tested (numerical representation)

file = strcat(strcat(date, '/', date,'_Subject_',num2str(subj_num),'_BinDevFreqs.txt'));
bin_freq_n = dlmread(file);
bin_freq_s = num2str(bin_freq_n,'%08i');
bin_freq_s = bin_freq_s';

fundamental_period = 80;   %events per cycle, AKA perturbations per period
fundamental_freq = 1/fundamental_period;    %cycles per event
dev_freq_dec = fundamental_freq * dev_harmonics;
freq_trials = zeros(num_trials, length(dev_harmonics)+1);
trial = 1;
for i = 1:8:8*num_trials
    bin_freq_trial = bin_freq_s(i:i+7);
    for j = 2:8
        if str2num(bin_freq_trial(j))
            freq_trials(trial, j) = dev_freq_dec(9-j);
        end
    end
    trial = trial+1;
end

end

