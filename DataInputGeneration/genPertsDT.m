function [dev_freqs, rhythm, sum_perts] = genPertsDT(bin_freq)
%GENPERTSDT takes in binary representation of the deviation frequencies being
%tested and outputs the decimal representation of each of the deviation 
%frequency/ies, whether or not a rhythm was present and the sinusodial 
%perturbations to the metronome


%% ------ Constants ------ %%
fundamental_period = 80;   %events per cycle, AKA perturbations per period
fundamental_freq = 1/fundamental_period;    %cycles per event
repeats = 2;            %number of times a complete trial will repeat

% these prime multiples are approx log spaced. 
dev_harmonics = [3 5 7 11 17 20 23]; % 23 29 ** 37 max for period = 80
% phase shifts were randomly chosen using rand(7,1)*2*pi
dev_phase_shift = [1.387 0.257 3.324 2.512 1.633 1.141 2.710];
dev_amplitude = [0.029 0.036 0.043 0.036 0.045 0.055 0.055];
%%
dev_freq_dec = fundamental_freq*dev_harmonics;

num_dev_freq = 7;
count = 1;
for i = 2:length(bin_freq)
    if str2num(bin_freq(i))
        dev_freqs(count) = dev_freq_dec(num_dev_freq);
        dev_amps(count) = dev_amplitude(num_dev_freq);
        dev_phase(count) = dev_phase_shift(num_dev_freq);
        dev_harmonics(count) = dev_harmonics(num_dev_freq); 
        count = count + 1;
    end
    num_dev_freq = num_dev_freq-1;
end

if str2num(bin_freq(1)) == 1
    rhythm = 1;
elseif str2num(bin_freq(1)) == 2
    rhythm = 2;
else
    rhythm = 0;
end


perts = zeros(length(dev_freqs),fundamental_period*repeats);
for j=1:length(dev_freqs)
    count = 1;
    for n=0:1:fundamental_period*repeats-1
        %the amplitude must be divided by sqrt(len) ONLY when dev_amps differ
        perts(j,count) = round(dev_amps(j)*sin(2*pi*dev_freqs(j)*n + dev_phase(j)), 4);
        count = count+1;
    end
end
sum_perts = sum(perts,1);

%% ---------- Calculating the FFT ---------- %%
%display sinusoidal sum perturbation to metronome
figure();
t = (1:1:fundamental_period);
plot(t, sum_perts(1:80),'*-');
str = sprintf('Deviation from the Steady Metronome: %.3f %.3f %.3f', dev_freqs);
title(str);
xlabel('Events (tap)');
ylabel('Deviation from the Metronome (s)');

%
repeats = 1;
L = fundamental_period*repeats;       %events per cycle
Y = fft(sum_perts(1:80));
figure;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = (0:(L/2))/L;
hold on;
plot(f, P1,'-*');
title('Magnitude Plot');
xlabel('Frequency (cycles per event)');
ylabel('Magnitude');
[M, ind] = max(P1);
display(P1(ind));
A3 = angle(Y);
A3 = A3(1:L/2+1)+pi/2;
display(A3(ind));
figure;
plot(f,A3, '--');
title('Phase Plot');
xlabel('Frequency (cycles per event)');
ylabel('Phase (rads)');





end