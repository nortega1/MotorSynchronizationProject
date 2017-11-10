function [Y, YY, P1, P3, A1, A3, f] = calcFFT(E, sum_buzz_dev, sum_tap_dev, plot_on)
%CALCFFT Calculates the fft of the deviation from the metronome for both
%the buzzing and tapping data. 
%
%NOTE: uses fundamental period and number of repeats to determine the
%length of the trial. AKA number of taps
%NOTE: subtracts the bias DC component to eliminate low frequency noise
%% ---------- Calculating the FFT ---------- %%
%display sinusoidal sum perturbation to metronome
L =  E.fundamental_period;  %length of trial, number of taps
N = (1:1:L);     %number of taps
sum_buzz_dev = sum_buzz_dev - mean(sum_buzz_dev); %subtracting the bias DC component
sum_tap_dev = sum_tap_dev - mean(sum_tap_dev);

if plot_on
    figure();
    plot(N, sum_buzz_dev, '*-b');
    hold on; 
    plot(N, sum_tap_dev, '-.r');
    xlabel('Count (tap)');
    ylabel('Deviation from the Metronome (s)');
    legend('Theoretical', 'Experimental');
end

Y = fft(sum_buzz_dev);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = (0:(L/2))/L;

A1 = angle(Y)+pi/2;
A1 = A1(1:L/2+1);

YY = fft(sum_tap_dev);
P4 = abs(YY/L);
P3 = P4(1:L/2+1);
P3(2:end-1) = 2*P3(2:end-1);
f = (0:(L/2))/L;
A3 = angle(YY)+pi/2;
A3 = A3(1:L/2+1);

if plot_on
    figure();
    plot(f, P1,'-*b'); hold on;
    plot(f, P3,'-.r');
    figure();
    plot(f, A1,'-*b'); hold on;
    plot(f, A3,'-.r');
end

end
