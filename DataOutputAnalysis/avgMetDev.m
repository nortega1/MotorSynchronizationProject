function [ avg_met_dev ] = avgMetDev( E, nz_buzz_times, nz_tap_times, plot_on )
%AVGMETDELAY Calculate the average metronome delay during a steady
%metronome
%Nicole Ortega @ 7/24/2017
met_period = E.met_period;
d = diff(diff(nz_buzz_times));
indx = find(d > .001 | d < -.001, 1, 'first'); 
steady_met = nz_buzz_times(1:indx);
steady_met = steady_met(end/2:end);
tap_dev = nan(length(steady_met),1);

for t = 1:length(steady_met)
    time = steady_met(t);
    h = time + met_period/2;
    l = time - met_period/2;
    
    taps = find(nz_tap_times < h & nz_tap_times > l);
    if length(taps) >= 2 || length(taps) <1
        continue;
    else
        tap_dev(t) = nz_tap_times(taps) - time;
    end
    
end

if plot_on
    figure;
    subplot(1,2,1);
    plot(tap_dev);
    title('Avg Steady Metronome Delay');
    xlabel('Metronome Count');
    ylabel('Delay (s)');
    subplot(1,2,2);
    histogram(tap_dev);
    title('Histogram of Avg Steady Metronome Delay');
    xlabel('Delay (s)');
    ylabel('Frequency');

end

avg_met_dev = nanmean(tap_dev);

end

