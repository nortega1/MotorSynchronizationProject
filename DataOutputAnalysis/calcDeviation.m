function [ buzz_dev, tap_dev, percentage_taps ] = calcDeviation( nz_met_times, nz_buzz_times, nz_tap_times, avg_met_dev, plot_on )
%CALCDEVIATION Calculate the deviation of buzzer and tapping from a steady
%metronome 
%Nicole Ortega @ 7/24/2017

    thrown_taps = 0; 
    %% --- Deviation of buzzer from steady metronome --- %%
    
    buzz_dev = round(nz_buzz_times,3) - round(nz_met_times,3);
    
    %% --- Plot time vs. buzzer deviation from steady metronome --- %%
    if plot_on
        figure();
        plot(nz_met_times, buzz_dev, '*-');
        title('Full Trial: Deviation from Steady Metronome vs. Time');
        xlabel('Time(s)');
        ylabel('Deviation from Steady Metronome (s)');
        hold on; 
    end
    
    %% --- Deviation of tapping from steady metronome --- %%
    tap_dev = nan(length(nz_buzz_times),1);

    for t = 1:length(nz_buzz_times)
        time = nz_met_times(t);
        h = time + .15 + avg_met_dev;
        l = time - .15 + avg_met_dev;

        taps = find(nz_tap_times < h & nz_tap_times > l);
        if length(taps) >= 2 || length(taps) < 1
            if(length(taps) >= 2)
                thrown_taps = thrown_taps + length(taps);
            end
            continue;
        else
            tap_dev(t) = nz_tap_times(taps) - time;
        end
    end
    if plot_on
        plot(nz_met_times, tap_dev, '*-');
    end
    
    %percentage of taps discarded over length of total taps
    percentage_taps = thrown_taps/length(nz_tap_times);

end

