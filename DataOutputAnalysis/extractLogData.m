function [ nz_met_times, nz_buzz_times, nz_tap_times] = ExtractLogData_rev0( E, logFile, dev_freq, plot_on)
%STORELOGDATA store raw data from logFile. Plot encoder counts vs. time
%Nicole Ortega @ 7/24/2017

    %% ---- Extract Data ---- %%
    met_period = E.met_period;
    thresh = E.thresh_val;
    
    buzz = logFile(1:end-2,1);                            % store
    time = logFile(1:end-2,3);
    encoder = logFile(1:end-2,4);

    %% --- Plotting Encoder Readings --- %%
    % Plot time on the x-axis and encoder on the y-axis 
    if plot_on
        figure();
        plot(time(5:end), encoder(5:end));
        str= sprintf('Encoder Counts vs. Time for a deviation of %.2f %.2f %.2f Hz', dev_freq);
        title(str);
        xlabel('Time (s)') % x-axis label 
        ylabel('Encoder Counts') % y-axis label
        hold on;
        hline = refline([0 thresh]);
        hline.Color = 'r';
    end
    
    %% --- Experimental Buzzer --- % 
    % During the experiment, the buzzer first follows a steady metronome,
    % then a perturbed metronome
    
    [bnum_rows, bnum_cols] = size(buzz);
    
    % Cleaning up buzzer data 
    buzz  = diff(buzz);                             % mark start of buzzer with a 1
    buzz(buzz == -1) = 0;                           % delete end of buzzer (-1)
    buzz(end+1) = 0;                                % add value
    last_buzz = find(buzz, 1, 'last');              % find last 1 and delete
    buzz(last_buzz) = 0;

    buzz_times = zeros(bnum_rows,1);

    % Record times when buzzer sounds
    % 1 - beginnning of buzzer, 0 - buzzer is off 
    % Experimental
    for i = 1:bnum_rows
        if buzz(i) == 1
            buzz_times(i) = time(i,:);             
        else
            buzz_times(i) = 0;
        end 
    end
    nz_buzz_times = nonzeros(buzz_times);           % nz = nonzero

    %% --- Create a Steady Metronome --- %%
    % Create a steady ~3 Hz (period of .333) metronome 
    
    nz_met_times = nz_buzz_times(1):met_period:nz_buzz_times(end)+met_period;
    nz_met_times = nz_met_times';
    
    [r1, c1] = size(nz_met_times);
    [r2, c2] = size(nz_buzz_times);
    
    while (r1 ~= r2) 
        nz_met_times = nz_met_times(1:(r1 - 1),1);
        [r1, c1] = size(nz_met_times);
    end

    %% ---- Tapping Experiment Times ---- %%
    % Using the threshold to determine when a subject has "tapped"
    
    encoder_post = encoder;
    rows = length(encoder_post);
    for i = 1:rows
        if encoder_post(i) < thresh 
            encoder_post(i) = 1;
        else 
            encoder_post(i) = 0;
        end
    end
    encoder_post = diff(encoder_post);
    encoder_post(encoder_post == -1) = 0;
    encoder_post(end+1) = 0;
    tap = encoder_post;

    tap_times = zeros(bnum_rows,1);

    % Record times when tapping occurs
    % 1 - beginnning of tapping, 0 - no tapping
    % Experimental
    for i = 1:bnum_rows
        if tap(i) == 1
            tap_times(i) = time(i,:); 
        else
            tap_times(i) = 0;
        end
    end

    % Extract nonzero tap times
    nz_tap_times = nonzeros(tap_times);

end

