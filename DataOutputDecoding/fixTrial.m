function [ logFile1 ] = fixTrial( logFile )
% FIXTRIAL takes in a corrupted file and removes corrupted areas. Oupts
% fixed trial
% Nicole Ortega (c) 10/2017
logFile1 = logFile(590:end,:);
time = diff(logFile1(:,3));
indx = find(time < 0);

while ~isempty(indx)
    logFile1(indx, :) = [];
    time = diff(logFile1(:,3));
    indx = find(time < 0);
end

end

