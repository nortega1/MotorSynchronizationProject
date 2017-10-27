function [ logFile1 ] = fixTrial( logFile )
%FIXTRIAL Summary of this function goes here
%   Detailed explanation goes here
logFile1 = logFile(590:end,:);
time = diff(logFile1(:,3));
indx = find(time < 0);

while ~isempty(indx)
    logFile1(indx, :) = [];
    time = diff(logFile1(:,3));
    indx = find(time < 0);
end

end

