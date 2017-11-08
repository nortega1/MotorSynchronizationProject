function [ logFile1 ] = fixTrial( logFile )
%FIXTRIAL takes in a log of data that has been corrupted and outputs an
%edited and uncorrupted version

logFile1 = logFile(590:end,:);
time = diff(logFile1(:,3));
indx = find(time < 0);

while ~isempty(indx)
    logFile1(indx, :) = [];
    time = diff(logFile1(:,3));
    indx = find(time < 0);
end

end

