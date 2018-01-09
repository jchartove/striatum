% http://matlaboratory.blogspot.co.uk/2015/05/event-spike-triggered-average.html

function [avg, nEvs] = evTrigAvg(events, stim, windowSize)

% Discard spikes in first window
events(1:windowSize)=0;

% Find events:
% Number of events
nEvs = sum(events);
% Indexes of events
evIdx = find(events);

% Preallocate average
avg = zeros(1,windowSize);
% For each event
for w = 1:nEvs
    % Find the indexes of the time window preceding the event
    wIdx = evIdx(w)-windowSize : evIdx(w)-1;
    % Add the stim from this window to the average
    avg = avg + stim(wIdx);
end
% Divide by number of events to complete average
avg=avg./sum(sum(events)); %edited to accomodate 2d data - jc
