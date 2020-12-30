function isi(directory)
cd(directory);
datadir = [directory, '*data.mat'];
datafiles = dir(datadir);

%i got bored this isn't a real function

ISI = [];
for k=1:n_trials
    spike_times = find(d(k,:) == 1);
    isi0 = diff(spike_times);
    ISI = [ISI, isi0];
end
end