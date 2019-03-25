function gvObj = gvRunPower

% make gv object
gvObj = gv();

% import ds data
gvObj.importDsData;

axismeta = gvObj.model.data.dsData.axis(1).axismeta;

%% power freq axis
% add new axes
gvObj.model.data.dsData.axis(end+1) = gvArrayAxis;

% add axis info
gvObj.model.data.dsData.axis(end).name = 'freq';
gvObj.model.data.dsData.axis(end).values = 1;

%% power data
% import power data
powerResults = dsImportResults(fullfile(pwd, 'results'), 'import_scope','custom', 'func','get_fft', 'as_cell',1);
nFreqs = length(powerResults{1});
freqs = 1:nFreqs;

%% sim IDS
simIDs = gvObj.model.data.dsData.data(end,:,:,:,:,:,:,:,:); % colons are for more vary dims HACK
simIDs = cell2mat(simIDs);

simIDsSize = size(simIDs);
if length(simIDsSize) < dsDataNdims
  simIDsSize(end+1:dsDataNdims) = dsDataSize(length(simIDsSize)+1:end);
end

simIDsSize(end+1) = nFreqs; % for later reshape

simIDs = simIDs(:); % reshape to col vector

%% get axis names and vals
axis_names = gvObj.model.data.dsData.exportAxisNames;
axis_vals = gvObj.model.data.dsData.exportAxisVals;
axis_vals{1} = {'power'};
axis_vals{end} = freqs;

%% make new data of right dims
data = zeros(length(simIDs), nFreqs); % col vector

nSims = length(powerResults);
for simID = 1:nSims
  if ~isempty(powerResults{simID})
    data(simIDs == simID, :) = powerResults{simID};
  else
    data(simIDs == simID, :) = nan;
  end
end

% reshape the data
data = reshape(data, simIDsSize);

%% add data
tempObj = gvArray;
tempObj = tempObj.importData(data, axis_vals, axis_names);

%% merge
gvObj.model.data.dsData.merge(tempObj);

%% add axis metadata
axismeta.dataType{end+1} = 'numeric';
gvObj.model.data.dsData.axis(1).axismeta = axismeta;

%% save gvArrayData
dynasimData = gvObj.model.data.dsData;
filePath = fullfile(pwd, 'gvArrayData.mat');
save(filePath, 'dynasimData') % save gvArray obj

%% run
gvObj.run();

end