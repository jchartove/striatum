function current = poisson_input(no_cells, fraction_shared)
no_e_inputs = 127*no_cells; %127 = number of AMPA input synapses per cell in Hjorth et al, times 10 cells
e_rate = 10; 

no_i_inputs = 93*no_cells; % = 93 times 10
i_rate = 10;

tau_i1 = 1; tau_ir = 0.5; tau_id = 5; tau_i = 10;
tau_e1 = 1; tau_er = 0.5; tau_ed = 2;

CE_e = repmat(eye(no_cells), 1, no_e_inputs/no_cells);
CE_i = repmat(eye(no_cells), 1, no_i_inputs/no_cells);    %Define connectivity from inputs to cells.

% EPSP for spikes at time t = 0.
epsp = tau_i*(exp(-max(t - tau_e1,0)/tau_ed) - exp(-max(t - tau_e1,0)/tau_er))/(tau_ed - tau_er);
epsp = epsp(epsp > eps);    %?
epsp = [zeros(1,length(epsp)) epsp]; %?

% IPSP for spikes at time t = 0.
ipsp = tau_i*(exp(-max(t - tau_i1,0)/tau_id) - exp(-max(t - tau_i1,0)/tau_ir))/(tau_id - tau_ir);
ipsp = ipsp(ipsp > eps);
ipsp = [zeros(1,length(ipsp)) ipsp];

e_size = 0.0025*k;
e_spikes = rand(floor((1-fraction_shared)*(no_e_inputs/no_cells))*no_cells,length(t));
e_spikes = e_spikes < e_rate*dt/1000;

i_size = 0.0025*k;
i_spikes = rand(floor((1-fraction_shared)*(no_i_inputs/no_cells))*no_cells,length(t));
i_spikes = i_spikes < i_rate*dt/1000;

%Define connectivity from inputs to cells.
CE_e = repmat(eye(no_cells), 1, floor((1-fraction_shared)*(no_e_inputs/no_cells)));
e_spike_arrivals = CE_e*e_spikes; % Calculating presynaptic spikes for each cell.
epsps = nan(size(e_spike_arrivals)); % Calculating EPSP experienced by each cell.

CE_i = repmat(eye(no_cells), 1, floor((1-fraction_shared)*(no_i_inputs/no_cells))); 
i_spike_arrivals = CE_i*i_spikes; % Calculating presynaptic spikes for each cell.
ipsps = nan(size(i_spike_arrivals)); % Calculating EPSP experienced by each cell.

if ~isempty(epsps)
    for c = 1:no_cells
        epsps(c,:) = e_size*conv(e_spike_arrivals(c,:),epsp,'same');
    end
else
    for c = 1:no_cells
        epsps(c,:) = zeros(no_cells,T);
    end
end


if ~isempty(ipsps)
    for c = 1:no_cells
        ipsps(c,:) = i_size*conv(i_spike_arrivals(c,:),ipsp,'same');
    end
else
    for c = 1:no_cells
        ipsps(c,:) = zeros(no_cells,T);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %shared input
CE_e_shared = ones(no_cells, floor((fraction_shared)*(no_e_inputs/no_cells)));
shared_e_spikes = rand(floor((fraction_shared)*(no_e_inputs/no_cells)),length(t));
shared_e_spikes  = shared_e_spikes < e_rate*dt/1000;
shared_e_arrivals = CE_e_shared*shared_e_spikes;
shared_e_input = nan(size(shared_e_arrivals)); % Calculating EPSP experienced by each cell. 
if ~isempty(shared_e_input)
    for c = 1:no_cells
        shared_e_input(c,:)= e_size*conv(shared_e_arrivals(c,:),epsp,'same'); 
    end
else
    shared_e_input = zeros(no_cells,T);
end

CE_i_shared = ones(no_cells, floor((fraction_shared)*(no_i_inputs/no_cells)));
shared_i_spikes = rand(floor((fraction_shared)*(no_i_inputs/no_cells)),length(t));
shared_i_spikes  = shared_i_spikes < i_rate*dt/1000;
shared_i_arrivals = CE_i_shared*shared_i_spikes;
shared_i_input = nan(size(shared_i_arrivals)); % Calculating IPSP experienced by each cell. 
if ~isempty(shared_e_input)
    for c = 1:no_cells
        shared_i_input(c,:)= i_size*conv(shared_i_arrivals(c,:),ipsp,'same');  
    end
else
    shared_i_input = zeros(no_cells,T);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:no_cells, epsps(i, :) = shared_e_input(1,:) + epsps(i,:); end
for i = 1:no_cells, ipsps(i, :) = shared_i_input(1,:) + ipsps(i,:); end

current = epsps-ipsps;
end