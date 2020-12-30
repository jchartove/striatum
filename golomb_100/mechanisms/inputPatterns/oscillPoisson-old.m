 function psps =oscillPoisson(no_cells, inputs_per_cell, rate, tau_i, tau_1, tau_d, tau_r, T, dt, fraction_shared, fraction_oscill, randstart, freq)

t = 0:dt:T;

tau_i = tau_i/freq; %normalized so amount of excitation is same regardless of frequency
% EPSP shape for spikes at time t = 0.
psp = tau_i*(exp(-max(t - tau_1,0)/tau_d) - exp(-max(t - tau_1,0)/tau_r))/(tau_d - tau_r);
psp = psp(psp > eps);    
psp = [zeros(1,length(psp)) psp]; 

%no_inputs = inputs_per_cell*no_cells;
%the previous line isn't commented in corrPoisson

%from here to line 45 is the new part

fraction_ind = 1-(fraction_shared+fraction_oscill);

%sin(t*dt/(2*pi)) will make a sinusoid with one cycle over the whole time
wave = sin(freq*T*t*dt/(2*pi*1000));
C_oscill = ones(no_cells, floor(fraction_oscill*inputs_per_cell)); %C stands for connectivity
num_o_inputs = floor(fraction_oscill*inputs_per_cell)*no_cells;
o_spikes = rand(num_o_inputs,length(t)); %rand matrix
o_spikes = o_spikes < rate*dt/1000; %rands occurring at frequency rate
%oscill = zeros(floor(fraction_oscill*inputs_per_cell),length(t));
%cycles = freq*(T/1000);
%wavelength = floor(length(t)/cycles);

%if randstart %this doesn't actually work
%    start = ceil(wavelength*rand);
%else
%    start = 1;
%end

for n = 1:num_o_inputs
    o_spikes(n,:) = conv(o_spikes(n,:),wave,'same');
end

oscill_arrivals = C_oscill*o_spikes;
%is the point of C_oscill just to get it in the right dimensionality

oscill_in = nan(size(oscill_arrivals)); %should be nan
% Calculating EPSP experienced by each cell.
if ~isempty(oscill_in)
    for c = 1:no_cells
        oscill_in(c,:) = conv(oscill_arrivals(c,:),psp,'same'); %convolve psp shape
    end
else
    oscill_in(c,:) = zeros(no_cells,T);
end

%old part resumes

%connectivity matrix: ones along diagonal
C = repmat(eye(no_cells), 1, floor(fraction_ind*inputs_per_cell));
spikes = rand(floor(fraction_ind*inputs_per_cell)*no_cells,length(t)); %rand matrix
spikes = spikes < rate*dt/1000; %rands occurring at frequency rate

spike_arrivals = C*spikes; % Calculating presynaptic spikes for each cell.

psps = nan(size(spike_arrivals)); %matrix of nans
% Calculating EPSP experienced by each cell.
if ~isempty(psps)
    for c = 1:no_cells
        psps(c,:) = conv(spike_arrivals(c,:),psp,'same'); %convolve psp shape
    end
 else
    for c = 1:no_cells
        psps(c,:) = zeros(no_cells,T);
    end
 end

 C_shared = ones(no_cells, floor((fraction_shared)*inputs_per_cell)); %appropriately sized ones matrix
 shared_spikes = rand(floor((fraction_shared)*inputs_per_cell),length(t)); %rand matrix
 
 shared_spikes = shared_spikes < rate*dt/1000; %rands occurring at frequency rate -> binary
 shared_arrivals = C_shared*shared_spikes; %binary arrival matrix
 shared_input = nan(size(shared_arrivals)); %matrix of nans
 
 % Calculating EPSP experienced by each cell. 
 if ~isempty(shared_input)
    for c = 1:no_cells
        shared_input(c,:)= conv(shared_arrivals(c,:),psp,'same'); %convolve psp shape
    end
 else
    shared_input = zeros(no_cells,T);
 end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
for i = 1:no_cells
    psps(i, :) = [shared_input(1,:) + psps(i,:) + oscill_in(i,:)]; 
end