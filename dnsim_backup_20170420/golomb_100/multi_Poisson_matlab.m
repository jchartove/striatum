function psps = multi_Poisson_matlab(no_cells, inputs_per_cell, rate, tau_i, tau_1, tau_d, tau_r, T, dt)

t = 0:dt:T;

% EPSP for spikes at time t = 0.
psp = tau_i*(exp(-max(t - tau_1,0)/tau_d) - exp(-max(t - tau_1,0)/tau_r))/(tau_d - tau_r);
psp = psp(psp > eps);    %?
psp = [zeros(1,length(psp)) psp]; %?

no_inputs = inputs_per_cell*no_cells;

spike_arrivals = poissrnd(rate*no_inputs, no_cells, length(t));

psps = nan(size(spike_arrivals)); % Calculating EPSP experienced by each cell.
for c = 1:no_cells
    psps(c,:) = conv(spike_arrivals(c,:),psp,'same');
end