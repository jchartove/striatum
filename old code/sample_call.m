input = rand(10, ceil(100/.005));
shared_input = rand(1, ceil(100/.005));
input = input < 1000*.005/1000;  %0.005
shared_input = shared_input < 1000*.005/1000;
for i = 1:10, input(i, :) = conv(double(input(i,:)), ones(1,200), 'same'); end
shared_input(1, :) = conv(double(input(1,:)), ones(1,200), 'same');
fraction_shared = 0.5;
for i = 1:10, input(i, :) = fraction_shared*shared_input + (1-fraction_shared)*input(i,:); end

[Vs,Vd,s,m,h,n,t] = ing_w_dendritic_gap_jxn(10, input*70, 100, [], zeros(10), .3*(rand(10) > .4));

plot(Vs')

%INPUTS:
%  10 = no_cells = number of simulated FSIs.
%  input*70 = I0  = Input to cell. If a scalar, treated as constant current, identical for each cell. If a
%  vector, must be length ceil(T0/0.005); treated as identical for each cell.
%  If different input is given to each cell, must be a matrix, of
%  dimensions (no_cells, ceil(T0/0.005)). 
%  100 = T0  = total time of simulation [ms].
%  [] = g_sd = conductance between dendrite and soma (default = leak
%  conductance/2).
%  zeros(10) = CS = inhibitory synapse connectivity matrix, including conductance (strength) for each synapse.
%  .3*(rand(10) > .4) = CG = gap junction connectivity matrix, including conductance (strength) for each gap junction.
% function [Vs,Vd,s,m,h,n,t] = ing_w_dendritic_gap_jxn(no_cells, I0, T0, g_sd, CS, CG)

%OUTPUTS:
%  Vs = voltage of I-cell soma.
%  Vd = voltage of I-cell dendrite.
%  s = inhibitory synapse.
%  m = activation gating variable of Na channel.
%  h = inactivation gating variable of Na channel.
%  n = gating variable of K channel.
%  t = time axis vector (useful for plotting).

%sample calls:
%no chemical synapses, increasing gap conductivity decreases firing rate

%if inputs are correlated, there will not be as much of a decrease in firing rate

%to find number of action potentials: number of time points when soma voltage is above 0 is a rough estimate...
firing_estimate = sum(sum(Vs > 0));
%i guess a more specific way to do it is to count the number of times a negative number is followed by a positive number

