#!/usr/bin/env python
from numpy import *
from matplotlib import *

def LIF_interneuron_network(no_cells,no_inputs,T0,e_rate,plot_opt):
    #UNTITLED2 Right-hand side of leaky-integrate-and-fire neuron, for integration with RK45.
    #T0 is in milliseconds.

    sim_name = "LIF_#dcells_#dinputs_#dms_#gHz"  + str(no_cells) + "_" + str(no_inputs) + "_" + str(T0) + "_" + str(e_rate)

    dt = .005

    bigT = floor(T0/dt)

    V_rest = -70 
    V_thresh = -52 
    V_reset = -85 
    V_spike = 0

    gE = 0.3 
    gI = 4

    tau_i1 = 1 
    tau_ir = 0.5 
    tau_id = 5 
    tau_i = 10 
    tau_r = 1
    tau_e1 = 1 
    tau_er = 0.5 
    tau_ed = 2 
    #tau_e = 20
    # dt_i = tau_i1/.005
    # dt_e = tau_e1/.005

    t = (1:bigT)*dt                     #Define time axis vector (useful for plotting).

    # spike_times = poiss_spike_times(T0/1000,1/e_rate,no_inputs) #Compute excitatory inputs.

    # Poisson input spikes at rate e_rate.
    spikes = random.rand(no_inputs,length(t))
    spikes = spikes[where(spikes < e_rate*dt/1000)]    #J: did i do this right

    # EPSP for spikes at time t = 0.   
    epsp = tau_i*(exp(-max(t - tau_e1,0)/tau_ed) - exp(-max(t - tau_e1,0)/tau_er))/(tau_ed - tau_er)
    #J: i am very confused about what the next two lines do
    epsp = epsp(epsp > eps)
    epsp = [zeros(1,length(epsp)) epsp]

    CE = rand(no_cells,no_inputs)    #Define connectivity from inputs to cells.
    CE = CE[where(CE<.1)]
    CE = CE*gE

    spike_arrivals = CE*spikes # Calculating presynaptic spikes for each cell.

    #J: what the heckie is this nan function
    epsps = nan(spike_arrivals.size) # Calculating EPSP experienced by each cell.
    for c in range(1, no_cells + 1):
        epsps[c,:] = convolve(spike_arrivals[c,:],epsp,'same')

    # for i in range(1, t.size + 1):
    #     cum_spikes[:,i] = sum(spike_times[where(spike_times< t[i])],axis=1)
    # 
    # e_spike_times = diff([cum_spikes cum_spikes[:,-1]],n=1,axis=1)   #J: is this a typo? wat

    CI = rand(no_cells,no_cells)     #Define connectivity between inhibitory cells.
    CI = CI[where(CI<.4)]
    CI = CI - diag(diag(CI))
    CI = CI*gI

    # IPSPs for spikes at time t = 0.
    ipsp = tau_i*(exp(-max(t - tau_i1,0)/tau_id) - exp(-max(t - tau_i1,0)/tau_ir))/(tau_id - tau_ir)
    #J: i do not understund
    ipsp = ipsp(ipsp > eps)
    ipsp = [zeros(1,length(ipsp)) ipsp]

    V = zeros(no_cells,bigT)                #Make empty variables to hold I-cell results.
    # s_i = zeros(no_cells,bigT)                #Make empty variables to hold the synapse results.
    # s_e = zeros(no_inputs,bigT)

    last_spike_times = -tau_r*ones(no_cells,1)
    # i_spike_arrivals = zeros(no_cells,bigT)
    ipsps = zeros(no_cells,bigT)
    # e_spike_times = t(end)*ones(no_inputs,1)

    V[:,1]= -70.0 + 10*rand(no_cells,1)  	#Set the initial conditions for I-cells.

    for i in range(1,T)                       #Integrate the equations.
    
        # LIF subthreshold dynamics.
        #st = V[:,i][where(V[:,i] < V_spike)]
        #V[st,i+1] = V[st,i] + dt*((V_rest - V[st,i])/tau_i - epsps[st,i]*V[st,i] + (CI[st,:]*ipsps[:,i])*(-70-V[st,i]))
        V[:,i+1] = V[:,i] + dt*((V_rest - V[:,i])/tau_i - epsps[:,i]*V[:,i] + (CI*ipsps[:,i])*(-70-V[:,i]))
    
        # LIF spiking.
        active_cells = (t(i+1) - last_spike_times) > tau_r
        spiking_cells = active_cells & (V[:,i+1] >= V_thresh)
        V(spiking_cells,i+1) = V_spike
        last_spike_times(spiking_cells) = t(i+1)
    
        #LIF reset.
        V(~active_cells, i+1) = V_reset
    
    #     # Handling spiking in I cells.
    #     if any(V[:,i+1] == V_spike)
    #         i_spike_arrivals(:, i+1) = CI*(V[:,i+1] == V_spike)
    #         isa_indices = find(i_spike_arrivals[:,i+1] > 0)
    # #         for c in range (1,length(isa_indices)+1)
    # #             ipsps(isa_indices(c),:) = conv(i_spike_arrivals(c,:), ipsp, 'same') 
    #         ipsps(isa_indices,:) = ipsps(isa_indices,:) + i_spike_arrivals(isa_indices,i+1)*conv([zeros(1,i) 1 zeros(1,T-i-1)],ipsp,'same')

    
    #     e_spike_times(spikes(:, max(i - dt_e,1)) == 1) = t(i)
    
    #     # Implements post-synaptic potentials as differential equation.    
    #     s_i_rhs = tau_i*(-exp(-max(t(i) - i_spike_times,0)/tau_id)/tau_id ...
    #         + exp(-max(t(i) - i_spike_times - tau_i1,0)/tau_ir)/tau_ir)/(tau_id - tau_ir) ...
    #         - (i_spike_times == t(end))*tau_i/(tau_id*tau_ir)
    #     
    #     s_e_rhs = tau_i*(-exp(-max(t(i) - e_spike_times,0)/tau_ed)/tau_ed ...
    #         + exp(-max(t(i) - e_spike_times - tau_e1,0)/tau_er)/tau_er)/(tau_ed - tau_er) ...
    #         - (e_spike_times == t(end))*tau_i/(tau_ed*tau_er)
    #     
    #     s_i[:,i+1] = s_i[:,i] + dt*s_i_rhs
    #     s_e[:,i+1] = s_e[:,i] + dt*s_e_rhs

    #     # Implements post-synaptic potentials as interrupting each other?
    #     s_i(i_spike_indices,:) = tau_i*(exp(-max(t - t(i+1)*ones(size(i_spike_indices)) - tau_i1,0)/tau_id) ...
    #         - exp(-max(t - t(i+1)*ones(size(i_spike_indices)) - tau_i1,0)/tau_ir))/(tau_id - tau_ir)
    #     
    #     s_e(e_spike_indices,:) = tau_i*(exp(-max(t - t(i+1)*ones(size(e_spike_indices)) - tau_e1,0)/tau_ed) ...
    #         - exp(-max(t - t(i+1)*ones(size(e_spike_indices)) - tau_e1,0)/tau_er))/(tau_ed - tau_er)

    #     # Implements PSPs as interrupting each other, differently.
    #     i_arg = max(t(i) - i_spike_times,0)
    #     s_i[:,i+1] = tau_i*(exp(-i_arg/tau_id) - exp(-i_arg/tau_ir))/(tau_id - tau_ir)
    #     
    #     e_arg = max(t(i) - e_spike_times,0)
    #     s_e[:,i+1] = tau_i*(exp(-e_arg/tau_ed) - exp(-e_arg/tau_er))/(tau_ed - tau_er)

        # Implements PSPs as sums over all past spikes.
        if any(V[:,i+1] >= V_thresh):
            i_spike_indices = find(V[:,i+1] >= V_thresh)
            ipsps(i_spike_indices,:) = ipsps(i_spike_indices,:) + repmat(conv([zeros(1,i) 1 zeros(1,T-i-1)],ipsp,'same'),length(i_spike_indices),1)
        end
    
    #     if any(spikes[:,i+1] == 1):
    #         e_spike_indices = find(spikes[:,i+1] == 1)
    #         s_e(e_spike_indices,(i+1):end) = s_e(e_spike_indices,(i+1):end) + epsp(e_spike_indices,1:T-i)

    save(sim_name,'V','ipsps','epsps','spikes','CI','CE','t')

    ## Plotting.

    if no_cells == 1 or plot_opt == 1:
        
        figure
    
        subplot(3,1,1)
        plot(t',V')
        axis tight
        box off
        title('I-cell Voltages')
    
        subplot(3,1,2)
        plot(t',ipsps')
        axis tight
        box off
        title('IPSPs')
    
        subplot(3,1,3)
        plot(t',epsps')
        axis tight
        box off
        title('EPSPs')
    
        save_as_pdf(gcf,[sim_name,'_Voltages'])

    if no_cells > 1 or plot_opt == 1:
    
        figure
    
        subplot(3,1,1)
    
        colormap('gray')
        imagesc(t,1:no_cells,1-(V>=V_thresh))
        title('Raster Plot')
    
        subplot(3,1,2)
    
        ax = plotyy(t',mean(V)',t',mean(ipsps)')
        axis tight
        box off
        title('Population Potentials')
        ylabel(ax(1),'Membrane Potential')
        ylabel(ax(2),'IPSP')
    
        subplot(3,1,3)
    
        plot(t',sum(V>=V_thresh)')
        axis tight
        box off
        title('Number of Spikes per Timestep')
    
        save_as_pdf(gcf,[sim_name,'_pop_mean'])
    
return [V, ipsps, epsps, spikes, CI, CE, t]
    
def main():
    [V, ipsps, epsps, spikes, CI, CE, t] = LIF_interneuron_network(no_cells,no_inputs,T0,e_rate,plot_opt)
    return
    
main()
