#!/usr/bin/env python
import numpy as np
import random

"""
Implement a model of Interneuron Network Gamma (ING).
Use the Hodgkin-Huxley equation for the individual neuron.
Use Borgers et al, PNAS, 2008 to model the synapses
  (see http://www.pnas.org/content/105/46/18023.abstract).

INPUTS:
  no_cells = number of simulated FSIs.
  I0  = Input to cell. If a scalar, treated as constant current, identical for each cell. If a
  vector, must be length ceil(T0/0.005); treated as identical for each cell.
  If different input is given to each cell, must be a matrix, of dimensions (no_cells, ceil(T0/0.005)). 
  T0  = total time of simulation [ms].
  g_sd = conductance between dendrite and soma (default = leak conductance/2).
  CS = inhibitory synapse connectivity matrix, including conductance (strength) for each synapse.
  CG = gap junction connectivity matrix, including conductance (strength) for each gap junction.

OUTPUTS:
  Vs = voltage of I-cell soma.
  Vd = voltage of I-cell dendrite.
  s = inhibitory synapse.
  m = activation gating variable of Na channel.
  h = inactivation gating variable of Na channel.
  n = gating variable of K channel.
  t = time axis vector (useful for plotting).
"""

def ing_w_dendritic_gap_jxn(no_cells, I0, T0, g_sd, CS, CG):
    dt = 0.005                      #The time step.
    T  = np.ceil(T0/dt)
    tauI = 12                       #Decay time of inhibition.
    gNa = 120  
    ENa=115                         #Sodium max conductance and reversal.    
    gK = 36  
    EK=-12                          #Potassium max conductance and reversal.
    gL = 0.3  
    ERest=10.6                      #Leak max conductance and reversal.
    if g_sd == 0:
        g_sd = gL/2                 #Conductance from dendrite to soma (half of leak conductance, as in Lewis & Rinzel 2003).
    littlet = range(1,T+1)*dt                    #Define time axis vector (useful for plotting).
  
    Vs = np.zeros(no_cells,T)          #Make empty variables to hold I-cell results.
    Vd = np.zeros(no_cells,T)
    m = np.zeros(no_cells,T)
    h = np.zeros(no_cells,T)
    n = np.zeros(no_cells,T)

    s = np.zeros(no_cells,T)               #Make empty variables to hold the synapse results.
  
    Vs[:,1]=-70+70*np.random.rand(no_cells,1)  	#Set the initial conditions for P-cells, I-cells, and synapses.
    Vd[:,1]=Vs[:,1]
    m[:,1]=np.random.rand(no_cells,1)
    h[:,1]=0.5*np.random.rand(no_cells,1)
    n[:,1]=0.35 + .4*np.random.rand(no_cells,1)
    s[:,1]=0.0 + 0.1*np.random.rand(no_cells,1)
      
    I0 = I0_arg_check(I0, no_cells, littlet)
  
    for i in range(1, T):                   #Integrate the equations.
        print "Time ", i
        Vs[:,i+1] = Vs[:,i] + dt*(gNa*(m[:,i]**3)*h[:,i]*(ENa-(Vs[:,i]+65)) + gK*(n[:,i]**4)*(EK-(Vs[:,i]+65)) + gL*(ERest-(Vs[:,i]+65)) + I0[:,i] + CS*s[:,i]*(-80-Vs[:,i]) + g_sd*(Vd[:,i]-Vs[:,i])) #Updating soma voltage: Sodium & potassium currents, Leak & applied currents, Synaptic & dendritic currents.                                    
        Vd[:,i+1] = Vd[:,i] + dt*(g_sd*(Vs[:,i]-Vd[:,i]) + (CG*np.diag(Vd[:,i])-np.diag(Vd[:,i])*CG)*np.ones(no_cells,1))         #Updating dendrite voltage.            
        m[:,i+1] = m[:,i] + dt*(alphaM(Vs[:,i])*(1-m[:,i]) - betaM(Vs[:,i])*m[:,i])                                    #Update m.
        h[:,i+1] = h[:,i] + dt*(alphaH(Vs[:,i])*(1-h[:,i]) - betaH(Vs[:,i])*h[:,i])                                    #Update h.
        n[:,i+1] = n[:,i] + dt*(alphaN(Vs[:,i])*(1-n[:,i]) - betaN(Vs[:,i])*n[:,i])                                    #Update n.
        s[:,i+1] = s[:,i] + dt*(((1+tanh(Vs[:,i]/10))/2)*(1-s[:,i])/0.5 - s[:,i]/tauI)                                  #Update s.
    
    return [Vs,Vd,s,m,h,n,littlet]

#Below, define the auxiliary functions alpha & beta for each gating variable.

def alphaM(V):
    return (2.5 - 0.1*(V + 65)) / (np.exp(2.5 - 0.1*(V + 65)) - 1)

def betaM(V):
    return 4*np.exp(-(V + 65)/18)

def alphaH(V):
    return 0.07*np.exp(-(V + 65)/20)

def betaH(V):
    return 1 / (np.exp(3.0 - 0.1*(V + 65)) + 1)

def alphaN(V):
    return (0.1 - 0.01*(V + 65)) / (np.exp(1 - 0.1*(V + 65)) - 1)

def betaN(V):
    return 0.125*np.exp(-(V + 65)/80)
    
def I0_arg_check(I0, no_cells, t):
    if np.isscalar(I0):
        I0_out = I0*np.ones(no_cells, length(t))
        return I0_out
    else:
        [r, c] = I0.shape
        if r == 1 or c == 1:
            if c == 1:
                I0 = np.transpose(I0)
            if max(r, c) == t.size:
                I0_out = tile(I0, (no_cells, 1))
                return I0_out
            else:
                print "I0 must be scalar, 1 x ceil(T0/0.005), or no_cells x ceil(T0/0.005)."
                return        
        elif r == no_cells or c == no_cells:
            if c == no_cells:
                I0 = np.transpose(I0)
                [r, c] = I0.shape
            if c != t.size:
                print "I0 must be scalar, 1 x ceil(T0/0.005), or no_cells x ceil(T0/0.005)."
                return
            else:
                I0_out = I0
                return I0_out
        else:
            print "I0 must be scalar, 1 x ceil(T0/0.005), or no_cells x ceil(T0/0.005)."
            return

def main():
    [Vs,Vd,s,m,h,n,t] = ing_w_dendritic_gap_jxn(no_cells, I0, T0, g_sd, CS, CG)
    return
    
main()