clear
% Model: wilson_activedend_10
cd '/project/crc-nak/jchartove/dnsim/golomb_100/wilson';

for numcells = [100]
spec=[];
T0 = 2000;
spec.nodes(1).label = 'soma';
spec.nodes(1).multiplicity = numcells;
spec.nodes(1).dynamics = {'v''=current'};
spec.nodes(1).mechanisms = {'soma_wilson_K','soma_wilson_Kdr','soma_wilson_Na','soma_input','soma_leak'};
spec.nodes(1).parameters = {'v_IC',-90+90*rand(1,numcells), 'Tfinal', T0};

spec.nodes(2).label = 'dend';
spec.nodes(2).multiplicity = numcells;
spec.nodes(2).dynamics = {'v''=current'};
spec.nodes(2).mechanisms = {'dend_wilson_K','dend_wilson_Kdr','dend_wilson_Na','dend_input','dend_leak','dend_iMultiPoissonExp'};
spec.nodes(2).parameters = {'v_IC',-90+90*rand(1,numcells), 'Tfinal', T0}; 

ncells = 2;  % number of MSN cells in the pool
g_gaba = 0.1/(ncells-1); % recurrent gaba conductance, normalized to the number of cells
g_m = 1.2; % 1.2; % 1.3; % 1.2 parkinsonian, 1.3 normal
V_ic = -63; % -63+5*randn(1,ncells); % Initial conditions for the phases

spec.nodes(3).label = 'MSN';
spec.nodes(3).multiplicity = ncells;
spec.nodes(3).dynamics = {'V''=(current)./cm'};
spec.nodes(3).mechanisms = {'naCurrentMSN','kCurrentMSN','mCurrentMSN','leakCurrentMSN','injectedCurrentMSN','noisyInputMSN'};
spec.nodes(3).parameters = {'cm',1,'V_IC',V_ic,'g_m',g_m,'Tfinal', T0}; % V_IC refers to the initial condition for the membrane potential

spec.connections(1,1).label = 'soma-soma';
spec.connections(1,1).mechanisms = {'soma_soma_iSYN'};
spec.connections(1,1).parameters = [];

spec.connections(1,2).label = 'soma-dend';
spec.connections(1,2).mechanisms = {'soma_dend_iCOM'};
spec.connections(1,2).parameters = [];

spec.connections(2,1).label = 'dend-soma';
spec.connections(2,1).mechanisms = {'dend_soma_iCOM'};
spec.connections(2,1).parameters = [];

spec.connections(2,2).label = 'dend-dend';
spec.connections(2,2).mechanisms = {'dend_dend_iGAP'};
spec.connections(2,2).parameters = [];

spec.connections(3,3).label = [spec.nodes(3).label,'-',spec.nodes(3).label];
spec.connections(3,3).mechanisms = {'gabaRecInputMSN'};
spec.connections(3,3).parameters = {'g_gaba',g_gaba};

spec.connections(1,3).label = 'soma-MSN';
spec.connections(1,3).mechanisms = {'soma_MSN_iSYN'};
spec.connections(1,3).parameters = [];

%dnsim(spec); % open model in DNSim GUI

% DNSim simulation and plots:
%data = runsim(spec,'timelimits',[0 T0],'dt',.01,'SOLVER','rk4','timesurfer_flag',0,'savedata_flag',0); % simulate DNSim models
model=buildmodel(spec); % parse DNSim spec structure

% scope = {'(soma,dend)','dend-dend','dend'};
% variable = {'gd','g_GAP','g_esyn'};
% values = {'[0:0.2:2]','[0:0.2:1]','[10]'};

% for n = 1:3
% %Sweep over parameter values:
% gd = (n-1)*4;
% rate = (gd+1)/4;
% gdval = strcat('[',num2str(gd),']');
% rateval = strcat('[',num2str(rate),']');
% values = {gdval,'[0:0.1:1]','[1]','[0]',rateval};

%val = strcat('[', num2str(numcells), ']')
%scope = {'soma-soma','dend', 'dend','(soma,dend)','soma-soma','(soma,dend)'};
%variable = {'numcells', 'tonic','rate','taub', 'tauD', 'tau_mult'}; %g_esyn
%values = {val, '[0,5]','[0,2]','[50:50:400]','[1:2:15]','[0.5:0.5:2]'};

%play with thetam and gd in the single cell at some point

%dopamine condition: tonic = 10, GJ = 0.4, inhib = 0.005
%low dopamine: tonic = 0, GJ = 0.08, inhib = 0.08

scope = {'(dend,dend-dend,soma-soma)','(soma,dend)','(dend-soma,soma-dend)'};
variable = {'DA','gd','gcom'}
values = {'[0,1]','[3:8]','[0.15,4]'}

[~,~,outdir]=simstudy(model,scope,variable,values,...
  'dt',.01,'SOLVER','rk4', 'memlimit', '64G', 'overwrite_flag', 1, 'timelimits',[0 T0], 'rootdir', '/projectnb/crc-nak/chartove/dnsim', ...
  'savedata_flag',1,'timesurfer_flag',0,'saveplot_flag',0,'plotvars_flag',0,'addpath','/project/crc-nak/jchartove/dnsim',...
  'cluster_flag',1);
end