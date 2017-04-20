clear
% Model: golomb_activedend_10
cd '/project/crc-nak/jchartove/dynasim/golomb_100';

for numcells = [100]
spec=[];
T0 = 2000;
spec.nodes(1).label = 'soma';
spec.nodes(1).multiplicity = numcells;
spec.nodes(1).dynamics = {'v''=current'};
spec.nodes(1).mechanisms = {'soma_golomb_K','soma_golomb_Kdr','soma_golomb_Na','soma_input','soma_leak'};
spec.nodes(1).parameters = {'v_IC',-90+90*rand(1,numcells), 'Tfinal', T0};

spec.nodes(2).label = 'dend';
spec.nodes(2).multiplicity = numcells;
spec.nodes(2).dynamics = {'v''=current'};
spec.nodes(2).mechanisms = {'dend_golomb_K','dend_golomb_Kdr','dend_golomb_Na','dend_input','dend_leak','dend_iMultiPoissonExp'};
spec.nodes(2).parameters = {'v_IC',-90+90*rand(1,numcells), 'Tfinal', T0}; 

ncells = 100;  % number of MSN cells in the pool
g_gaba = 0.1/(ncells-1); % recurrent gaba conductance, normalized to the number of cells
g_m = 1.2; % 1.2; % 1.3; % 1.2 parkinsonian, 1.3 normal
%V_ic = -63;
vrand = 63*rand(1,ncells);

spec.nodes(3).label = 'D1';
spec.nodes(3).multiplicity = ncells;
spec.nodes(3).dynamics = {'V''=(current)./cm'};
spec.nodes(3).mechanisms = {'naCurrentMSN','kCurrentMSN','mCurrentMSN','leakCurrentMSN','injectedCurrentD1','noisyInputMSN'};
spec.nodes(3).parameters = {'cm',1,'V_IC',-63,'g_m',g_m,'Tfinal', T0}; % V_IC refers to the initial condition for the membrane 

spec.nodes(4).label = 'D2';
spec.nodes(4).multiplicity = ncells;
spec.nodes(4).dynamics = {'V''=(current)./cm'};
spec.nodes(4).mechanisms = {'naCurrentMSN','kCurrentMSN','mCurrentMSN','leakCurrentMSN','injectedCurrentD2','noisyInputMSN'};
spec.nodes(4).parameters = {'cm',1,'V_IC',-63,'g_m',g_m,'Tfinal', T0}; % V_IC refers to the initial condition for the membrane potential


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

spec.connections(4,4).label = [spec.nodes(4).label,'-',spec.nodes(4).label];
spec.connections(4,4).mechanisms = {'gabaRecInputMSN'};
spec.connections(4,4).parameters = {'g_gaba',g_gaba};

spec.connections(1,3).label = 'soma-D1';
spec.connections(1,3).mechanisms = {'soma_MSN_iSYN'};
spec.connections(1,3).parameters = {'gsyn',6*g_gaba};

spec.connections(1,4).label = 'soma-D2';
spec.connections(1,4).mechanisms = {'soma_MSN_iSYN'};
spec.connections(1,4).parameters = {'gsyn',6*g_gaba};

%dnsim(spec); % open model in DNSim GUI

% DNSim simulation and plots:
%data = runsim(spec,'timelimits',[0 T0],'dt',.01,'SOLVER','rk4','timesurfer_flag',0,'savedata_flag',0); % simulate DNSim models

%model=buildmodel(spec); % parse DNSim spec structure

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

% scope = {'(soma,dend)','(soma,dend)','(soma,dend)'};
% variable = {'gna','gkdr','gd'};
% values = {'[100:50:250]','[200:50:450]','[6:7]'};

% scope = {'soma-soma','dend-dend', 'dend'};
% variable = {'gsyn','g_GAP', 'fraction_shared'};
% values = {'[0:0.005:0.01]','[0:0.005:0.01]','[0:0.2:1]'};

% scope = {'dend-dend','dend-dend','soma-soma'};
% variable = {'gcon','g_GAP','i_con'};
% values = {'[1]','[0:0.001:0.01,0.02:0.01:0.1,0.2:0.1:0.5]','[0]'};

% scope = {'(soma,dend)','dend-dend','dend','dend-dend','soma-soma'};
% variable = {'tau_mult','g_GAP','fraction_shared','gcon','i_con'};
% values = {'[1:5:21]','[0:0.5:5]','[0:0.2:1]','[1]','[0]'};

% scope = {'(soma,dend)','dend-dend','soma-soma','dend'};
% variable = {'tau_mult','g_GAP','i_con','fraction_shared'};
% values = {'[1:5:21]','[0:0.1:0.5]','[0,0.3]','[0:0.2:1]'};

%val = strcat('[', num2str(numcells), ']')
%scope = {'soma-soma','dend', 'dend','(soma,dend)','soma-soma','(soma,dend)'};
%variable = {'numcells', 'tonic','rate','taub', 'tauD', 'tau_mult'}; %g_esyn
%values = {val, '[0,5]','[0,2]','[50:50:400]','[1:2:15]','[0.5:0.5:2]'};

% scope = {'soma-soma','soma-soma','dend-dend'};
% variable = {'tauD','gsyn','g_GAP'};
% values = {'[1:6]','[0:0.01:0.3]','[0,0.05]'};

%scope = {'dend','(soma,dend)','dend-dend','soma-soma'};
%variable = {'tonic', 'gd','g_GAP','gsyn'};
%values = {'[0:2:20]','[4:8]', '[0.02,0.37]','[0.005,0.1,0.3]'};

% val = strcat('[', num2str(numcells), ']')
% scope = {'soma-soma','(soma,dend)','soma-soma','(soma,dend)','dend', 'dend'};
% variable = {'numcells','taub', 'tauD', 'tau_mult','tonic', 'rate'};
% values = {val,'[50:50:400]','[1:2:15]','[0.5:0.5:2]','[10]','[0,2]'};

%scope = {'dend-dend','soma-soma','dend'};
%variable = {'g_GAP','gsyn','fraction_shared'};
%values = {'[0.12]','[0.2]','[0.6]'};

 %scope = {'(soma,dend)','dend-dend'};
 %variable = {'thetam','g_GAP'};
 %values = {'[-22]','[0]'};

%play with thetam and gd in the single cell at some point

% scope = {'(soma,dend)','(soma,dend)','(soma,dend)','dend-dend','soma-soma'};
% variable = {'gna','gkdr','gd','g_GAP','gsyn'};
% values = {'[100,112]','[200,225]','[8]','[0:0.02:0.1]','[0,0.002]'};


% scope = {'(soma,dend)','(soma,dend)','(soma,dend)','dend-dend','soma-soma'};
% variable = {'gd','gkdr','gna','gGAP','gsyn'};
% values = {'[6:2:10]','[200:100:400]','[100:100:300]','[0:0.02:0.1]','[0:0.002:0.01]'};

%dopamine condition: tonic = 10, GJ = 0.4, inhib = 0.005
%low dopamine: tonic = 0, GJ = 0.08, inhib = 0.08

% scope = {'(dend,dend-dend,soma-soma)','dend','soma-soma','soma-MSN'};
% variable = {'DA','tonic','gsyn','gsyn'};
% values = {'[0]','[0:10]','[0:0.005:0.08]','[0.1]'};

% scope = {'MSN','MSN'}
% variable = {'injectedCurrent','sigma_noise'}
% values = {'[1:0.5:5]','[4:10]'}

% scope = {'soma-MSN','soma-MSN'}
% variable = {'ko','i_con'}
% values = {'[0:10:100]','[0:0.1:1]'}

% scope = {'MSN','MSN','soma-MSN','soma-MSN'}
% variable = {'freq','injectedCurrent','i_con','gsyn'}
% values = {'[1:80]','[1.0:0.2:1.8]','[0]','[0]'}

scope = {'(dend,dend-dend,soma-soma,D1,D2)','(soma-D1,soma-D2)','(D1,D2)','(soma,dend)'};
variable = {'DA','gsyn','g_m','taub'};
values = {'[0:1]','[0.01]','[1.3]','[120]'};

%scope = {'(dend,dend-dend,soma-soma)','(soma-D1,soma-D2)','(D1,D2)','(D1,D2)'};
%variable = {'DA','gsyn','g_m','injectedCurrent'};
%values = {'[0:0.1:1]','[0.01]','[1.2,1.3]','[1.1:0.01:1.3]'};

%[0:0.0005:0.01]

% scope = {'(dend,dend-dend,soma-soma)','dend','soma-soma'};
% variable = {'DA','tonic','tauD'}
% values = {'[0,1]','[0:20]','[9:13]'}

[~,~]=SimulateModel(spec,scope,variable,values,...
  'dt',.01,'SOLVER','rk4', 'memlimit', '64G', 'overwrite_flag', 1, 'timelimits',[0 T0], 'rootdir', '/projectnb/crc-nak/chartove/dnsim', ...
  'savedata_flag',1,'timesurfer_flag',0,'saveplot_flag',0,'plotvars_flag',0,'addpath','/project/crc-nak/jchartove/dnsim',...
  'cluster_flag',1);
end