eqns={ 'dV/dt = (Iapp + @current )/Cm;I=0; Cm=1; V(0)=-63 + 63.*rand(1,Npop)';};

% for 
    numcells = [100]
spec=[];
T0 = 4000;

ncells = 100;  % number of MSN cells in the pool
g_gaba = 0.1/(ncells-1); % recurrent gaba conductance, normalized to the number of cells
g_m = 1.2; % 1.2; % 1.3; % 1.2 parkinsonian, 1.3 normal
vrand = 63*rand(1,ncells);

spec.nodes(1).name = 'D1';
spec.nodes(1).size = ncells;
spec.nodes(1).equations = eqns;
spec.nodes(1).mechanism_list = {'naCurrentMSN','kCurrentMSN','mCurrentMSN','leakCurrentMSN','injectedCurrentD1','noisyInputMSN'};
spec.nodes(1).parameters = {'cm',1,'V_IC',-63,'g_m',g_m,'Tfinal', T0, 'Iapp',0}; 

spec.nodes(2).name = 'D2';
spec.nodes(2).size = ncells;
spec.nodes(2).equations = eqns;
spec.nodes(2).mechanism_list = {'naCurrentMSN','kCurrentMSN','mCurrentMSN','leakCurrentMSN','injectedCurrentD2','noisyInputMSN'};
spec.nodes(2).parameters = {'cm',1,'V_IC',-63,'g_m',g_m,'Tfinal', T0, 'Iapp',0}; 

spec.connections(1).direction = 'D1->D1';
spec.connections(1).mechanism_list = {'gabaRecInputMSN'};
spec.connections(1).parameters = {'g_gaba',g_gaba};

spec.connections(2).direction = 'D2->D2';
spec.connections(2).mechanism_list = {'gabaRecInputMSN'};
spec.connections(2).parameters = {'g_gaba',g_gaba};

vary={
  %'(D1, D2)', 'DA',	[-1:2];
  '(D2)','injectedCurrent',[1:0.01:1.3]
};

namearray = cellfun(@num2str,vary,'UniformOutput',0);
namestr = strjoin(reshape(namearray, 1, []));
cd '/projectnb/crc-nak/chartove/dynasim/'; 
memlimit = '64G';
cluster_flag = 1;
overwrite_flag = 0;
save_data_flag = 1;
save_results_flag = 1;
verbose_flag = 1;
compile_flag = 0;
disk_flag = 0;
downsample_factor = 10;

dsSimulate(spec,...
              'analysis_functions', {@gvFRsoma, @gvCalcPower},...
              'save_data_flag',save_data_flag,'study_dir','spn_only_granular',...
              'cluster_flag',cluster_flag,'verbose_flag',verbose_flag,...
              'overwrite_flag',overwrite_flag,'tspan',[0 T0],...
              'save_results_flag',save_results_flag,'solver','rk4',...
              'memlimit',memlimit,'compile_flag',compile_flag,...
				'copy_run_file_flag',1, 'copy_mech_files_flag',1, ...
              'disk_flag',disk_flag,'downsample_factor',downsample_factor,...
              'vary',vary, 'dt', .01);