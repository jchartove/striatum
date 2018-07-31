clear

eqns={
  'dV/dt=Iapp+@current';
};

 numcells = [2]
spec=[];
T0 = 4000;
spec.nodes(1).name = 'soma';
spec.nodes(1).size = numcells;
spec.nodes(1).equations = eqns;
spec.nodes(1).mechanism_list = {'somaGolombK','somaGolombKdr','somaInput','somaGolombNa','somaLeak'}; 
spec.nodes(1).parameters = {'v_IC',-90+90*rand(1,numcells), 'Tfinal', T0, 'Iapp',0};

spec.nodes(2).name = 'dend';
spec.nodes(2).size = numcells;
spec.nodes(2).equations = eqns;
spec.nodes(2).mechanism_list = {'dendGolombK','dendGolombKdr','dendGolombNa','dendInput','dendLeak','dendiMultiPoissonExp'};
spec.nodes(2).parameters = {'v_IC',-90+90*rand(1,numcells), 'Tfinal', T0, 'Iapp',0}; 

spec.connections(1).direction = 'soma->soma';
spec.connections(1).mechanism_list = {'somaSomaiSYN'};
spec.connections(1).parameters = {'Tfinal', T0};

spec.connections(2).direction = 'soma->dend';
spec.connections(2).mechanism_list = {'iCOM'};
spec.connections(2).parameters = {'gCOM',.15};

spec.connections(3).direction = 'dend->soma';
spec.connections(3).mechanism_list = {'iCOM'};
spec.connections(3).parameters = {'gCOM', .15};

spec.connections(4).direction = 'dend->dend';
spec.connections(4).mechanism_list = {'dendDendiGAP'};
spec.connections(4).parameters = {'Tfinal', T0};

vary={
  '(dend)',			'tonic',	[0:1:10];
  '(dend)',			'rate',	[0:0.2:4];
  '(dend)',			'DA',	[0];
  '(dend-dend)',			'g_GAP',	[0];
  '(soma-soma)',			'gsyn',	[0];
};

namearray = cellfun(@num2str,vary,'UniformOutput',0);
namestr = strjoin(reshape(namearray, 1, []));
cd '/projectnb/crc-nak/chartove/dynasim/'; %try to cd to this directory and leave data_dir blank
memlimit = '64G';
cluster_flag = 1;
overwrite_flag = 1;
save_data_flag = 1;
save_results_flag = 1;
verbose_flag = 1;
compile_flag = 0;
disk_flag = 0;
downsample_factor = 10;

[~,~]=dsSimulate(spec,...
              'analysis_functions', {@get_fft},...
              'save_data_flag',save_data_flag,'study_dir','new_noise_sim_again',...
              'cluster_flag',cluster_flag,'verbose_flag',verbose_flag,...
              'overwrite_flag',overwrite_flag,'tspan',[0 T0],...
              'save_results_flag',save_results_flag,'solver','rk4',...
              'memlimit',memlimit,'compile_flag',compile_flag,...
              'disk_flag',disk_flag,'downsample_factor',downsample_factor,...
              'vary',vary, 'dt', .01);