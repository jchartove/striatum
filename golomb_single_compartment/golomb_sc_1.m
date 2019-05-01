clear

eqns={
  'dV/dt=Iapp+@current';
};

 numcells = [1]
spec=[];
T0 = 2000;
spec.nodes(1).name = 'FSI';
spec.nodes(1).size = numcells;
spec.nodes(1).equations = eqns;
spec.nodes(1).mechanism_list = {'FSIGolombK','FSIGolombKdr','FSIInput','FSIGolombNa','FSILeak'}; 
spec.nodes(1).parameters = {'v_IC',-90, 'Tfinal', T0, 'Iapp',0};


vary={
  '(FSI)',			'tonic',	[0:10];
  '(FSI)',			'vl',	[-70:5:-55];
  '(FSI)',			'gd',	[3];
};

namearray = cellfun(@num2str,vary,'UniformOutput',0);
namestr = strjoin(reshape(namearray, 1, []));
cd '/projectnb/crc-nak/chartove/dynasim/'; %try to cd to this directory and leave data_dir blank
cluster_flag = 1;
overwrite_flag = 1;
save_data_flag = 1;
save_results_flag = 1;
verbose_flag = 1;
compile_flag = 0;
disk_flag = 0;
downsample_factor = 10;

[~,~]=dsSimulate(spec,...
              'analysis_functions', {@gvFRsoma, @gvCalcSpikePower},...
              'save_data_flag',save_data_flag,'study_dir','single_compartment_leak_reversal_low_gd_v2',...
              'cluster_flag',cluster_flag,'verbose_flag',verbose_flag,...
              'overwrite_flag',overwrite_flag,'tspan',[0 T0],...
              'save_results_flag',save_results_flag,'solver','rk4',...
              'compile_flag',compile_flag,...
			  'copy_run_file_flag',1, 'copy_mech_files_flag',1, ...
              'disk_flag',disk_flag,'downsample_factor',downsample_factor,...
              'vary',vary, 'dt', .01, ...
				  'plot_functions',{@dsPlot,@dsPlot,@dsPlot,@dsPlot},...
              'plot_options',{{'plot_type','waveform','format','png'},...
							{'plot_type','rastergram','format','png'},...
							{'plot_type','density','format','png'},...
                              {'plot_type','power','format','png',...
                               'freq_limits',[0 100]}});