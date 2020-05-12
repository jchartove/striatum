function gv_PLV_D1 = gvPLVD1(data,varargin)
	%[sim,~] = dsImport(data);
	stim = decimate(data.model.fixed_variables.dend_iPeriodicPulsesBen_s2(:,1),10);
foi = data.model.parameters.dend_iPeriodicPulsesBen_PPfreq;
[gv_PLV, gv_PLV_D1] = PLV_F(data.D1_V, stim, foi, data.time);
end