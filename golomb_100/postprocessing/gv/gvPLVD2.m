function gv_PLV_D2 = gvPLVD2(data,varargin)
	%[sim,~] = dsImport(data);
	stim = decimate(data.model.fixed_variables.dend_iPeriodicPulsesBen_s2(:,1),10);
foi = data.model.parameters.dend_iPeriodicPulsesBen_PPfreq;
[gv_PLV, gv_PLV_D2] = PLV_F(data.D2_V, stim, foi, data.time);
end