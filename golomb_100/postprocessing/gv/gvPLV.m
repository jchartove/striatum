function gv_PLV = gvPLV(data,varargin)
	%[sim,~] = dsImport(data);
	stim = decimate(data.model.fixed_variables.dend_iPeriodicPulsesBen_s2(:,1),10);
foi = data.model.parameters.dend_iPeriodicPulsesBen_PPfreq;
[gv_PLV, gv_PLV_s] = PLV_F(data.soma_V, stim, foi, data.time);
end