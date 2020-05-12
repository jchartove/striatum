function fr_D2 = gvFRD2(data,varargin)

	fr = dsCalcFR(data,'bin_size',1, 'bin_shift',1);
	fr_D2 = {fr.D2_V_FR(1)};
end