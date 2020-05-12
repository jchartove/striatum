function fr_D1 = gvFRD1(data,varargin)

	fr = dsCalcFR(data,'bin_size',1, 'bin_shift',1);
	fr_D1 = {fr.D1_V_FR(1)};
end