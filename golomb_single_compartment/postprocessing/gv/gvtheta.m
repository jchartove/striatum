function outdata = gvtheta(data,varargin)
    y_array = dsImportResults(pwd,@get_fft);
    y = cell2mat(y_array(data.simulator_options.sim_id));
    outdata = {sum(y(4:7))};
end
