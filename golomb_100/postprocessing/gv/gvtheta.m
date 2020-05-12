function theta = gvtheta(data,varargin)
    y_array = dsImportResults(pwd,@gvCalcPower);
    y = cell2mat(y_array(data.simulator_options.sim_id));
    theta = {sum(y((3<(y(:,2)) & (y(:,2))<7),1))};
end