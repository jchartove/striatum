function fr_soma = gvFRsoma(data,varargin)

	fr = dsCalcFR(data,'bin_size',1, 'bin_shift',1);
	fr_soma = {fr.soma_V_FR(1)};
    
%     T_total = size(soma_V,1)-1;
%     T_start = T_total*0.25;
%     new_T = T_total + 1; %removed a factor of 10. this is all very silly at this point
%     numcells = size(soma_V,2);
%       time = zeros(1,size(soma_V,1));
%     for j = 1:new_T
%         time(j) = (j-1)*simulator_options.dt;
%     end
% 
% 
% if datatype == 1
%             data = soma_V;
%         elseif datatype == 2
%             data = D1_V;
%             
% v_new = data(T_start:T_total+1,:);
%         T_new = length(time)-T_start-1; %this is dumb but saves me time rewriting
% 
%         spike_indicator = zeros(numcells,T_new+1); 
% 
%         spikes = zeros(1,numcells);
% 
%         for t = 1:T_new
%             spike_indicator(:,t) = (v_new(t,:)<0) & (v_new(t+1,:) >= 0);
%             s = (v_new(t,:)<0) & (v_new(t+1,:) >= 0);
%             spikes = spikes + s;
%         end
% 
%         avgfr = mean(spikes)/(T_new/(100/dt)); %was this 100? why was this 100?
end