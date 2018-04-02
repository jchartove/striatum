function posthoc_gvbands 
	dsAnalyze(pwd, ...
    'analysis_functions', {@get_fft,@gvtheta,@gvdelta},...
    'analysis_options', {},...
    'plot_functions',{@dsPlot, @dsPlot}, 'plot_options',...
    {
        {'plot_type','power','variable','soma_V'},...
        {'plot_type','power','variable','D1_V'}...
    },...
    'load_all_data_flag',1, 'verbose_flag',1, 'parfor_flag',1, 'save_results_flag',1,'close_fig_flag',1);
	%gv.Run
end

%seperate function for every variable including normalized
%move some preprocessing to main function
%also functions for MSNs

function outdata = gvtheta(data,varargin)
	%y = get_fft(data);
    %simID = 1;
    dsImportRsults(pwd,@get_fft)
    outdata = {sum(y(4:7))};
end

function outdata = gvdelta(data,varargin)
	%y = get_fft(data);
	outdata = sum(y(1:3));
end

 %   totalp = sum(y(1:150)); %total power. below: eeg bands
%	ap = sum(y(8:12));
    %for the broader peaks, also find peak location
%    [~,lowpeak] = max(y(1:12));
%    bp = sum(y(13:35));
%    [~,bpeak] = max(y(13:35));
%    bpeak = bpeak + 12;
%    gplow = sum(y(36:65));
%    [~,glopeak] = max(y(36:65));
%    glopeak = glopeak + 35;
%    gphigh = sum(y(66:100));
%    [~,ghipeak] = max(y(66:100));
%    ghipeak = ghipeak + 65;
%    hfop = sum(y(101:150));
%    [~,hfopeak] = max(y(101:150));
%    hfopeak = hfopeak + 100;
%    [~,gpeak] = max(y(36:100));
%    gpeak = gpeak + 35;
%    [~,hipeak] = max(y(66:150));
%    hipeak = hipeak + 65;