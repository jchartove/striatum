for n = 1:length(data)
    %h= dsPlot(data(n),'visible', 'off');
	dsAnalyze(data(n),{@gvfr},'result_file',['power_sim' num2str(n)]);
    %print(h,['foo_sim' num2str(n)],'-dpng')
end

function data = gvfr(data, varargin)
	data = dsCalcFR(data,'bin_size',1, 'bin_shift',1);
	data = data.result.soma_V_FR(1);
end

function posthoc(data)
for n = 1:121
    %h= dsPlot(data(n),'visible', 'off');
	dsAnalyze(data(n),{@gvfr},'result_file',['study_sim' num2str(n) '_analysis1_gvfr']);
    %studyinfo.simulations(n).result_functions = {@gvfr};
    %studyinfo.simulations(n).result_files = {['study_sim' num2str(n) '_analysis1_gvfr.mat']};
    %print(h,['foo_sim' num2str(n)],'-dpng')
end
end
%save('studyinfo2.mat','studyinfo')

function outdata = gvfr(data, varargin)
	outdata = dsCalcFR(data,'bin_size',1, 'bin_shift',1);
	outdata = {outdata.soma_V_FR(1)};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%command history

gv.Run
load('study_sim1_data.mat')
data = dsImport(pwd)
n = 1
dsAnalyze(data(n),{@dsCalcFR},'result_file',['power_sim' num2str(n)],'function_options',{{ 'bin_size',1, 'bin_shift',1}});
load('power_sim1.mat')
close all
posthoc(data)
dbstop if error
posthoc(data)
dbstop
posthoc(data)
result = evalFnWithArgs(fInd, data, func, options, varargin{:});
posthoc
posthoc(data)
gv.Run
load('studyinfo.mat')
studyinfo.simulations.result_functions = {@gvfr}
for n = 1:length(data)
%h= dsPlot(data(n),'visible', 'off');
%dsAnalyze(data(n),{@gvfr},'result_file',['study_sim' num2str(n) '_analysis1_gvfr']);
studyinfo.simulations(n).result_functions = {@gvfr};
studyinfo.simulations(n).result_files = {['study_sim' num2str(n) '_analysis1_gvfr.mat']};
%print(h,['foo_sim' num2str(n)],'-dpng')
end
save('studyinfo2.mat','studyinfo')
dsImportResults('studyinfo.mat',@gvfr)
a = dsImportResults('studyinfo.mat',@gvfr)
gv.Run
clear
gv.Run
load('old.mat')
n = 2
studyinfo.simulations(n).result_functions
for n = 1:length(data)
%h= dsPlot(data(n),'visible', 'off');
%dsAnalyze(data(n),{@gvfr},'result_file',['study_sim' num2str(n) '_analysis1_gvfr']);
studyinfo.simulations(n).result_functions = {@gvfr};
studyinfo.simulations(n).result_files = {['study_sim' num2str(n) '_analysis1_gvfr.mat']};
%print(h,['foo_sim' num2str(n)],'-dpng')
end
for n = 1:121
%h= dsPlot(data(n),'visible', 'off');
%dsAnalyze(data(n),{@gvfr},'result_file',['study_sim' num2str(n) '_analysis1_gvfr']);
studyinfo.simulations(n).result_functions = {@gvfr};
studyinfo.simulations(n).result_files = {['study_sim' num2str(n) '_analysis1_gvfr.mat']};
%print(h,['foo_sim' num2str(n)],'-dpng')
end
save('studyinfo2.mat','studyinfo')
gv.Run
clear
data = dsImport(pwd)
posthoc(data)
gv.Run