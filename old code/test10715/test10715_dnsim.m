% Model: test10715
cd /project/crc-nak/jchartove/test10715;
spec=[];
spec.nodes(1).label = 'FSI';
spec.nodes(1).multiplicity = 10;
spec.nodes(1).dynamics = {'V''=current/c'};
spec.nodes(1).mechanisms = {'FSI_iK_FSI','FSI_iLeak_FSI','FSI_iNa_FSI'};
spec.nodes(1).parameters = {'v',50,'c',1,'v(0)',50,'V',50};
spec.connections(1,1).label = 'FSI-FSI';
spec.connections(1,1).mechanisms = {'FSI_FSI_iinh_FSI','FSI_FSI_com','FSI_FSI_iGAP'};
spec.connections(1,1).parameters = [];
%dnsim(spec); % open model in DNSim GUI

% DNSim simulation and plots:
data = runsim(spec,'timelimits',[0 100],'dt',.02,'SOLVER','euler'); % simulate DNSim models
plotv(data,spec,'varlabel','V'); % quickly plot select variables
%visualizer(data); % ugly interactive tool hijacked to visualize sim_data

% Sweep over parameter values:
model=buildmodel(spec); % parse DNSim spec structure
simstudy(model,{'FSI'},{'N'},{'[1 2]'},'timelimits',[0 100],'dt',.02,'SOLVER','euler'); % N = # of cells


% Manual simulation and plots:
%{
%-----------------------------------------------------------
% Auxiliary variables:
	FSI_FSI_iinh_FSI_width = inf;
	FSI_FSI_iinh_FSI_Nmax = max((10),(10));
	FSI_FSI_iinh_FSI_srcpos = linspace(1,FSI_FSI_iinh_FSI_Nmax,(10))'*ones(1,(10));
	FSI_FSI_iinh_FSI_dstpos = (linspace(1,FSI_FSI_iinh_FSI_Nmax,(10))'*ones(1,(10)))';
	FSI_FSI_iinh_FSI_netcon = rand((10),(10))<0.33;;
	FSI_FSI_com_UB       = max((10),(10));
	FSI_FSI_com_Xpre     = linspace(1,FSI_FSI_com_UB,(10))'*ones(1,(10));
	FSI_FSI_com_Xpost    = (linspace(1,FSI_FSI_com_UB,(10))'*ones(1,(10)))';
	FSI_FSI_com_mask     = abs(FSI_FSI_com_Xpre-FSI_FSI_com_Xpost)<=(0.5);
	FSI_FSI_iGAP_UB      = max((10),(10));
	FSI_FSI_iGAP_Xpre    = linspace(1,FSI_FSI_iGAP_UB,(10))'*ones(1,(10));
	FSI_FSI_iGAP_Xpost   = (linspace(1,FSI_FSI_iGAP_UB,(10))'*ones(1,(10)))';
	FSI_FSI_iGAP_mask    = rand((10),(10))<0.33;;

% Anonymous functions:
	FSI_iK_FSI_aN        = @(IN) (.1-.01.*(IN+65))./(exp(1-.1.*(IN+65))-1);
	FSI_iK_FSI_bN        = @(IN) .125.*exp(-(IN+65)./80);          
	FSI_iK_FSI_IK_FSI    = @(IN,FSI_iK_FSI_nK) (36).*FSI_iK_FSI_nK.^4.*((IN+65)-(-12));
	FSI_iLeak_FSI_ILeak_FSI = @(IN) (0.3).*((IN+65)-(10.6));          
	FSI_iNa_FSI_aM       = @(IN) (2.5-.1.*(IN+65))./(exp(2.5-.1.*(IN+65))-1);
	FSI_iNa_FSI_bM       = @(IN) 4.*exp(-(IN+65)./18);             
	FSI_iNa_FSI_aH       = @(IN) .07.*exp(-(IN+65)./20);           
	FSI_iNa_FSI_bH       = @(IN) 1./(exp(3-.1.*(IN+65))+1);        
	FSI_iNa_FSI_INa_FSI  = @(IN,FSI_iNa_FSI_mNa,FSI_iNa_FSI_hNa) (120).*FSI_iNa_FSI_mNa.^3.*FSI_iNa_FSI_hNa.*((IN+65)-(115));
	FSI_FSI_iinh_FSI_Iinh = @(IN,FSI_FSI_iinh_FSI_s) ((2.5).*(FSI_FSI_iinh_FSI_netcon*FSI_FSI_iinh_FSI_s).*(IN-(-80)));
	FSI_FSI_com_ICOM     = @(V1,V2) (0.2).*sum(((V1*ones(1,size(V1,1)))'-(V2*ones(1,size(V2,1)))).*FSI_FSI_com_mask,2);
	FSI_FSI_iGAP_IGAP    = @(V1,V2) (0.2).*sum(((V1*ones(1,size(V1,1)))'-(V2*ones(1,size(V2,1)))).*FSI_FSI_iGAP_mask,2);

% ODE Handle, ICs, integration, and plotting:
ODEFUN = @(t,X) [((-((36).*X(11:20).^4.*((X(1:10)+65)-(-12))))+((-((0.3).*((X(1:10)+65)-(10.6))))+((-((120).*X(21:30).^3.*X(31:40).*((X(1:10)+65)-(115))))+((-(((2.5).*(FSI_FSI_iinh_FSI_netcon*X(41:50)).*(X(1:10)-(-80)))))+((((0.2).*sum(((X(1:10)*ones(1,size(X(1:10),1)))'-(X(1:10)*ones(1,size(X(1:10),1)))).*FSI_FSI_com_mask,2)))+((((0.2).*sum(((X(1:10)*ones(1,size(X(1:10),1)))'-(X(1:10)*ones(1,size(X(1:10),1)))).*FSI_FSI_iGAP_mask,2)))+0))))))/(1);((.1-.01.*(X(1:10)+65))./(exp(1-.1.*(X(1:10)+65))-1)).*(1-X(11:20))-(.125.*exp(-(X(1:10)+65)./80)).*X(11:20);((2.5-.1.*(X(1:10)+65))./(exp(2.5-.1.*(X(1:10)+65))-1)).*(1-X(21:30))-(4.*exp(-(X(1:10)+65)./18)).*X(21:30);(.07.*exp(-(X(1:10)+65)./20)).*(1-X(31:40))-(1./(exp(3-.1.*(X(1:10)+65))+1)).*X(31:40);(1 + tanh(X(1:10)./10)./2).*(1-X(41:50))./0.5 - X(41:50)./(12);];
IC = [0           0           0           0           0           0           0           0           0           0         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1];

[t,y]=ode23(ODEFUN,[0 100],IC);   % numerical integration
figure; plot(t,y);           % plot all variables/functions
try legend('FSI\_V','FSI\_iK\_FSI\_nK','FSI\_iNa\_FSI\_mNa','FSI\_iNa\_FSI\_hNa','FSI\_FSI\_iinh\_FSI\_s'); end
%-
%}
