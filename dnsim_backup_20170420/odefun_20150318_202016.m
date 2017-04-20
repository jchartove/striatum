function [Y,T] = odefun_20150318_202016
tspan=[0 1000]; dt=0.01;
T=tspan(1):dt:tspan(2); nstep=length(T);
fileID = 1; nreports = 5; enableLog = 1:(nstep-1)/nreports:nstep;enableLog(1) = [];
fprintf('\nSimulation interval: %g-%g\n',tspan);fprintf('Starting integration (euler, dt=%g)\n',dt);FSI_FSI_iinh_FSI_width = inf;
FSI_FSI_iinh_FSI_Nmax = max((10),(10));
FSI_FSI_iinh_FSI_srcpos = linspace(1,FSI_FSI_iinh_FSI_Nmax,(10))'*ones(1,(10));
FSI_FSI_iinh_FSI_dstpos = (linspace(1,FSI_FSI_iinh_FSI_Nmax,(10))'*ones(1,(10)))';
FSI_FSI_iinh_FSI_netcon = rand((10),(10))<0.3;;
FSI_FSI_com_UB = max((10),(10));
FSI_FSI_com_Xpre = linspace(1,FSI_FSI_com_UB,(10))'*ones(1,(10));
FSI_FSI_com_Xpost = (linspace(1,FSI_FSI_com_UB,(10))'*ones(1,(10)))';
FSI_FSI_com_mask = abs(FSI_FSI_com_Xpre-FSI_FSI_com_Xpost)<=(0.5);
FSI_FSI_iGAP_UB = max((10),(10));
FSI_FSI_iGAP_Xpre = linspace(1,FSI_FSI_iGAP_UB,(10))'*ones(1,(10));
FSI_FSI_iGAP_Xpost = (linspace(1,FSI_FSI_iGAP_UB,(10))'*ones(1,(10)))';
FSI_FSI_iGAP_mask = rand((10),(10))<0.33;;
FSI_iK_FSI_aN = @(IN) (.1-.01.*(IN+65))./(exp(1-.1.*(IN+65))-1);
FSI_iK_FSI_bN = @(IN) .125.*exp(-(IN+65)./80);
FSI_iK_FSI_IK_FSI = @(IN,FSI_iK_FSI_nK) (36).*FSI_iK_FSI_nK.^4.*((IN+65)-(-12));
FSI_iLeak_FSI_ILeak_FSI = @(IN) (0.3).*((IN+65)-(10.6));
FSI_iNa_FSI_aM = @(IN) (2.5-.1.*(IN+65))./(exp(2.5-.1.*(IN+65))-1);
FSI_iNa_FSI_bM = @(IN) 4.*exp(-(IN+65)./18);
FSI_iNa_FSI_aH = @(IN) .07.*exp(-(IN+65)./20);
FSI_iNa_FSI_bH = @(IN) 1./(exp(3-.1.*(IN+65))+1);
FSI_iNa_FSI_INa_FSI = @(IN,FSI_iNa_FSI_mNa,FSI_iNa_FSI_hNa) (120).*FSI_iNa_FSI_mNa.^3.*FSI_iNa_FSI_hNa.*((IN+65)-(115));
FSI_FSI_iinh_FSI_Iinh = @(IN,FSI_FSI_iinh_FSI_s) ((2.5).*(FSI_FSI_iinh_FSI_netcon*FSI_FSI_iinh_FSI_s).*(IN-(-80)));
FSI_FSI_com_ICOM = @(V1,V2) (0.2).*sum(((V1*ones(1,size(V1,1)))'-(V2*ones(1,size(V2,1)))).*FSI_FSI_com_mask,2);
FSI_FSI_iGAP_IGAP = @(V1,V2) (0.2).*sum(((V1*ones(1,size(V1,1)))'-(V2*ones(1,size(V2,1)))).*FSI_FSI_iGAP_mask,2);
FSI_V = zeros(10,nstep);
FSI_V(:,1) = [0  0  0  0  0  0  0  0  0  0];
FSI_iK_FSI_nK = zeros(10,nstep);
FSI_iK_FSI_nK(:,1) = [0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1];
FSI_iNa_FSI_mNa = zeros(10,nstep);
FSI_iNa_FSI_mNa(:,1) = [0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1];
FSI_iNa_FSI_hNa = zeros(10,nstep);
FSI_iNa_FSI_hNa(:,1) = [0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1];
FSI_FSI_iinh_FSI_s = zeros(10,nstep);
FSI_FSI_iinh_FSI_s(:,1) = [0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1         0.1];
tstart = tic;
for k=2:nstep
  t=T(k-1);
  F=((-((36).*FSI_iK_FSI_nK(:,k-1).^4.*((FSI_V(:,k-1)+65)-(-12))))+((-((0.3).*((FSI_V(:,k-1)+65)-(10.6))))+((-((120).*FSI_iNa_FSI_mNa(:,k-1).^3.*FSI_iNa_FSI_hNa(:,k-1).*((FSI_V(:,k-1)+65)-(115))))+((-(((2.5).*(FSI_FSI_iinh_FSI_netcon*FSI_FSI_iinh_FSI_s(:,k-1)).*(FSI_V(:,k-1)-(-80)))))+((((0.2).*sum(((FSI_V(:,k-1)*ones(1,size(FSI_V(:,k-1),1)))'-(FSI_V(:,k-1)*ones(1,size(FSI_V(:,k-1),1)))).*FSI_FSI_com_mask,2)))+((((0.2).*sum(((FSI_V(:,k-1)*ones(1,size(FSI_V(:,k-1),1)))'-(FSI_V(:,k-1)*ones(1,size(FSI_V(:,k-1),1)))).*FSI_FSI_iGAP_mask,2)))+0))))))/c;
  FSI_V(:,k) = FSI_V(:,k-1) + dt*F;
  F=((.1-.01.*(FSI_V(:,k-1)+65))./(exp(1-.1.*(FSI_V(:,k-1)+65))-1)).*(1-FSI_iK_FSI_nK(:,k-1))-(.125.*exp(-(FSI_V(:,k-1)+65)./80)).*FSI_iK_FSI_nK(:,k-1);
  FSI_iK_FSI_nK(:,k) = FSI_iK_FSI_nK(:,k-1) + dt*F;
  F=((2.5-.1.*(FSI_V(:,k-1)+65))./(exp(2.5-.1.*(FSI_V(:,k-1)+65))-1)).*(1-FSI_iNa_FSI_mNa(:,k-1))-(4.*exp(-(FSI_V(:,k-1)+65)./18)).*FSI_iNa_FSI_mNa(:,k-1);
  FSI_iNa_FSI_mNa(:,k) = FSI_iNa_FSI_mNa(:,k-1) + dt*F;
  F=(.07.*exp(-(FSI_V(:,k-1)+65)./20)).*(1-FSI_iNa_FSI_hNa(:,k-1))-(1./(exp(3-.1.*(FSI_V(:,k-1)+65))+1)).*FSI_iNa_FSI_hNa(:,k-1);
  FSI_iNa_FSI_hNa(:,k) = FSI_iNa_FSI_hNa(:,k-1) + dt*F;
  F=(1 + tanh(FSI_V(:,k-1)./10)./2).*(1-FSI_FSI_iinh_FSI_s(:,k-1))./0.5 - FSI_FSI_iinh_FSI_s(:,k-1)./(12);;
  FSI_FSI_iinh_FSI_s(:,k) = FSI_FSI_iinh_FSI_s(:,k-1) + dt*F;
  if any(k == enableLog)
    elapsedTime = toc(tstart);
    elapsedTimeMinutes = floor(elapsedTime/60);
    elapsedTimeSeconds = rem(elapsedTime,60);
    if elapsedTimeMinutes
        fprintf(fileID,'Processed %g of %g ms (elapsed time: %g m %.3f s)\n',T(k),T(end),elapsedTimeMinutes,elapsedTimeSeconds);
    else
        fprintf(fileID,'Processed %g of %g ms (elapsed time: %.3f s)\n',T(k),T(end),elapsedTimeSeconds);
    end
  end
end
Y=cat(1,FSI_V,FSI_iK_FSI_nK,FSI_iNa_FSI_mNa,FSI_iNa_FSI_hNa,FSI_FSI_iinh_FSI_s)';
