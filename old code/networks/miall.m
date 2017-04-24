function [mitable,mipeaks,mipeaks2,milocs,milocs2,mibool,miinhib,nmitable,nmipeaks,nmipeaks2,nmilocs,nmilocs2,nmibool,nmiinhib,vitable,vipeaks,vipeaks2,vilocs,vilocs2,vibool,viinhib,cetable,cepeaks,cepeaks2,celocs,celocs2,cebool,ceinhib] = miall(traces)
n = size(traces,1);
t = size(traces,2);
fileID = fopen('miall.txt','w');

mitable = zeros(n,n,51);
mipeaks = zeros(n,n);
milocs = zeros(n,n);
mipeaks2 = zeros(n,n);
milocs2 = zeros(n,n);

nmitable = zeros(n,n,51);
nmipeaks = zeros(n,n);
nmilocs = zeros(n,n);
nmipeaks2 = zeros(n,n);
nmilocs2 = zeros(n,n);

vitable = zeros(n,n,51);
vipeaks = zeros(n,n);
vilocs = zeros(n,n);
vipeaks2 = zeros(n,n);
vilocs2 = zeros(n,n);

cetable = zeros(n,n,51);
cepeaks = zeros(n,n);
celocs = zeros(n,n);
cepeaks2 = zeros(n,n);
celocs2 = zeros(n,n);
for i = 1:n
    for j = 1:n
        i
        j
        if i ~= j
            for k = 1:50 %maximum allowable lag is 10s
                [mitable(i,j,k),nmitable(i,j,k),vitable(i,j,k),cetable(i,j,k)] = mnc(traces(j,k:end),traces(i,1:(t-k)+1)); 
                fprintf(fileID,'%d, %d, %d, %d, %d, %d, %d \r\n',i,j,k,mitable(i,j,k),nmitable(i,j,k),vitable(i,j,k),cetable(i,j,k));
            end
        end
        [maxcorr,maxlag] = max(mitable(i,j,:));
        [mincorr,minlag] = min(mitable(i,j,:));
        mipeaks(i,j) = maxcorr;
        mipeaks2(i,j) = mincorr;
        milocs(i,j) = maxlag;
        milocs2(i,j) = minlag;
		
		[maxcorr,maxlag] = max(nmitable(i,j,:));
        [mincorr,minlag] = min(nmitable(i,j,:));
        nmipeaks(i,j) = maxcorr;
        nmipeaks2(i,j) = mincorr;
        nmilocs(i,j) = maxlag;
        nmilocs2(i,j) = minlag;
		
		[maxcorr,maxlag] = max(vitable(i,j,:));
        [mincorr,minlag] = min(vitable(i,j,:));
        vipeaks(i,j) = maxcorr;
        vipeaks2(i,j) = mincorr;
        vilocs(i,j) = maxlag;
        vilocs2(i,j) = minlag;
		
		[maxcorr,maxlag] = max(vitable(i,j,:));
        [mincorr,minlag] = min(vitable(i,j,:));
        vipeaks(i,j) = maxcorr;
        vipeaks2(i,j) = mincorr;
        vilocs(i,j) = maxlag;
        vilocs2(i,j) = minlag;
    end
end
linear = reshape(mipeaks,[1,numel(mipeaks)]);
mibool = (mipeaks > mean(linear) + 1.43*std(linear));

linear2 = reshape(mipeaks2,[1,numel(mipeaks2)]);
bool2 = (mipeaks2 < mean(linear2) - 1.43*std(linear2));
miinhib = (bool2-mibool)>0;

linear = reshape(nmipeaks,[1,numel(nmipeaks)]);
nmibool = (nmipeaks > mean(linear) + 1.43*std(linear));

linear2 = reshape(nmipeaks2,[1,numel(nmipeaks2)]);
bool2 = (nmipeaks2 < mean(linear2) - 1.43*std(linear2));
nmiinhib = (bool2-nmibool)>0;

linear = reshape(vipeaks,[1,numel(vipeaks)]);
vibool = (vipeaks > mean(linear) + 1.43*std(linear));

linear2 = reshape(vipeaks2,[1,numel(vipeaks2)]);
bool2 = (vipeaks2 < mean(linear2) - 1.43*std(linear2));
viinhib = (bool2-vibool)>0;

linear = reshape(cepeaks,[1,numel(cepeaks)]);
cebool = (cepeaks > mean(linear) + 1.43*std(linear));

linear2 = reshape(cepeaks2,[1,numel(cepeaks2)]);
bool2 = (cepeaks2 < mean(linear2) - 1.43*std(linear2));
ceinhib = (bool2-cebool)>0;

fclose('all');
end

