function mask = genmask(Npre,Npost,con,cond,dir,aut,ko)
	mask = rand(Npost,Npre)<con;
	if not(dir)
		mask = triu(mask);
		mask = mask + mask.';
	end
	
	if aut
		mask = mask - diag(diag(mask));
	end
	
	if dir
		filename = strcat('synmask_', mat2str(clock),'.mat')
	else
		mask = mask - diag(diag(mask));
		filename = strcat('gjmask_', mat2str(clock),'.mat')
	end
	
	mask(:,(Npre-ko):end)=0
	save(strcat('/projectnb/crc-nak/chartove/dnsim/masks/', filename))
	pwd
end