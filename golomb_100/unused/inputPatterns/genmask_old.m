function mask = genmask(Npre,Npost,con,cond,dir)
	mask = rand(Npost,Npre)<con
	if dir
		mask = mask - diag(diag(mask));
		filename = strcat('synmask_', mat2str(clock),'.mat');
	else
		mask = triu(mask);
		mask = mask + mask.';
		mask = mask - diag(diag(mask));
		filename = strcat('gjmask_', mat2str(clock),'.mat');
	end
	save(filename)
	pwd
end