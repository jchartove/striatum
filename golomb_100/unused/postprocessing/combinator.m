function combinator(directory)
    %opengl('save', 'software')
	%matchstr example: '*gcon0_*icon0-*numcells10_*spectrum.mat'
    %matchstr example: '*fractionshared0_*gGAP0_*spectrum.mat'
    %matchstr example: '*rate',num2str(n,'%0.1f'),'*numcells',num2str(h)); %'*taub',num2str(m,'%03d'),
    %for n = 110:10:200
        for n = [1,2,5,10,50,80,100,120,130]
                datatype = 'gpe';
                matchstr = strcat('*',num2str(n),'Hz*times.mat');
                cd(directory);
                datadir = [directory, matchstr];
                datafiles = dir(datadir);

                for file = datafiles'
                    load(file.name, 'train_prestim', 'train_instim', 'train_poststim'); 
                    newname = strrep(file.name,'-','_');
                    newname = strrep(newname,'.','_');
                    newname = strrep(newname,'+','_');
                    S_pre.(newname) = train_prestim;
                    S_in.(newname) = train_instim;
                    S_post.(newname) = train_poststim;
                    clearvars train_prestim train_instim train_poststim;
                end
                filename = strrep(matchstr, '*', '_');
                filename = strcat(filename, datatype, '_combined.mat');

                %if i were smart i could do this without duplicate code
                s_pre = struct2cell(S_pre);
                pre_len=cellfun('length',s_pre);
                pre_max=max(pre_len);
                final_pre = zeros(length(s_pre),pre_max);
                for m = 1:length(s_pre)
                    row_m = cell2mat(s_pre(m));
                    long_m = padarray(row_m,[(pre_max-length(row_m)) 0],'post');
                    final_pre(m,:) = final_pre(m,:) + long_m';
                end
                [pre_cells,pre_times] = find(final_pre);
                handle1 = scatter(pre_times,pre_cells);
                imgtitle = strcat(filename,'_pre.fig');
                title(imgtitle);
                ylabel('Cells');
                xlabel('Time (ms)');
                saveas(handle1, imgtitle, 'fig');
                
                s_in = struct2cell(S_in);
                in_len=cellfun('length',s_in);
                in_max=max(in_len);
                final_in = zeros(length(s_in),in_max);
                for m = 1:length(s_in)
                    row_m = cell2mat(s_in(m));
                    long_m = padarray(row_m,[(in_max-length(row_m)) 0],'post');
                    final_in(m,:) = final_in(m,:) + long_m';
                end
                [in_cells,in_times] = find(final_in);
                handle2 = scatter(in_times,in_cells);
                imgtitle = strcat(filename,'_in.fig');
                title(imgtitle);
                ylabel('Cells');
                xlabel('Time (ms)')
                saveas(handle2, imgtitle, 'fig');
                
                s_post = struct2cell(S_post);
                post_len=cellfun('length',s_post);
                post_max=max(post_len);
                final_post = zeros(length(s_post),post_max);
                for m = 1:length(s_post)
                    row_m = cell2mat(s_post(m));
                    long_m = padarray(row_m,[(post_max-length(row_m)) 0],'post');
                    final_post(m,:) = final_post(m,:) + long_m';
                end
                [post_cells,post_times] = find(final_post);
                handle3 = scatter(post_times,post_cells);
                imgtitle = strcat(filename,'_post.fig');
                title(imgtitle);
                ylabel('Cells');
                xlabel('Time (ms)');
                saveas(handle3, imgtitle, 'fig');

                save(filename);

                fclose('all');
                clearvars S_pre S_in S_post;
                close all;
       end
    %end
end