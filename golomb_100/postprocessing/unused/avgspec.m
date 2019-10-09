function avgspectrum = avgspec(directory)
	%matchstr example: '*gcon0_*icon0-*numcells10_*spectrum.mat'
    %matchstr example: '*fractionshared0_*gGAP0_*spectrum.mat'
    %matchstr example: '*fractiongamma0-*fractionshared0_*spectrum.mat'
    for h = (1)
        for n = (0.5:0.5:2)
                matchstr = strcat('*taumult',num2str(n,'%0.1f'),'*numcells1'); %'*taub',num2str(m,'%03d'),
                matchstr = strrep(matchstr,'.','pt');
                matchstr = strcat(matchstr,'_*spectrum.mat')
                cd(directory);
                datadir = [directory, matchstr];
                datafiles = dir(datadir);
                
             
                avgspectrum = zeros(1,151); %this should not be hardcoded
                for file = datafiles'
                    load(file.name, 'y')
                    newname = strrep(file.name,'-','_');
                    newname = strrep(newname,'.','_');
                    newname = strrep(newname,'+','_');
                    S.(newname) = y;
                    avgspectrum = avgspectrum + y;
                    clearvars y;
                end
                avgspectrum = avgspectrum/length(datafiles);
                filename = strrep(matchstr, '*', '_');
                filename = strcat(filename,'avg.mat');

                handle3 = figure;
                plot(avgspectrum);
                ylabel('Power');
                xlabel('Freq [Hz]')
                xlim([0 100])
                set(gca,'XTick',[0:5:100]);
                imgtitle = strcat(filename,'.png');
                title(imgtitle);

                saveas(handle3, imgtitle, 'png');

                s = cell2mat(struct2cell(S));
                handle4 = figure;
                imagesc(s);
                ylabel('taub (D-current)')
                xlabel('Freq [Hz]');
                colorbar;
                xlim([0 100])
                set(gca,'XTick',[0:5:100]);
                imgtitle = strcat(filename,'color.png');
                title(imgtitle);

                saveas(handle4, imgtitle, 'png');

                save(filename);

                fclose('all');
                clearvars avgspectrum;
                clearvars S;
                clearvars s;
                close all;
            end
    end
end