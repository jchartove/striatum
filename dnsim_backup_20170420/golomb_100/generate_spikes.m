function [avgfr,spike_pairs, spike_indicator] = generate_spikes(data, v_new, filenew, time, T_start, dt, numcells)
        T_new = length(time)-T_start-1; %this is dumb but saves me time rewriting
        if numcells > 1
            lfp = mean(data');
        else
            lfp = data';
        end

        spike_indicator = zeros(numcells,T_new+1); 
        synch_indicator = zeros(numcells, numcells, T_new);

        spikes = zeros(1,numcells);

        for t = 1:T_new
            spike_indicator(:,t) = (v_new(t,:)<0) & (v_new(t+1,:) >= 0);
            s = (v_new(t,:)<0) & (v_new(t+1,:) >= 0);
            spikes = spikes + s;
        end

        avgfr = mean(spikes)/(T_new/(1000/dt));

        spike_pairs = 0;
%         if numcells<=10 % this calculation is so expensive i'm sorry
% 		
%             %rectangle for convolutions
%             rect = 3*ones(1,5/dt);
%             rect = rect(rect > eps);
%             rect = [zeros(1,length(rect)) rect];
% 
%             wide_spikes = zeros(numcells,T_new+1);
%             for c = 1:numcells
%                 wide_spikes(c,:) = conv(spike_indicator(c,:),rect,'same');
%             end
% 
%             for b = 1:numcells
%                 for c = 1:numcells
%                     if b ~= c
%                         wide_sum = wide_spikes(b,:) + wide_spikes(c,:);
%                         foo = wide_sum > 5;
%                         wide_synch = diff(foo, [], 2);
%                         synch_indicator(b,c,:) = wide_synch == 1;
%                         spike_pairs = spike_pairs + sum(synch_indicator(b,c,:));
%                         c
%                     end
%                 end
%             end
%         end

        %%%%%%%%%%%%%%%%%%%% spike plots
        handle0 = figure;
        hist(spikes/(T_new/(1000/dt)))
        xlabel('Firing rate');
        imgtitle = strcat(filenew,'hist.png')
        title(imgtitle);
        saveas(handle0, imgtitle, 'png');

        handle1 = figure;
        plot(time,data);
        hold on;
        plot(time, lfp,'k','LineWidth',2);
        hold off;
        xlabel('Time');
        ylabel('Voltage');
        imgtitle = strcat(filenew,'spikes.png')
        title(imgtitle);
        saveas(handle1, imgtitle, 'png');

        xlim([T_start*dt+100 (T_start*dt)+200]);
        imgtitle = strcat(filenew,'spikes_zoom.png')
        title(imgtitle);
        saveas(handle1, imgtitle, 'png');

        %%%%%%%%%%%%%%%%%%%% rasters
        handle2 = figure;
        imagesc(data');
		xlim([T_start length(time)])
        colorbar;
        xlabel('Time');
        ylabel('Cell number');
        imgtitle = strcat(filenew,'raster.png')
        title(imgtitle);
        saveas(handle2, imgtitle, 'png');

        xlim([T_start+(100/dt) T_start+(200/dt)]);
        imgtitle = strcat(filenew,'raster_zoom.png')
        title(imgtitle);
        saveas(handle2, imgtitle, 'png');
end