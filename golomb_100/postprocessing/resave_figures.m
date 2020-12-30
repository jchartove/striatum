figure_names = dir('*.fig');

for f = 1:length(figure_names)

    open(figure_names(f).name)
    children = get(gcf, 'Children');
    
    xtick = children(end).XTick;
    try
        xticklabel = cellfun(@str2num, children(end).XTickLabel);
    end
    try
        xticklabel = str2num(children(end).XTickLabel);
    end
    
    ytick = children(end).YTick;
    try 
        yticklabel = cellfun(@str2num, children(end).YTickLabel);
    end
    try
        yticklabel = str2num(children(end).YTickLabel);
    end
    
    grandchildren = get(children(end), 'Children');
    cdata = grandchildren(end).CData;
    cdata(isnan(cdata)) = 0;
    figure_names(f).cdata = cdata;
    
    newfig = figure;
    h = pcolor(flipud(cdata));
    shading interp
    set(gca, 'XTick', xtick, 'XTickLabel', xticklabel, 'YTick', ytick, 'YTickLabel', sort(abs(yticklabel)))
    
    axis xy
    
    supersizeme(2)
    
    xlabel('Input frequency (Hz)');
    ylabel('Input strength (\mu A)')
    
    saveas(newfig, [figure_names(f).name, '_NEW.fig'])
    saveas(newfig, [figure_names(f).name, '_NEW.png'])
    
end

save('figs.mat');