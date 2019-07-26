function [b,pgon] = errorghost(data,xaxis,color)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
top = mean(data)+std(data);
bot = mean(data)-std(data);
pgon = polyshape([xaxis fliplr(xaxis)],[top fliplr(bot)],'Simplify', false);
b = plot(pgon,'HandleVisibility','off', 'FaceColor',color);
hold on
b.EdgeAlpha = 0;
b.FaceAlpha = .15;
end

