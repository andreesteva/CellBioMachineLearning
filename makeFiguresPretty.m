figHandles = get(0,'Children');

for i =1:length(figHandles)
    figureHandle = figHandles(i);
    figureAxes = findall(figureHandle,'type','axes');
    % make all text in the figure to size 14 and bold
    set(findall(figureHandle,'type','text'),'fontSize',60,'fontWeight','bold') % Set text size
    set(get(gca,'xlabel'), 'FontSize', 20); % Set xlabel and ylabel fonts
    set(get(gca,'ylabel'), 'FontSize', 20); % Set xlabel and ylabel fonts
    set(figureAxes, 'fontsize', 15, 'linewidth', 2);
    
end