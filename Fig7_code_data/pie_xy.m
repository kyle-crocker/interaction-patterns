function h = pie_xy(X,xpos,ypos,scale)
%This is a very narrowly designed function which plots takes a three
%element vector and turns it into a pie chart which gets plotted on an x-y
%axis. The size is controlled using the scale variable. The first fraction
%is blue, the second red, and the third grey
X = X + 1e-5;
X = X/sum(X);
h = pie(X,repmat({''},size(X,2)));

colormap([0 0 1;      %// blue
          1 0 0;      %// red
          .5 .5 .5])  %// grey

%Get aspect of current figure
xlimits = xlim;
ylimits = ylim;

%Move to xpos,ypos and rescale by scale
for k = 1:length(h) % Walk the vector of patch handles
    if strcmp(get(h(k),'Type'),'patch') % Patch graphics
        XData = get(h(k),'XData');  % Extract current data
        YData = get(h(k),'YData');
        set(h(k),'XData',XData*scale + xpos); % Insert modified data
        set(h(k),'YData',YData*scale + ypos);
        set(h(k),'EdgeColor','none');
    end
end

end

