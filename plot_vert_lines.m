function plot_vert_lines(x,y_min,y_max)
% function plot_vert_lines(x,y_min,y_max)
%
% for each value x plot vertical lines extending from y_min to y_max
% If x is a matrix, each row in x is plotted in a different color
colors=['k','r','g','b','m','c','w'];
clr_k=1;
for m=1:size(x,1)
    for n=1:size(x,2)
        line([x(m,n),x(m,n)],[y_min,y_max],'linewidth',1,'color',colors(clr_k));
    end
    clr_k+=1;
    if (clr_k > length(colors))
        clr_k=1;
    end
end
