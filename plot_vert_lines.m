function plot_vert_lines(x,y_min,y_max)
% PLOT_VERT_LINES
%
% for each value x plot vertical lines extending from y_min to y_max
% If x is a matrix, each row in x is plotted in a different color
colorinc=1/(size(x,1)+1);
coloraccum=colorinc;
for m=1:size(x,1)
    for n=1:size(x,2)
        line([x(m,n),x(m,n)],[y_min,y_max],'linewidth',2,'color',[coloraccum,coloraccum,coloraccum]);
    end
    coloraccum+=colorinc;
end
