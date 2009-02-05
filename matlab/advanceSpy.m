function advanceSpy(matrix, pointsize, nContours, fig)
    % ADVANCESPY : matrix visualization, but shows some sense of magnitude and
    % also positivity
    % portions copied from spy.m
    % Seth R. Johnson
    
    if (nargin == 0)
        disp('You must enter a matrix name to spy')
        return
    end
    if (nargin < 2)
        pointsize = 5;
    end
    if (nargin < 3)
        nContours = 0;
    end
    if (nargin == 4) 
        figure(fig)
    end
    clf

    halfmaxval = max(max(matrix))/2 + .00000001;
    tinyval = max(max(matrix))/100;
    maxnegval = min(min(matrix));
    if (maxnegval >= 0), maxnegval = 0; end
    
    %plot positive really small vals
    [i,j] = find(matrix <= tinyval & matrix > 0);
    plotMatrix(i,j,pointsize,[0.75 0.75 0.75])
    hold on
    
    %plot positive small vals
    [i,j] = find(matrix <= halfmaxval & matrix > tinyval);
    plotMatrix(i,j,pointsize,[0.5 0.5 1.00])
    
    %plot positive large values
    [i,j] = find(matrix > halfmaxval);
    plotMatrix(i,j,pointsize,[0.0 0.0 1.00])
    
    if (maxnegval < 0)
        %plot small negative values
        [i,j] = find(matrix < 0 & matrix > maxnegval/2);
        plotMatrix(i,j,pointsize,[1.0 0.5 0.5])
        
        %plot large negative values
        [i,j] = find(matrix < maxnegval/2);
        plotMatrix(i,j,pointsize,[1.0 0.0 0.00])
    end
    if (nContours > 0)
        contour(matrix, nContours)
    end
    
 %   title(sprintf('Half of maximum is %.3g',halfmaxval))
    hold off
    
    [m,n] = size(matrix);
    set(gca,'xlim',[0 n+1],'ylim',[0 m+1],'ydir','reverse', ...
       'grid','none','plotboxaspectratio',[n+1 m+1 1]);
end

function plotMatrix(i,j,size,color)
    if (isempty(i)), return; end
    
    plot(j,i,'s',...
        'MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',size);
end
    
