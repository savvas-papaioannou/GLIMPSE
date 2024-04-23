 function [X,Y,N] = clipping_fcn_v5(x,y,n,theta) %,flg,flgIndex
    % Function for clipping sections that intersect.
    % inputs:
    %           x - vector with x cordinates
    %           y - vector with y cordinates
    %           intX - 
    %           intY -
    %           intSegments -
    % outputs:
    %           X - 
    %           Y - 
    %           N - 
    
%     concaveResult = getConcavePoints(x,y);
%     concaveResult.angle = (concaveResult.angle*360)/(2*pi);


    %% Reindex the cordinates starting point based on wind direction
    if theta>=0 && theta<=pi
        [~,index_new] = min(x);
    else
        [~,index_new] = max(x);
    end
    
    xrdx = [x(index_new:n-1);x(1:index_new)];
    yrdx = [y(index_new:n-1);y(1:index_new)];
    
    %% Check which line segments intersect
    [intX,intY,intSegments,intN] = checkLinesIntersect(xrdx,yrdx);
    X=[];Y=[];
    if intN == 0
        N=-1;
        return
    end

    %% Flag internal loop intersection segments
    flg = ones(1,intN);
    for i=1:intN-1
        for j=i+1:intN
            if flg(j) && segmentIncluded(intSegments{j}(1),intSegments{j}(4),intSegments{i}(1),intSegments{i}(4))
                flg(j) = 0;
            end
        end
    end
    updintN = sum(flg);
    flgIndex = find(flg);
    updintX = intX(flgIndex);
    updintY = intY(flgIndex);
    for i=1:updintN
        updintSegments{i} = intSegments{flgIndex(i)};
    end

%     intX = intX(flag);
%     intY = intY(flag);
%     intSegments = intSegments

    %% Remove loops
    for i=1:updintN
        if i==1
            X = xrdx(1:updintSegments{i}(1));
            Y = yrdx(1:updintSegments{i}(1));
        else
            X = [X;xrdx(updintSegments{i-1}(4):updintSegments{i}(1))];
            Y = [Y;yrdx(updintSegments{i-1}(4):updintSegments{i}(1))];
        end
        X = [X;updintX(i)];
        Y = [Y;updintY(i)];
    end
    X = [X;xrdx(updintSegments{i}(4):n)];
    Y = [Y;yrdx(updintSegments{i}(4):n)];
    N = length(X);
    
end


function res = segmentIncluded(B1,B2,A1,A2)
    if (B1>=A1) && (B2<=A2)
        res = true;
    else
        res = false;
    end
end


