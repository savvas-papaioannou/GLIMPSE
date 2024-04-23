function [X,Y,N] = clipping_fcn(x,y,n)
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
    
    [intX,intY,intSegments] = checkLinesIntersect(x,y);
    if numel(intX) == 0
        X=[];Y=[];N=-1;
        return
    end

    z=length(intX);
    X=[];Y=[];
    for i=1:z
        if i==1
            X = x(1:intSegments{i}(1));
            Y = y(1:intSegments{i}(1));
        else
            X = [X;x(intSegments{i-1}(4):intSegments{i}(1))];
            Y = [Y;y(intSegments{i-1}(4):intSegments{i}(1))];
        end
        X = [X;intX(i)];
        Y = [Y;intY(i)];
    end
    X = [X;x(intSegments{i}(4):n)];
    Y = [Y;y(intSegments{i}(4):n)];

    %         if intSegments{i}(4) ~= n
    %             X = [X;x(intSegments{i}(4):length(x))];
    %             Y = [Y;y(intSegments{i}(4):length(y))];
    %         end
    N = length(X);
end