function res = getConcavePoints(x,y)
    % Function for adding points in high curved areas.
    % inputs:
    %           x - vector with x cordinates
    %           y - vector with y cordinates
    % outputs:
    %           res.index - index of concave points
    %           res.angle - angle of concave points

    n=length(x);
    i=1; res.index = []; res.angle = [];

    for k=1:n
        % Set points A->B->C angle is in B
        % point k=1 is the same with k=n
        if k==1
            A = [x(n-1) y(n-1)];
            C = [x(k+1) y(k+1)];
        elseif k==n
            A = [x(k-1) y(k-1)];
            C = [x(2) y(2)];
        else
            A = [x(k-1) y(k-1)];
            C = [x(k+1) y(k+1)];
        end
        B = [x(k) y(k)];
        % calculate vectors
        a = B - A;  
        b = B - C;
        % find the internal counter clockwise angle
        th =  atan2(a(2), a(1)) - atan2(b(2), b(1));
        if th < 0
            th = 2*pi + th;
        end
        if (th > pi)
            res.index(i) = k;
            res.angle(i) = th;
            i=i+1;
        end
    end
%     n = n-1;
%     total = (n-2) * 180;
%     total_sum = sum(th2d(1:n));
end