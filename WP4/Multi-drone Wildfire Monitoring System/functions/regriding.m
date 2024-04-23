function [X,Y,N] = regriding(x,y,n,T)
% Regriding function that adds extra points to high curvature areas to
% increase the points density and avoid errors.
% inputs:
%           x - vector with x cordinates
%           y - vector with y cordinates
%           n - vectors size
%           T - threshold
% outputs:
%           X - vector with x cordinates after regriding
%           Y - vector with x cordinates after regriding
%           N - vectors size

    %% initializations
    i=1; X = []; Y = [];
    
%     l(k) = norm([x(k) y(k)] - [x(n-1) y(n-1)]);
%     a = [x(k) y(k)] - [x(n-1) y(n-1)];
%     b = [x(k) y(k)] - [x(k+1) y(k+1)];
%     th(k) = acos(dot(a,b)/(norm(a)*norm(b)));
    
    for k=1:n
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
        l(k) = norm(B - A);
        a = B - A;
        b = B - C;
        th(k) = acos(dot(a,b)/(norm(a)*norm(b)));
        if k>1
            X(i) = A(1);
            Y(i) = A(2);
            i=i+1;
            [X,Y,i] = addPoints(B,A,th(k),th(k-1),l(k),T,i,X,Y);
        end
    end
    X(i)=x(k);
    Y(i)=y(k);
    X=X';
    Y=Y';
    N=i;
end

function [X,Y,i] = addPoints(B,A,th1,th2,l,T,i,X,Y)
% Recursive function for adding points in high curved areas.
% inputs:
%           B - vector [x(k) y(k)]
%           A - vector [x(k-1) y(k-1)]
%           th1 - angle of k point in rad
%           th2 - angle of k-1 point in rad
%           l - length at between k - k-1 points.
%           T - threshold
%           i - index
% outputs:
%           X - vector with x cordinates after regriding
%           Y - vector with x cordinates after regriding
%           i - index
    if max([cos(th1/2) cos(th2/2)]) > (T/l)^2
        [X,Y,i] = addPoints((B+A)/2,A,pi,th2,l/2,T,i,X,Y);
        X(i) = (B(1)+A(1))/2;
        Y(i) = (B(2)+A(2))/2;
        i=i+1;
        [X,Y,i] = addPoints(B,(B+A)/2,th1,pi,l/2,T,i,X,Y);
    end
end