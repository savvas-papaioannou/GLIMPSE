function [intX,intY,intSegments,intN] = checkLinesIntersect(x,y)
    n = length(x);
    i=0; intX=[]; intY=[];intSegments=[];
    for k=1:n-2
        for m=k+1:n-1
            [res,intP] = doLineSegmentsIntersect([x(k) y(k)],[x(k+1) y(k+1)], [x(m) y(m)],[x(m+1) y(m+1)]);
            if res
                i=i+1;
                intSegments{i} = [k,k+1,m,m+1];
                intSegDiff(1)= m+1-k;
                intX(i)=intP(1);
                intY(i)=intP(2);
            end
        end
    end
    intN=i;
end

function [res,intP] = doLineSegmentsIntersect(p,p2,q,q2)
    %% Checks if the two segments are intersecting and provides the intersection point
    %  The solution is based on https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect 
    % inputs:
    %           p - point 1 (x,y) of line segment A
    %           p2 - point 2 (x,y) of line segment A
    %           p - point 1 (x,y) of line segment A
    %           p2 - point 2 (x,y) of line segment A
    % outputs:
    %           res - true if lines intersect and false of they don't
    %           intP - intersection point, NaN if there is no intersection
    %           point
    r = p2-p;
    s = q2-q;
    uNumerator = crossProduct((q-p),r);
    Denominator = crossProduct(r,s);

    if Denominator == 0 && uNumerator == 0 || Denominator == 0
		% segments are collinear or parallel
        res = false;
        intP = NaN;
    else
        t = crossProduct((q-p),s)/Denominator;
        u = uNumerator/Denominator;
        if (t > 0) && (t < 1) && (u > 0) && (u < 1)
            res = true;
            intP = p + t*r;
        else
            res = false;
            intP = NaN;
        end
    end
end