function res = myatan(X,Y)
    res = atan2(Y,X);
    if res<0
        res = res+2*pi;
    end
end