function res = checkIntersection2(X,Y,Xemin,Xemax,Yemin,Yemax,Xamin,Xamax,Yamin,Yamax)
    e=0.000001;
    res1 = X>Xemin && X<Xemax && Y>Yemin && Y<Yemax;
    res2 = X>Xamin+e && X<Xamax-e && Y>Yamin+e && Y<Yamax-e;
    res = res1 && res2;
end
