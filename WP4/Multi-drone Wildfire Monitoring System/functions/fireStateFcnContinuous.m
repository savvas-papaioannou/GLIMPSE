function Q = fireStateFcnContinuous(q,u)
%fireStateFcnContinuous Evaluate differentials
R = u(1);
U = u(2);
theta = u(3);
xs = u(4);
ys = u(5);

LB = 0.936*exp(0.2566*U) + 0.461*exp(-0.1548*U) - 0.397;
HB = (LB + (LB^2 - 1)^0.5) / (LB - (LB^2 - 1)^0.5);

a = (0.5*(R+(R/HB)))/LB;
b = (R+(R/HB))/2;
c = b - R/HB;

% a = 0.02;
% b = 0.04;
% c = 0.03;

lowPart = sqrt(b^2*(xs*cos(theta)-ys*sin(theta))^2 + a^2*(xs*sin(theta)+ys*cos(theta))^2);
xupPart = a^2*cos(theta)*(xs*sin(theta)+ys*cos(theta)) - b^2*sin(theta)*(xs*cos(theta)-ys*sin(theta));
yupPart = -a^2*sin(theta)*(xs*sin(theta)+ys*cos(theta)) - b^2*cos(theta)*(xs*cos(theta)-ys*sin(theta));

if lowPart==0 || xupPart==0
    xRatio = 0;
else
    xRatio = xupPart/lowPart;
end

if lowPart==0 || yupPart==0
    yRatio = 0;
else
    yRatio = yupPart/lowPart;
end

Xt = xRatio + c*sin(theta);
Yt = yRatio + c*cos(theta);
Q = [Xt Yt];
end