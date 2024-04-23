function [v,w] = control(p,p_goal,vref)
    global x y z;
    ua = p_goal-p;
    phi = myatan(ua(x),ua(y));
    theta = myatan(ua(z),sqrt(ua(x)^2+ua(y)^2));
    vx = vref*cos(phi)*sin(theta);
    vy = vref*sin(phi)*sin(theta);
    vz = vref*cos(theta);
    v = [vx;vy;vz];
    w = [theta;phi];
end

