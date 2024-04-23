function [FOV] = calculateFOV(p,AFOV)
    x=1; y=2; z=3;
    FOVxdim = 2*p(z)*tan(AFOV(x)/2);
    FOVydim = 2*p(z)*tan(AFOV(y)/2);
    
    FOVx = [p(x)-FOVxdim/2 p(x)+FOVxdim/2 p(x)+FOVxdim/2 p(x)-FOVxdim/2 p(x)-FOVxdim/2];
    FOVy = [p(y)+FOVydim/2 p(y)+FOVydim/2 p(y)-FOVydim/2 p(y)-FOVydim/2 p(y)+FOVydim/2];
    FOVz = [0 0 0 0 0];

    FOV = [FOVx;FOVy;FOVz];
    
    % in case of gimbal
%     FOVx0 = [-FOVxdim/2 FOVxdim/2 FOVxdim/2 -FOVxdim/2 -FOVxdim/2];
%     FOVy0 = [FOVydim/2 FOVydim/2 -FOVydim/2 -FOVydim/2 +FOVydim/2];
%     for i=1:length(FOVx0)
%         FOVx(i) = FOVx0(i)*cos(theta)-FOVy0(i)*sin(theta)+p(1);
%         FOVy(i) = FOVx0(i)*sin(theta)+FOVy0(i)*cos(theta)+p(2);
%     end
end