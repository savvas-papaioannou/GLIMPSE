function q = fireStateFcn(q,u) 
% fireStateFcn Discrete-time approximation of simplified FARSITE model 
%
% qt1 = fireStateFcn(qt,ut)
%
% Inputs:
%    qt - States q[t] [x;y] coordinates of firefront j [m]
%    ut - States u[t] [R; U; theta; xs, ys, Dt] 
%       R-> surface fire spread rate [m/min]
%       U-> effective midflame windspeed [m/s],
%       theta-> wind azimuth angle direction of maximum spread rate [radians]
%       xs
%       ys
%       Dt-> time step [min]
%
% Outputs:
%    qt1 - Propagated states q[t+1]
%
% @Constantinos Heracleous, Jan 2022, KIOS CoE

% Euler integration of continuous-time dynamics x'=f(x) with sample time dt
Dt = u(6); % [min] Sample time
q = q + Dt*fireStateFcnContinuous(q,u);
end