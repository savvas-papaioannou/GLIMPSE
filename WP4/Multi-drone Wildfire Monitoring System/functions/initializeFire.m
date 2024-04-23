function fire = initializeFire(sim)
    fire.gridSize.Xaxis = 22000;% X dim in [m]
    fire.gridSize.Yaxis = 22000;% Y dim in [m]
    fire.N0 = 80;               % number of initial fire fronts points
    fire.N = fire.N0;           % total number of firefront points
    fire.T = 35;                % [m] threshold for regriding
    
    
    fire.U_mean = 2;          % mean value of wind speed (m/s)
    fire.U_sigma = 0.1*fire.U_mean;         % std of wind speed (m/s)
    fire.U = abs(normrnd(fire.U_mean,fire.U_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));
    fire.theta_mean = deg2rad(45);       % mean value of wind angle (0-360 degrees)
    fire.theta_sigma = deg2rad(5);       % std of wind angle
    fire.theta = abs(normrnd(fire.theta_mean,fire.theta_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));

    fire.R = zeros(fire.gridSize.Xaxis,fire.gridSize.Yaxis);

    % Rate of spread calculation via normal distribution in all area
%     fire.R_mean = 180;         % mean value of fire spreading (m/min) 180m/min = 3m/s = 10.8km/h
%     fire.R_sigma = 80;         % std of fire spreading (m/min)
%     fire.R = abs(normrnd(fire.R_mean,fire.R_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));

    % Rate of spread calculation via constant cells
    fire.R_mean = [60 100 55 70 200 250 150 180 80 150 180 100 70 90 110 130];
    fire.R_sigma = 0.1*fire.R_mean;

    Nx = [1 floor(fire.gridSize.Xaxis/3) floor(fire.gridSize.Xaxis/2) floor(2*fire.gridSize.Xaxis/3) fire.gridSize.Xaxis];
    Ny = [1 floor(fire.gridSize.Yaxis/2) floor(fire.gridSize.Yaxis/2) floor(2*fire.gridSize.Yaxis/3) fire.gridSize.Yaxis];
    
    m=0;
    for i=1:length(Nx)-1
        for j=1:length(Ny)-1
            m=m+1;
            fire.R(Nx(i):Nx(i+1),Ny(j):Ny(j+1)) = abs(normrnd(fire.R_mean(m),fire.R_sigma(m),Nx(i+1)-Nx(i)+1,Ny(j+1)-Ny(j)+1));
        end
    end

    % Rate of spread calculation via random boxes    
%     N = 50;
%     R_mean_max = 350; % m/min
%     R_sigma_max = 0.1*R_mean_max;
% 
%     for i=1:N
%         fire.R_mean = R_mean_max*rand;
%         fire.R_sigma = R_sigma_max*rand;
%         x_min = randi(fire.gridSize.Xaxis);
%         x_max= randi([x_min fire.gridSize.Xaxis]);
%         y_min = randi(fire.gridSize.Yaxis);
%         y_max= randi([x_min fire.gridSize.Yaxis]);
%         fire.R(x_min:x_max,y_min:y_max) = abs(normrnd(fire.R_mean,fire.R_sigma,x_max-x_min+1,y_max-y_min+1));
%     end

    
    % initial fire fronts curve
    s = linspace(0,2*pi,fire.N0)';
    a = 5000;
    b = 7000;
    c = 5000;
    fire.center.x = 5000;
    fire.center.y = 7000;
    fire.q{1} = [sim.dt*a*cos(s)+fire.center.x sim.dt*b*sin(s)+sim.dt*c+fire.center.y];
    fire.perimeter = [];
    fire.area = [];
    fire.qm{1} = [nan nan];

    % Estimation Settings
    fire.qest{1} = [nan nan];
    fire.est=false;
    fire.Nest = 0;

    fire.change_percent = 0.8;
    % Rate of spread calculation via constant cells
    fire.Rhis_mean = fire.change_percent*fire.R_mean;
    fire.Rhis_sigma = 0.1*fire.Rhis_mean;

    Nx = [1 floor(fire.gridSize.Xaxis/3) floor(fire.gridSize.Xaxis/2) floor(2*fire.gridSize.Xaxis/3) fire.gridSize.Xaxis];
    Ny = [1 floor(fire.gridSize.Yaxis/2) floor(fire.gridSize.Yaxis/2) floor(2*fire.gridSize.Yaxis/3) fire.gridSize.Yaxis];
    
    m=0;
    for i=1:length(Nx)-1
        for j=1:length(Ny)-1
            m=m+1;
            fire.Rhis(Nx(i):Nx(i+1),Ny(j):Ny(j+1)) = abs(normrnd(fire.Rhis_mean(m),fire.Rhis_sigma(m),Nx(i+1)-Nx(i)+1,Ny(j+1)-Ny(j)+1));
        end
    end
    fire.Rfus = fire.Rhis;

    fire.Uhis_mean = fire.change_percent*fire.U_mean;            % mean value of historic wind speed (m/min)
    fire.Uhis_sigma = 0.5*fire.Uhis_mean;         % std of historic wind speed (m/min)
    fire.Uhis = abs(normrnd(fire.U_mean,fire.U_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));

    fire.thetahis_mean = fire.change_percent*fire.theta_mean;       % mean value of historic wind angle (0-360 degrees)
    fire.thetahis_sigma = deg2rad(5);       % std of historic wind angle
    fire.thetahis = abs(normrnd(fire.thetahis_mean,fire.thetahis_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));

    fire.Ufus = fire.Uhis;
    fire.thetafus = fire.thetahis;

    % 1:measuring firefronts
    % 2:measuring firefront & rate of spread
    % 3:measuring firefront & weather
    % 4:measuring firefront & rate of spread & weather
    
    fire.mode = 1;

    % Prediction Settings
    fire.qpre{1} = fire.q{1};
    fire.Npre = fire.N;

end