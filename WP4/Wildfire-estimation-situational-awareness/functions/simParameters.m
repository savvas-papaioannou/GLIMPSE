function sim = simParameters()
    %% Main simulation parameters
    
    % Time
    sim.T = 60;                     % simulation time in minutes
    sim.dt = 10/60;                  % time step in minutes (1/60 is 1 sec)
    sim.Steps = sim.T/sim.dt;       % simulation steps
    
    % Scenario
    sim.seed = 1;                   % seed for random variables
    
    % Fire
    sim.fire.gridsizeXaxis = 2000;  % fire environment X grid dim [m]
    sim.fire.gridsizeYaxis = 2000;  % fire environment Y grid dim [m]
    sim.fire.gridResolution = 20;   % fire environment in [m]

    sim.fire.matrixsizeX = round(sim.fire.gridsizeXaxis/sim.fire.gridResolution); % input matrices size X dim
    sim.fire.matrixsizeY = round(sim.fire.gridsizeYaxis/sim.fire.gridResolution); % input matrices size Y dim

    sim.fire.InitialX = 600;       % fire initial X position [m]
    sim.fire.InitialY = 600;       % fire initial Y position [m]
    
    sim.fire.N0 = 15;               % initial number of fire fronts points
    sim.fire.T = 50;                % [m] threshold for regriding
    
    sim.fire.weather.numberofchange = 4;
    sim.fire.weather.Change = floor(linspace(1,sim.Steps,sim.fire.weather.numberofchange+1)); % weather change windows
    sim.fire.weather.idx = zeros(sim.Steps,1);
    sim.fire.weather.idx = 1;

    sim.fire.U_mean = [2 4 6 5];                    % mean value of wind speed (m/s)
    sim.fire.U_sigma = 0.1*sim.fire.U_mean;         % std of wind speed (m/s)
    
    sim.fire.theta_mean = deg2rad([45 0 45 90]);    % mean value of wind angle (0-360 degrees)
    sim.fire.theta_sigma = deg2rad([5 5 5 5]);      % std of wind angle
    
    
    % sim.fire.R_mean = 250;            % mean value of fire spreading (m/min) 250m/min = 2m/s = 7.2km/h
    % sim.fire.R_sigma = 0.00001;       % std of fire spreading (m/min)
    
    sim.fire.R_mean_max_value = 45;     % mean value of fire spreading (m/min) 60m/min = 1m/s = 3.6km/h
    sim.fire.RcellNum = 20;             % constant cell number for rate of spread it should be >=1
    sim.Rcoef = 0.1;

    sim.firehist.scenarioID = 3;
    sim.firehist.scenarios = [0.6 1.05 1.4];
    
    sim.error_hist = zeros(sim.Steps,1);
    sim.rmse_error_hist = zeros(sim.Steps,1);


end