function sim = simParametersIsotropic()
    %% Main simulation parameters
    
    % Time
    sim.T = 2*60;                   % simulation time in minutes
    sim.dt = 1;
    sim.fire.dt = 1;                     % time step in minutes (1/60 is 1 sec)
    sim.Steps = sim.T/sim.dt;       % simulation steps
    sim.minsTosec = 1/60;
    
    % Scenario
    sim.seed = 1;                   % seed for random variables
    
    % Fire
    sim.fire.gridsizeXaxis = 1000;  % fire environment X grid dim [m]
    sim.fire.gridsizeYaxis = 1000;  % fire environment Y grid dim [m]
    sim.fire.gridResolution = 10;   % fire environment in [m]

    sim.fire.matrixsizeX = round(sim.fire.gridsizeXaxis/sim.fire.gridResolution); % input matrices size X dim
    sim.fire.matrixsizeY = round(sim.fire.gridsizeYaxis/sim.fire.gridResolution); % input matrices size Y dim

    sim.fire.InitialX = 500;       % fire initial X position [m]
    sim.fire.InitialY = 500;       % fire initial Y position [m]

    sim.fire.InitialGuessX = 470;       % fire initial X position [m]
    sim.fire.InitialGuessY = 180;       % fire initial Y position [m]
    
    sim.fire.N0 = 30;               % initial number of fire fronts points
    sim.fire.T = 50;                % [m] threshold for regriding
    
    sim.fire.weather.numberofchange = 1;    % number of weather change should be >=1
    sim.fire.weather.Change = floor(linspace(1,sim.Steps,sim.fire.weather.numberofchange+1)); % weather change windows
    sim.fire.weather.idx = zeros(1,sim.Steps);
    sim.fire.weather.idx(1) = 1;

    sim.fire.U_mean = [0 4 6 5];                    % mean value of wind speed (m/s)
    sim.fire.U_sigma = 0;%0.1*sim.fire.U_mean;         % std of wind speed (m/s)
    
    sim.fire.theta_mean = deg2rad([0 0 45 90]);    % mean value of wind angle (0-360 degrees)
    sim.fire.theta_sigma = 0;%deg2rad([5 5 5 5]);      % std of wind angle
    
    
    % sim.fire.R_mean = 250;            % mean value of fire spreading (m/min) 250m/min = 2m/s = 7.2km/h
    % sim.fire.R_sigma = 0.00001;       % std of fire spreading (m/min)
    
    %%%% Change the constant value 35 set in initializeFire
    sim.fire.R_mean_max_value = 50;     % mean value of fire spreading (m/min) 60m/min = 1m/s = 3.6km/h
    sim.fire.RcellNum = 1;             % constant cell number for rate of spread it should be >=1
    sim.Rcoef = 0.0;                    % std of fire spreading (m/min)

    sim.firehist.scenarioID = 2;
    sim.firehist.scenarios = [0.6 1 1.4];
    
    sim.error_hist = zeros(sim.Steps,1);
    sim.rmse_error_hist = zeros(sim.Steps,1);

    % EnKF parameters
    sim.N_e = 5;

    sim.Sigma = 5;                  % std of measurment noise
    sim.rmse_error_obs = zeros(sim.Steps,1);

    sim.initialPositionSigma = 15;

end