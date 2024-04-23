function sim = simParametersAnisotropic(runMethod,CoverThreshold,propagationHorizon_mins)
    %% Main simulation parameters
    
    % Time
    sim.Tmins = 110;                    %110 total simulation time in minutes
    sim.T =  sim.Tmins*60;              % total simulation time in secs
    sim.dt = 1;                         % time step for simulation in secs
    sim.uav.dt = sim.dt;                % time step for uav in secs
    sim.fire.dt = 10;                   % time step for fire in secs
    sim.minTosec = 60;
    sim.secTomin = 1/60;
    sim.dtmins = sim.dt*sim.secTomin;
    sim.fire.dtmins = sim.fire.dt*sim.secTomin;
    sim.UAViterationsPerFire = sim.fire.dt/sim.uav.dt;

    sim.uav.t = 0:sim.uav.dt:sim.T;
    sim.fire.t = 0:sim.fire.dt:sim.T;

    sim.Steps = sim.T/sim.dt+1;         % fire simulation steps
    sim.fireSteps = sim.T/sim.fire.dt+1;% uav fire simulation steps

    % Scenario
    sim.seed = 1;                   % seed for random variables

    % Method
    tstart_mins = 10;                    % [mins] Initial start of future perimeter estimation
    tstart = tstart_mins*60;
    sim.startEnKF = tstart/sim.fire.dt+1;
    sim.RunEnKF = sim.startEnKF;
    sim.updateStep_mins = 5;             % [mins] Update Step.
    sim.updateStep = sim.updateStep_mins*60;
    sim.UpdateRate = sim.updateStep/sim.fire.dt;      % Time for update.

    thorizon_mins = 30;                  % [mins] Method look-ahead horizon
    thorizon = thorizon_mins*60;
    sim.Horizon = thorizon/sim.fire.dt+1;


    sim.propagationHorizon_mins = propagationHorizon_mins;
    tpropagationHorizon_mins = sim.propagationHorizon_mins;
    tpropagationHorizon = tpropagationHorizon_mins*60;
    sim.propagationHorizon = tpropagationHorizon/sim.fire.dt+1;
    sim.weights = [1 0.5 0.5];

    %sim.Cycles = floor((sim.T-(sim.RunEnKF+sim.Horizon))/sim.UpdateRate);     % Method run cycles
    sim.fireCutoff = (sim.T-thorizon)/sim.fire.dt+1;                          % [mins] Method look-ahead horizon cutoff
    sim.Cutoff = (sim.T-thorizon)/sim.dt;
    
    sim.runMethod = runMethod;       % value to run the method
    sim.CoverThreshold = CoverThreshold;   % selected threshold to cover an area

    % UAV
    sim.uav.flyAltitude = 100;
    sim.uav.N = 1;
    sim.uav.initialState = [1000;2000;0;0;0;0]; % [x;y;z;vx;vy;vz] initial UAV position
    sim.uav.m = 6.5;
    sim.uav.rho = 0.2;
    sim.uav.errTh = 30;


    % Fire
    sim.fire.gridsizeXaxis = 5000;  % fire environment X grid dim [m]
    sim.fire.gridsizeYaxis = 5000;  % fire environment Y grid dim [m]
    sim.fire.gridResolution = 50;   % fire environment in [m]
    sim.fire.gridResolutionMetrics = 5; % fire environment in [m] for metrics

    sim.fire.matrixsizeX = ceil(sim.fire.gridsizeXaxis/sim.fire.gridResolution); % input matrices size X dim
    sim.fire.matrixsizeY = ceil(sim.fire.gridsizeYaxis/sim.fire.gridResolution); % input matrices size Y dim

    sim.fire.InitialX = 1500;       % fire initial X position [m]
    sim.fire.InitialY = 650;       % fire initial Y position [m]


    sim.fire.InitialGuessX = 1570;       % fire initial X position [m]
    sim.fire.InitialGuessY = 670;       % fire initial Y position [m]
    
    sim.fire.N0 = 30;               % initial number of fire fronts points
    sim.fire.T = 70;                % [m] threshold for regriding
    
    sim.fire.weather.numberofchange = 1;    % number of weather change should be >=1
    sim.fire.weather.Change = floor(linspace(1,sim.fireSteps+sim.Horizon,sim.fire.weather.numberofchange+1)); % weather change windows
    sim.fire.weather.idx = zeros(1,sim.fireSteps+sim.Horizon);
    sim.fire.weather.idx(1) = 1;

    sim.fire.U_mean = [1 4 6 5];                    % mean value of wind speed (m/s)
    sim.fire.U_sigma = [10 10 10 10];         % std of wind speed (m/s)
    
    sim.fire.theta_mean = deg2rad([20 0 45 90]);    % mean value of wind angle (0-360 degrees)
    sim.fire.theta_sigma = deg2rad([20 5 5 5]);      % std of wind angle
    
    % sim.fire.R_mean = 250;            % mean value of fire spreading (m/min) 250m/min = 2m/s = 7.2km/h
    % sim.fire.R_sigma = 0.00001;       % std of fire spreading (m/min)
    
    %% Change the constant value 35 set in initializeFire
    sim.fire.R_mean_max_value = 70;    % mean value of fire spreading (m/min) 60m/min = 1m/s = 3.6km/h
    sim.fire.R_sigma_max_value = 80;   % std of fire spreading (m/min)
    sim.fire.RcellNum = 20;            % constant cell number for rate of spread it should be >=1

    %%
    sim.firehist.scenarioID = 2;
    sim.firehist.scenarios = [0.6 1 1.2];

    %% EnKF parameters
    sim.N_e = 150;
    sim.Sigma = 10;                  % std of measurment noise
    sim.rmse_error_obs = zeros(sim.Steps,1);

    sim.initialPositionSigma = 100; %100;
end