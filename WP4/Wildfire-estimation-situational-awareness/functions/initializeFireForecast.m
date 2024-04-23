function fireforecast = initializeFireForecast(fusion,sim)

    fireforecast.U = struct([]);
    fireforecast.theta = struct([]);
    for i=1:sim.fire.weather.numberofchange
        fireforecast.U{i} = abs(normrnd(fusion.U_mean{i},fusion.U_sigma{i}));
        fireforecast.theta{i} = abs(normrnd(fusion.theta_mean{i},fusion.theta_sigma{i}));
    end
    fireforecast.R = abs(normrnd(fusion.R_mean,fusion.R_sigma));

    s = linspace(0,2*pi,sim.fire.N0)';
    a = 300;
    b = 300;
    c = 300;
    
    X = normrnd(sim.fire.InitialGuessX,sim.initialPositionSigma);
    Y = normrnd(sim.fire.InitialGuessY,sim.initialPositionSigma);
    fireforecast.q{1} = [sim.dt*a*cos(s)+X...
            sim.dt*b*sin(s)+sim.dt*c+Y];

    fireforecast.N = zeros(sim.Steps,1);
    fireforecast.N(1) = sim.fire.N0;
    fireforecast.perimeterSize = zeros(sim.Steps,1);
    fireforecast.areaSize = zeros(sim.Steps,1);
end