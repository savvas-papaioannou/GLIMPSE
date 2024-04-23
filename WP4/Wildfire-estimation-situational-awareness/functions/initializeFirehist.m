function firehist = initializeFirehist(fusion,sim,fire)

    firehist.U = struct([]);
    firehist.theta = struct([]);
    for i=1:sim.fire.weather.numberofchange
        % firehist.U{i} = fire.U{i}.*sim.firehist.scenarios(sim.firehist.scenarioID);
        % firehist.theta{i} = fire.theta{i}.*sim.firehist.scenarios(sim.firehist.scenarioID);
        firehist.U{i} = abs(normrnd(fusion.U_mean{i},fusion.U_sigma{i})).*sim.firehist.scenarios(sim.firehist.scenarioID);
        firehist.theta{i} = abs(normrnd(fusion.theta_mean{i},fusion.theta_sigma{i})).*sim.firehist.scenarios(sim.firehist.scenarioID);
    end
    firehist.R = abs(normrnd(fusion.R_mean,fusion.R_sigma)).*sim.firehist.scenarios(sim.firehist.scenarioID);
    

    s = linspace(0,2*pi,sim.fire.N0)';
    a = 300*sim.firehist.scenarios(sim.firehist.scenarioID);
    b = 300*sim.firehist.scenarios(sim.firehist.scenarioID);
    c = 300*sim.firehist.scenarios(sim.firehist.scenarioID);
    firehist.q{1} = [sim.dt*a*cos(s)+sim.fire.InitialGuessX...
        sim.dt*b*sin(s)+sim.dt*c+sim.fire.InitialGuessY];

    %firehist.q{1} = fire.q{1} + normrnd(0,5,[length(fire.q{1}),2]);

    firehist.N = zeros(sim.Steps,1);
    firehist.N(1) = sim.fire.N0;
    firehist.perimeterSize = zeros(sim.Steps,1);
    firehist.areaSize = zeros(sim.Steps,1);

    % firehist.R = fire.R;
    % firehist.U = fire.U;
    % firehist.theta = fire.theta;

    % firehist.R = fusion.R_mean;
    % firehist.U = fusion.U_mean;
    % firehist.theta = fusion.theta_mean;
    % firehist.q{1} = fire.q{1};

end