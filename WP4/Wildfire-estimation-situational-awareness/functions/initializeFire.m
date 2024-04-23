function fire = initializeFire(sim)
    rng(sim.seed,'twister');
    
    fire.U = struct([]);
    fire.theta = struct([]);
    for i=1:sim.fire.weather.numberofchange
        fire.U{i} = getWindMagnitude(sim.fire.U_mean(i),sim.fire.U_sigma(i),sim);
        fire.theta{i} = getWindAngle(sim.fire.theta_mean(i),sim.fire.theta_sigma(i),sim);
    end
    fire.R = zeros(sim.fire.matrixsizeX,sim.fire.matrixsizeY);

    % Rate of spread calculation via normal distribution
    % fire.R = abs(normrnd(sim.fire.R_mean,sim.fire.R_sigma,...
    %   sim.fire.gridsizeXaxis,sim.fire.gridSizeYaxis));
    
    % Rate of spread calculation via constant cells
    fire.R_mean = randi(sim.fire.R_mean_max_value,sim.fire.RcellNum);
    fire.R_sigma = randi(sim.fire.R_sigma_max_value,sim.fire.RcellNum);
    %fire.R_sigma = sim.Rcoef*fire.R_mean;

    Nx = round(linspace(1,sim.fire.matrixsizeX,sim.fire.RcellNum+1));
    Ny = round(linspace(1,sim.fire.matrixsizeY,sim.fire.RcellNum+1));
    
    for i=1:length(Nx)-1
        for j=1:length(Ny)-1
            fire.R(Nx(i):Nx(i+1),Ny(j):Ny(j+1)) = ...
                abs(normrnd(fire.R_mean(i,j),fire.R_sigma(i,j),Nx(i+1)-Nx(i)+1,Ny(j+1)-Ny(j)+1));
        end
    end

    % Visualize Rate of spread
    % [X,Y] = meshgrid(1:sim.fire.matrixsizeX,1:sim.fire.matrixsizeY);
    % Z = fire.R;
    % s = surf(X,Y,Z);
    % % Enhancements for better visualization
    % shading interp; % Interpolate colors across faces and edges
    % colorbar;   % Display a color bar to indicate the scale of values
    
    % initial fire fronts curve
    s = linspace(0,2*pi,sim.fire.N0)';
    a = 300/sim.minTosec;
    b = 300/sim.minTosec;
    c = 300/sim.minTosec;
    fire.q{1} = [sim.dt*a*cos(s)+sim.fire.InitialX...
        sim.dt*b*sin(s)+sim.dt*c+sim.fire.InitialY];
    fire.N = zeros(sim.Steps,1);
    fire.N(1) = sim.fire.N0;
    fire.perimeterSize = zeros(sim.Steps,1);
    fire.areaSize = zeros(sim.Steps,1);
end