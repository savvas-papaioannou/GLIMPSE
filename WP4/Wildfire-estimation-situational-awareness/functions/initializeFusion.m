function fusion = initializeFusion(fire,sim)
    % Nx = round(linspace(1,sim.fire.matrixsizeX,sim.fire.RcellNum+1));
    % Ny = round(linspace(1,sim.fire.matrixsizeY,sim.fire.RcellNum+1));
    % 
    % for i=1:length(Nx)-1
    %     for j=1:length(Ny)-1
    %         fusion.R_mean(Nx(i):Nx(i+1),Ny(j):Ny(j+1)) = fire.R_mean(i,j);
    %         fusion.R_sigma(Nx(i):Nx(i+1),Ny(j):Ny(j+1)) = fire.R_sigma(i,j);
    %     end
    % end
    %



    expansion_matrix = ones(sim.fire.matrixsizeX,sim.fire.matrixsizeY);
    for i=1:sim.fire.weather.numberofchange
        fusion.U_mean{i} = sim.fire.U_mean(i).*expansion_matrix;
        fusion.U_sigma{i} = sim.fire.U_sigma(i).*expansion_matrix;

        fusion.theta_mean{i} = sim.fire.theta_mean(i).*expansion_matrix;
        fusion.theta_sigma{i} = sim.fire.theta_sigma(i).*expansion_matrix;
    end

    fire.R = zeros(sim.fire.matrixsizeX,sim.fire.matrixsizeY);

    % Rate of spread calculation via normal distribution
    % fire.R = abs(normrnd(sim.fire.R_mean,sim.fire.R_sigma,...
    %   sim.fire.gridsizeXaxis,sim.fire.gridSizeYaxis));

    % Rate of spread calculation via constant cells
    R_mean = randi(sim.fire.R_mean_max_value,sim.fire.RcellNum);
    R_sigma = randi(sim.fire.R_sigma_max_value,sim.fire.RcellNum);

    Nx = round(linspace(1,sim.fire.matrixsizeX,sim.fire.RcellNum+1));
    Ny = round(linspace(1,sim.fire.matrixsizeY,sim.fire.RcellNum+1));

    for i=1:length(Nx)-1
        for j=1:length(Ny)-1
            fusion.R_mean(Nx(i):Nx(i+1),Ny(j):Ny(j+1)) = R_mean(i,j);
            fusion.R_sigma(Nx(i):Nx(i+1),Ny(j):Ny(j+1)) = R_sigma(i,j);
        end
    end



    % load("fusion0.01.mat");

    % fusion.R_mean = fire.R;
    % fusion.R_sigma = zeros(sim.fire.matrixsizeX,sim.fire.matrixsizeY);
    % fusion.U_mean{1} = fire.U{1};
    % fusion.U_sigma{1} = zeros(sim.fire.matrixsizeX,sim.fire.matrixsizeY);
    % fusion.theta_mean{1} = fire.theta{1};
    % fusion.theta_sigma{1} = zeros(sim.fire.matrixsizeX,sim.fire.matrixsizeY);
end