
%%  Wildfire Estimation System
%   Author: Constantinos Heracleous
%   Last Updated: 07/03/2024
%   KIOS CoE
clc
clear variables
addpath('functions')
clfall

%% Main
sim = simParametersAnisotropic(true,0.4,30); % options: simParameters simParametersIsotropic simParametersAnisotropic
fire = initializeFire(sim);
fusion = initializeFusion(fire,sim);
%firehist = initializeFirehist(fusion,sim,fire);
uav = initializeUAV(sim);

firesimulation.forecast = struct([]);
firesimulation.analysis = struct([]);
firesimulation.future = struct([]);
firesimulation.temp = struct([]); 
fire.yo = struct([]);
mapping.Oburn = struct([]);
planning.path = struct([]);
planning.waypoints = struct([]);

for i=1:sim.N_e
    firesimulation.forecast{i} = initializeFireForecast(fusion,sim);
end

%% Weather change update
for kf=2:(sim.fireSteps+sim.Horizon)
    sim = updateParams(sim, kf);
end

dbox = waitbar(0,'Running... 0%');

counter = 0;
kf = 2;

tic

for k=2:sim.Steps
    %% waitbar update
    waitbar(k/sim.Steps,dbox,sprintf('Running... %0.0f%%',k/sim.Steps*100));

    %% UAV simulation

    % UAV control
    if uav.deploy
        % Check if current waypoint is reached
        if norm(uav.states(1:3,k-1) - uav.current_waypoints(uav.current_waypoint_index, :)') < sim.uav.errTh && uav.current_waypoint_index < size(uav.current_waypoints, 1)
            uav.current_waypoint_index = uav.current_waypoint_index + 1; % Move to next waypoint
            for i=1:3
                uav.mpc1.OV(i).Min = min(uav.states(i,k-1),uav.current_waypoints(uav.current_waypoint_index, i)) - sim.uav.errTh;
                uav.mpc1.OV(i).Max = max(uav.states(i,k-1),uav.current_waypoints(uav.current_waypoint_index, i)) + sim.uav.errTh;
            end
        end
        
        % Set the reference trajectory to the current waypoint
        ref_trajectory = [uav.current_waypoints(uav.current_waypoint_index, :)'; uav.states(4:6,k-1)];
        old_status = mpcverbosity('off');
        uav.inputs(:,k) = mpcmove(uav.mpc1, uav.mpc_state, uav.states(:,k-1), ref_trajectory); % Compute control input
        mpcverbosity(old_status);
        uav.states(:,k) = uav.A * uav.states(:,k-1) + uav.B * uav.inputs(:,k); % Update system state
    else
        uav.states(:,k) = uav.states(:,k-1);
    end
    uav.FOV{k} = calculateFOV(uav.states(:,k),uav.AFOV);
    uav.waypoints{k} = uav.current_waypoints;


    % Fusion
    if uav.deploy
        fusion = makeFusion(fusion,uav.FOV{k},fire,sim);
    end

    %% Fire simulation 
    
    % Increment counter for running fire simulation 
    counter = counter + 1;
    
    if counter >= sim.UAViterationsPerFire
        
        % True Fire
        fire = fireModel(sim,fire,kf);

        fire.yo{kf} = fire.q{kf} + normrnd(0,sim.Sigma,size(fire.q{kf}));
        fire.yo{kf}(end,:) = fire.yo{kf}(1,:);

        % reshape(fire.q{kf}',1,[])
        
        %% Free Run
        %firehist = fireModel(sim,firehist,kf);
    
        %% EnKF
        for i=1:sim.N_e
            firesimulation.forecast{i} = updateFireEstimationParameters(firesimulation.forecast{i},fusion,sim);
            firesimulation.forecast{i} = fireModel(sim,firesimulation.forecast{i},kf);
        end
    
        if kf==sim.RunEnKF
            sim.RunEnKF = sim.RunEnKF + sim.UpdateRate;
            [enkf.largerEnsembleSize,enkf.largerEnsembleIndex] = max(arrayfun(@(x) x{1}.N(kf), firesimulation.forecast));   % Gets the size of ensembles and their max index
            enkf.fireMeasurmentsSize = size(fire.yo{kf},1);
            if enkf.largerEnsembleSize >= enkf.fireMeasurmentsSize  % Check which is larger vector
                for i = 1:sim.N_e
                    if i~=enkf.largerEnsembleIndex
                        firesimulation.forecast{i}.q{kf} = makeinterpolation(firesimulation.forecast{enkf.largerEnsembleIndex}.q{kf},firesimulation.forecast{i}.q{kf});
                        firesimulation.forecast{i}.N(kf) = size(firesimulation.forecast{i}.q{kf},1);
                    end
                end
                fire.yo{kf} = makeinterpolation(firesimulation.forecast{enkf.largerEnsembleIndex}.q{kf},fire.yo{kf});
            else
                for i = 1:sim.N_e
                    firesimulation.forecast{i}.q{kf} = makeinterpolation(fire.yo{kf},firesimulation.forecast{i}.q{kf});
                    firesimulation.forecast{i}.N(kf) = size(firesimulation.forecast{i}.q{kf},1);
                end
            end
            % Run EnKF
            [enkf.x_fbar,enkf.x_abar,enkf.x_ensemble] = ensembleKalmanFilter(sim.N_e,firesimulation.forecast,fire.yo{kf},sim.Sigma,kf);
            firesimulation.q_fbar{kf} = reshape(enkf.x_fbar,2,[])';
            for i=1:sim.N_e
                firesimulation.analysis{i}.q{kf} = clipping_fcn_Matlab(reshape(enkf.x_ensemble(:,i),2,[])');
                firesimulation.forecast{i}.q{kf} = firesimulation.analysis{i}.q{kf};
                firesimulation.forecast{i}.N(kf) = size(firesimulation.analysis{i}.q{kf},1);
                %% Future Fire Propagation
                firesimulation.future{i} = updateFireEstimationParameters(firesimulation.forecast{i},fusion,sim);
                firesimulation.future{i}.q{kf} = firesimulation.analysis{i}.q{kf};
                %firesimulation.future{i}.q{kf} = fire.q{kf};
                if kf<=sim.fireCutoff
                    for kfut=kf+1:kf+sim.Horizon
                        firesimulation.future{i} = fireModel(sim,firesimulation.future{i},kfut);

                        % save the fire in the selected propagation horizon
                        if kfut == kf+sim.propagationHorizon
                            firesimulation.temp{i}.q = firesimulation.future{i}.q{kfut};
                            firesimulation.temp{i}.N = size(firesimulation.future{i}.q{kfut},1);
                        end
                    end
                    firesimulation.futuresave{i}.q{kf} = firesimulation.future{i}.q{kfut};
                    firesimulation.futuresave{i}.N(kf) = size(firesimulation.futuresave{i}.q{kf},1);
                end
            end
            firesimulation.q_abar{kf} = clipping_fcn_Matlab(reshape(enkf.x_abar,2,[])');
            firesimulation.Oanalysis{kf} = getBinaryOccupancyMatrix(firesimulation.q_abar{kf},sim.fire.gridsizeXaxis,...
                sim.fire.gridsizeYaxis,sim.fire.gridResolution);
            if kf<=sim.fireCutoff
                mapping.Oburn{kf} = zeros(sim.fire.matrixsizeX,sim.fire.matrixsizeY);
                [enkf.largerEnsembleSize,enkf.largerEnsembleIndex] = max(arrayfun(@(x) x{1}.N(kf), firesimulation.futuresave));
                firesimulation.q_sbar{kf} = firesimulation.futuresave{enkf.largerEnsembleIndex}.q{kf};

                [enkf.propagationlargerEnsembleSize,enkf.propagationlargerEnsembleIndex] = max(arrayfun(@(x) x{1}.N, firesimulation.temp));
                firesimulation.q_pbar{kf} = firesimulation.temp{enkf.propagationlargerEnsembleIndex}.q;

                for i = 1:sim.N_e
                    if i~=enkf.largerEnsembleIndex
                        firesimulation.futuresave{i}.q{kf} = makeinterpolation(firesimulation.futuresave{enkf.largerEnsembleIndex}.q{kf},firesimulation.futuresave{i}.q{kf});
                        firesimulation.futuresave{i}.N(kf) = size(firesimulation.futuresave{i}.q{kf},1);
                        firesimulation.q_sbar{kf} = firesimulation.q_sbar{kf} + firesimulation.futuresave{i}.q{kf};
                    end
                    Oburn_temp = getBinaryOccupancyMatrix(firesimulation.futuresave{i}.q{kf},sim.fire.gridsizeXaxis,...
                        sim.fire.gridsizeYaxis,sim.fire.gridResolution);
                    Oburn_temp = Oburn_temp - (Oburn_temp & firesimulation.Oanalysis{kf});
                    mapping.Oburn{kf} = mapping.Oburn{kf} + Oburn_temp;

                    
                    if i~=enkf.propagationlargerEnsembleIndex
                        firesimulation.temp{i}.q = makeinterpolation(firesimulation.temp{enkf.propagationlargerEnsembleIndex}.q,firesimulation.temp{i}.q);
                        firesimulation.q_pbar{kf} = firesimulation.q_pbar{kf} + firesimulation.temp{i}.q;
                    end
                end
                firesimulation.q_sbar{kf} = firesimulation.q_sbar{kf}./sim.N_e;
                firesimulation.q_pbar{kf} = firesimulation.q_pbar{kf}./sim.N_e;
                mapping.Oburn{kf} = mapping.Oburn{kf}/sim.N_e;
                % Sigma_R = normalize(fusion.R_sigma,'range');
                % Sigma_U = normalize(fusion.U_sigma{1},'range');
                % Sigma_theta = normalize(fusion.theta_sigma{1},'range');
                Sigma_R = fusion.R_sigma;
                Sigma_U = fusion.U_sigma{1};
                Sigma_theta = fusion.theta_sigma{1};
                
                Sigma_tot = sim.weights(1)*Sigma_R + sim.weights(2)*Sigma_U + sim.weights(3)*Sigma_theta;
                mapping.CoverMap{kf} = normalize(mapping.Oburn{kf} .* Sigma_tot,'range');
                mapping.CoverMapBinary{kf} = mapping.CoverMap{kf} >= sim.CoverThreshold;
                
                planning.waypointsCovered{kf} = uav.current_waypoint_index;
                if uav.current_waypoint_index==0
                    planning.waypointsTocover{kf} = 0;
                else
                    planning.waypointsTocover{kf} = size(planning.waypoints{kf-sim.UpdateRate},1);
                end
            
                if sim.runMethod
                    planning.path{kf} = findPath(firesimulation.Oanalysis{kf},mapping.CoverMapBinary{kf},uav.states(1:2,k),sim.fire.gridResolution);
                    planning.waypoints{kf} = [planning.path{kf}*sim.fire.gridResolution ones(size(planning.path{kf},1),1)*sim.uav.flyAltitude];

                    uav.waypoints{k} = [uav.states(1:3,k)';planning.waypoints{kf}];
                    % deploy UAV
                    if ~uav.deploy
                        uav.deploy = true;
                        sim.start = k;
                    end
                    % prepare uav for the new path
                    uav.current_waypoints = planning.waypoints{kf};
                    uav.current_waypoint_index = 1;
                    for i=1:3
                        uav.mpc1.OV(i).Min = min(uav.states(i,k),uav.current_waypoints(uav.current_waypoint_index, i)) - sim.uav.errTh;
                        uav.mpc1.OV(i).Max = max(uav.states(i,k),uav.current_waypoints(uav.current_waypoint_index, i)) + sim.uav.errTh;
                    end
                end
            end
        end
        counter = 0; % Reset counter for fire simulation
        kf = kf + 1;
    end

    %% Calculate error

    %sim.error_hist(k) = getError(fire.q{k},firehist.q{k});
    % sim.rmse_error_hist(k) = getRMSE(fire.q{k},firehist.q{k});
    % sim.rmse_error_obs(k) = getRMSE(fire.q{k},fire.yo{k});
end

close(dbox);

step=sim.updateStep/sim.fire.dt;
m=1;
for kp=sim.startEnKF:step:kf-sim.Horizon+1
    % Get Metrics
    [firesimulation.metrics.sdiCurrent(m),firesimulation.metrics.adiCurrent(m),...
        firesimulation.metrics.jcCurrent(m),firesimulation.metrics.ssCurrent(m)] = ....
        getMetrics(fire.q{kp},firesimulation.q_abar{kp},sim);
    [firesimulation.metrics.sdiFuture(m),firesimulation.metrics.adiFuture(m),...
        firesimulation.metrics.jcFuture(m),firesimulation.metrics.ssFuture(m)] = ....
        getMetrics(fire.q{kp+sim.Horizon-1},firesimulation.q_sbar{kp},sim);
    firesimulation.metrics.time(m) = (kp-1)*sim.fire.dtmins;

    [firesimulation.metrics.sdiPropagation(m),firesimulation.metrics.adiPropagation(m),...
        firesimulation.metrics.jcPropagation(m),firesimulation.metrics.ssPropagation(m)] = ....
        getMetrics(fire.q{kp+sim.propagationHorizon-1},firesimulation.q_pbar{kp},sim);

    % [firehist.metrics.sdiCurrent(m),firehist.metrics.adiCurrent(m),...
    %     firehist.metrics.jcCurrent(m),firehist.metrics.ssCurrent(m)] = ....
    %     getMetrics(fire.q{kp},firehist.q{kp},sim);
    % [firehist.metrics.sdiFuture(m),firehist.metrics.adiFuture(m),...
    %     firehist.metrics.jcFuture(m),firehist.metrics.ssFuture(m)] = ....
    %     getMetrics(fire.q{kp+sim.Horizon-1},firehist.q{kp+sim.Horizon-1},sim);
    % firehist.metrics.time(m) = (kp-1)*sim.dt;

    m=m+1;
end

uav.FOVCoords.x = cell2mat(cellfun(@(c) c(1,:), uav.FOV(1:k), 'UniformOutput', false)')';
uav.FOVCoords.y = cell2mat(cellfun(@(c) c(2,:), uav.FOV(1:k), 'UniformOutput', false)')';
toc

if sim.runMethod
    makeanimation(sim,fire,uav,sim.Cutoff)
end

offline_plots

res.fire = fire;
res.metrics = firesimulation.metrics;
res.firesimulation.q_fbar = firesimulation.q_fbar;
res.firesimulation.q_abar = firesimulation.q_abar;
res.firesimulation.q_sbar = firesimulation.q_sbar;
res.firesimulation.q_pbar = firesimulation.q_pbar;
res.firesimulation.Oanalysis = firesimulation.Oanalysis;
res.fusion = fusion;
res.uav.states = uav.states;
res.uav.inputs = uav.inputs;
res.uav.waypoints = uav.waypoints;
res.uav.FOV = uav.FOV;
res.uav.FOVCoords = uav.FOVCoords;
res.mapping = mapping;
res.planning = planning;

fileName = ['res' num2str(sim.CoverThreshold) '_' num2str(sim.propagationHorizon_mins) '.mat'];
save(fileName, 'res');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%