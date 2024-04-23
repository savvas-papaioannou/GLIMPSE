function uav = initializeUAV(sim)

    %   UAV initialization function

    xi = sim.uav.dt / sim.uav.m;
    uav.A = [eye(3) sim.uav.dt * eye(3); zeros(3) sim.uav.rho * eye(3)];    % Define A matrix
    uav.B = [zeros(3); xi * eye(3)];                                % Define B matrix
    uav.C = [eye(6)];                                               % Output matrix (assuming full-state feedback)
    uav.D = zeros(6, 3);                                            % D matrix

    plant = ss(uav.A, uav.B, uav.C, uav.D, sim.uav.dt); % State-space model with sampling time dt
    
    old_status = mpcverbosity('off');
    %% create MPC controller object with sample time
    uav.mpc1 = mpc(plant, sim.uav.dt);
    %% specify prediction horizon
    uav.mpc1.PredictionHorizon = 10;
    %% specify control horizon
    uav.mpc1.ControlHorizon = 2;
    %% specify nominal values for inputs and outputs
    uav.mpc1.Model.Nominal.U = [0;0;0];
    uav.mpc1.Model.Nominal.Y = [0;0;0;0;0;0];
    %% specify constraints for MV and MV Rate
    % mpc1.MV(1).RateMin = -40;
    % mpc1.MV(1).RateMax = 40;
    % mpc1.MV(2).RateMin = -40;
    % mpc1.MV(2).RateMax = 40;
    % mpc1.MV(3).RateMin = -40;
    % mpc1.MV(3).RateMax = 40;
    %% specify constraints for MV and MV Rate
    uav.mpc1.MV(1).Min = -150;
    uav.mpc1.MV(1).Max = 150;
    uav.mpc1.MV(2).Min = -150;
    uav.mpc1.MV(2).Max = 150;
    uav.mpc1.MV(3).Min = -150;
    uav.mpc1.MV(3).Max = 150;
    %% specify constraints for OV
    uav.mpc1.OV(4).Min = -23;
    uav.mpc1.OV(4).Max = 23;
    uav.mpc1.OV(5).Min = -23;
    uav.mpc1.OV(5).Max = 23;
    uav.mpc1.OV(6).Min = -6;
    uav.mpc1.OV(6).Max = 6;
    %% specify overall adjustment factor applied to weights
    beta = 4.0552;
    %% specify weights
    uav.mpc1.Weights.MV = [0 0 0]*beta;
    uav.mpc1.Weights.MVRate = [0.1 0.1 0.1]/beta;
    uav.mpc1.Weights.OV = [1 1 1 1 1 1]*beta;
    uav.mpc1.Weights.ECR = 100000;
    %% use custom output disturbance model
    mpc1_ModelOD = [];
    setoutdist(uav.mpc1, 'model', mpc1_ModelOD);
    noise = load('noise_model.mat');
    %% use custom measurement noise model
    uav.mpc1.Model.Noise = noise.mpc1_ModelMN_1;


    uav.mpc_state = mpcstate(uav.mpc1); % Initialize MPC state
    mpcverbosity(old_status);
    
    % UAV camera system spesifications
    uav.Sx = 24;        % horizontal sensor size [mm]
    uav.Sy = 24;        % vertical sensor size [mm]
    uav.f = 24;         % focal length [mm]

    uav.states = zeros(6,sim.Steps);
    uav.states(:,1) = sim.uav.initialState;
    uav.inputs = zeros(3,sim.Steps);
    
    uav.deploy = false;
    uav.current_waypoint_index = 0;
    uav.current_waypoints = [];
    uav.waypoints = struct([]);

    uav.AFOV = [2*atan(uav.Sx/(2*uav.f)) 2*atan(uav.Sy/(2*uav.f))]; % angular field of view [rad]
    uav.FOV{1} = calculateFOV(uav.states(1:3,1),uav.AFOV);
end