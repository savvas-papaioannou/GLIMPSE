function uav = initializeUAV(sim,fireX,fireY)
    %   UAV initialization function

    p0_radius = 1000; % distance radius of the UAVs from the fire center [m]
    p0_angle = pi/4;
%     for i=1:sim.numberofUAVs
%         p0(:,i) = normrnd(mean_p0(i),sigma_p0,2,1); %initial position for UAVs
%     end
    p0 = [ones(1,sim.numberofUAVs)*p0_radius*sin(p0_angle)+fireX;...
        ones(1,sim.numberofUAVs)*p0_radius*cos(p0_angle)+fireY];
%    [~,I] = sort(vecnorm(p0 - [fireX;fireY])); %sort by distance to the fire
    
    for i=1:sim.numberofUAVs
        uav{i}.p0 = [p0(:,i);0]; %initial position of i-th UAV [x;y;z]
    
        uav{i}.deployAltitude = 50; %deploy altitude
        uav{i}.deployTime = i+1;     % deploy time in min

        uav{i}.pgoal = ones(3,sim.Steps).*[uav{i}.p0];
        uav{i}.vref = 900;          %uav reference speed m/min    
        uav{i}.p(:,1) = uav{i}.p0;
        uav{i}.targetError = uav{i}.vref*sim.dt;
        uav{i}.target = [nan;nan];
        uav{i}.emptyFOV = 0;
        uav{i}.emptyFOVThreshold = 30;


        % uav phase: 1-depot, 2-infield, 3-fire_tracking
        uav{i}.phase = 1;
        % uav modes: 1-ready, 2-deploy, 3-path_planing, 4-path_following,
        % 5-retired, 6-charging
        uav{i}.mode = 1;

        % uav control
        uav{i}.decommission = false;
        uav{i}.ready = false;

        uav{i}.v(:,1) = [0;0;0];   % deploy speed vector
        uav{i}.w(:,1) = [0;0];     % angles theta and phi

        % UAV camera system spesifications
        uav{i}.Sx = 24;        % horizontal sensor size [mm]
        uav{i}.Sy = 24;        % vertical sensor size [mm]
        uav{i}.f = 24;         % focal length [mm]
        uav{i}.s = 2;          % scale factor

        uav{i}.AFOV = [2*atan(uav{i}.Sx/(2*uav{i}.f)) 2*atan(uav{i}.Sy/(2*uav{i}.f))]; % angular field of view [rad]
        uav{i}.FOV{1} = calculateFOV(uav{i}.p(:,1),uav{i}.AFOV);
    end
end