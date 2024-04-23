
%%  Multi-drone Wildfire Monitoring System
%   Author: Constantinos Heracleous
%   Last Updated: 31/05/2023
%   KIOS CoE
clc
clear variables
clfall

global x y z;
x=1; y=2; z=3;
addpath('functions')

%% Main
sim.T = 90;                     % simulation time in minutes
sim.t = 0;                      % simulation running time
sim.dt = 1/60;                  % time step in minutes
sim.Steps = sim.T/sim.dt;       % simulation steps
sim.numberofUAVs = 5;           % number of UAVs
%sim.initialUAVposition = [8000 3000 2000 5000 3000; 3000 8000 10000 5000 4000 6000];
%sim.initialUAVposition = [10000 10000 10000 10000 10000 10000; 10000 10000 10000 10000 10000 10000];
% sim.UAVdistanceradius = 1000;   % distance radius of the UAVs from the fire center [m]
% sim.UAVangle = 0;
% sim.initialUAVposition = [ones(1,sim.numberofUAVs)*sim.UAVdistanceradius*cos(sim.UAVangle);...
%     ones(1,sim.numberofUAVs)*sim.UAVdistanceradius*sin(sim.UAVangle)];



fire = initializeFire(sim);
uav = initializeUAV(sim,fire.center.x,fire.center.y);
% plt = initializePlot(fire,uav);
hw = waitbar(0,'Running... 0%');

% figure(1)
% surf(1:fire.gridSize.Xaxis,1:fire.gridSize.Yaxis,fire.R)
tic

for k=1:sim.Steps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % waitbar update
%     if mod(k,floor(sim.Steps/10))<1e-2
        str1 = sprintf('Running... %0.0f%%',k/sim.Steps*100);
        waitbar(k/sim.Steps,hw,str1);
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Wildfire Simulation Part

    % Wind Direction and Magnitude Change Scenario
    if k==20/sim.dt
        fire.theta_mean = deg2rad(285);
        fire.theta = abs(normrnd(fire.theta_mean,fire.theta_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));
        fire.thetahis = abs(normrnd(fire.theta_mean,fire.theta_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));
        fire.thetafus = fire.thetahis;
        %         fire.U_mean = 30;            % mean value of wind speed (m/min)
%         fire.U_sigma = 3;         % std of wind speed (m/min)
%         fire.U = abs(normrnd(fire.U_mean,fire.U_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));
    end

    if k==40/sim.dt
        fire.theta_mean = deg2rad(75);
        fire.theta = abs(normrnd(fire.theta_mean,fire.theta_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));
        fire.thetahis = abs(normrnd(fire.change_percent*fire.theta_mean,fire.theta_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));
        fire.thetafus = fire.thetahis;
%         fire.U_mean = 40;            % mean value of wind speed (m/min)
%         fire.U_sigma = 4;         % std of wind speed (m/min)
%         fire.U = abs(normrnd(fire.U_mean,fire.U_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));
    end

    if k==60/sim.dt
        fire.theta_mean = deg2rad(25);
        fire.theta = abs(normrnd(fire.theta_mean,fire.theta_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));
        fire.thetahis = abs(normrnd(fire.change_percent*fire.theta_mean,fire.theta_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));
        fire.thetafus = fire.thetahis;
%         fire.U_mean = 80;            % mean value of wind speed (m/min)
%         fire.U_sigma = 8;         % std of wind speed (m/min)
%         fire.U = abs(normrnd(fire.U_mean,fire.U_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));
    end

%     if k==84/sim.dt
%         fire.theta_mean = deg2rad(285);
%         fire.theta = normrnd(fire.theta_mean,fire.theta_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis);
%         fire.U_mean = 80;            % mean value of wind speed (m/min)
%         fire.U_sigma = 8;         % std of wind speed (m/min)
%         fire.U = abs(normrnd(fire.U_mean,fire.U_sigma,fire.gridSize.Xaxis,fire.gridSize.Yaxis));
%     end

    %% Fire Evolution Calculation
    fire = fireEvolution(sim,fire,k);

    %% Clipping
    [X,Y,N] = clipping_fcn_v5(fire.q{k+1}(:,x),fire.q{k+1}(:,y),...
        fire.N,fire.theta_mean);
    if N~=-1
        fire.N = N;
        fire.q{k+1} = [];
        for i=1:fire.N
            fire.q{k+1}(i,x) = X(i);
            fire.q{k+1}(i,y) = Y(i);
        end
    end

    %% Regriding
    [X,Y,N] = regriding(fire.q{k+1}(:,x),fire.q{k+1}(:,y),fire.N,fire.T);
    if N~=fire.N
        fire.N = N;
        fire.q{k+1} = [];
        for i=1:fire.N
            fire.q{k+1}(i,x) = X(i);
            fire.q{k+1}(i,y) = Y(i);
        end
    end

    %% Clipping & Regriding estimated fire
    if fire.est
        [X,Y,N] = clipping_fcn_v5(fire.qest{k+1}(:,x),fire.qest{k+1}(:,y),...
        fire.Nest,fire.theta_mean);
        if N~=-1
            fire.Nest = N;
            fire.qest{k+1} = [];
            for i=1:fire.Nest
                fire.qest{k+1}(i,x) = X(i);
                fire.qest{k+1}(i,y) = Y(i);
            end
        end
    
        [X,Y,N] = regriding(fire.qest{k+1}(:,x),fire.qest{k+1}(:,y),fire.Nest,fire.T);
        if N~=fire.Nest
            fire.Nest = N;
            fire.qest{k+1} = [];
            for i=1:fire.Nest
                fire.qest{k+1}(i,x) = X(i);
                fire.qest{k+1}(i,y) = Y(i);
            end
        end
    end

    %% Clipping & Regriding predicted fire
    [X,Y,N] = clipping_fcn_v5(fire.qpre{k+1}(:,x),fire.qpre{k+1}(:,y),...
    fire.Npre,fire.theta_mean);
    if N~=-1
        fire.Npre = N;
        fire.qpre{k+1} = [];
        for i=1:fire.Npre
            fire.qpre{k+1}(i,x) = X(i);
            fire.qpre{k+1}(i,y) = Y(i);
        end
    end

    [X,Y,N] = regriding(fire.qpre{k+1}(:,x),fire.qpre{k+1}(:,y),fire.Npre,fire.T);
    if N~=fire.Npre
        fire.Npre = N;
        fire.qpre{k+1} = [];
        for i=1:fire.Npre
            fire.qpre{k+1}(i,x) = X(i);
            fire.qpre{k+1}(i,y) = Y(i);
        end
    end

    fire.perimeter(k) = sum(sqrt((fire.q{k}(2:end,x)-fire.q{k}(1:end-1,x)).^2 + ...
        (fire.q{k}(2:end,y)-fire.q{k}(1:end-1,y)).^2));

    % Fire vectors initialization
    if (k>1)
        fire.qm{k} = fire.qm{k-1};
        if (~fire.est)
            fire.qest{k} = fire.qest{k-1};
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% UAVs Monitoring and Control Part
    for i=1:sim.numberofUAVs
        %% FOV calculation
        uav{i}.FOV{k} = calculateFOV(uav{i}.p(:,k),uav{i}.AFOV);
        xfov = round(uav{i}.FOV{k}(x,:));
        yfov = round(uav{i}.FOV{k}(y,:));

        switch uav{i}.phase(k)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Depot phase
            case 1
                % UAV control
                if uav{i}.mode(k) == 1  % ready
                    if round(sim.t(k),2)>=uav{i}.deployTime
                        uav{i}.phase(k) = 2;
                        uav{i}.mode(k) = 2;
                        fprintf([num2str(sim.t(k),'%.1f') ' min. - UAV' num2str(i) ' deployed.\n'])
                        if (~fire.est)
                            fire.est = true;
                            fire.qm{k} = fire.q{k};
                            fire.qest{k+1} = fire.q{k+1};
                            fire.Nest = fire.N;
                        end
                    end
                else                    % charging
                % 
                end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Infield phase
            case 2
                % UAV control
                if uav{i}.mode(k) == 2  % deploy
                    uav{i}.pgoal(:,k) = [fire.center.x;fire.center.y;uav{i}.deployAltitude];
                    q_detected = inpolygon(fire.q{k}(:,x),fire.q{k}(:,y),uav{i}.FOV{k}(x,:),uav{i}.FOV{k}(y,:));
                    if sum(q_detected) > 0
                        uav{i}.phase(k) = 3;
                        uav{i}.mode(k) = 3;
                        fprintf([num2str(sim.t(k),'%.1f') ' min. - UAV' num2str(i) ' has reached fire.\n']);
                    end
                else                    % retired
                    
                end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Fire Tracking phase
            case 3
                % UAV control
                if uav{i}.mode(k) == 3  % path_planing
                    fprintf([num2str(k) ' ' num2str(sim.t(k),'%.1f') ' min. - UAV' num2str(i) ' initialize path planing... ']);
                    % find possible waypoints
                    [Xpway,Ypway,Zpway,Xsgmt,Ysgmt,Xsgmt_plot,Ysgmt_plot,dist_way] =...
                        findPossibleWaypoints(fire.qm{k}(:,x),fire.qm{k}(:,y),fire.qest{k}(:,x),fire.qest{k}(:,y),uav{i}.AFOV(x),1);
                    A = generateAdjecencyMatrix(Xpway,Ypway);
                    Npway=length(Xpway); % pway - possible waypoints
                    lpway = vecnorm([Xpway(1:Npway);Ypway(1:Npway)]-[Xpway(2:Npway) Xpway(1);Ypway(2:Npway) Ypway(1)]); % length of possible waypoints
                    
                    dpway = [];
                    lpuav = [];
                    index = [];
                    
                    numberofFTUAVs=0; %number of firetracking UAVs
                    for l=1:sim.numberofUAVs
                        if uav{l}.phase(k) == 3 %|| uav{l}.phase(k) == 2
                            numberofFTUAVs=numberofFTUAVs+1;
                            FTUAVs(numberofFTUAVs) = l; %firetracking UAVs
                        end
                    end

                    for l=1:numberofFTUAVs
                        [lpuav(FTUAVs(l),:),index(FTUAVs(l),:)] = sort(sqrt(sum(([uav{FTUAVs(l)}.p(x,k);uav{FTUAVs(l)}.p(y,k)]-[Xpway;Ypway]).^2)));
                    end
    

                    % update UAVs posible waypoints start indexes to avoid
                    % two or more UAVs to have the same starting point and 
                    % fail to partition the perimeter accordingly
                    loop=0;
                    while loop<100
                        loop=loop+1;
                        flag=true;
                        for l=1:numberofFTUAVs-1
                            for j=l+1:numberofFTUAVs
                                if index(FTUAVs(l),1) == index(FTUAVs(j),1)
                                    flag=false;
                                    if lpuav(FTUAVs(l),1)<=lpuav(FTUAVs(j),1)
                                        index(FTUAVs(j),:) = circshift(index(FTUAVs(j),:),-1);
                                        lpuav(FTUAVs(j),:) = circshift(lpuav(FTUAVs(j),:),-1);
                                    else
                                        index(FTUAVs(l),:) = circshift(index(FTUAVs(l),:),-1);
                                        lpuav(FTUAVs(l),:) = circshift(lpuav(FTUAVs(l),:),-1);
                                    end
                                end
                            end
                        end
                        if flag
                            break;
                        end
                    end
                    start_node = index(:,1)';
    
                    for l=1:numberofFTUAVs
                        dpway_ccw=zeros(1,Npway);
                        dpway_cw=zeros(1,Npway);
                        %[~,start_node(l)] = min(sqrt(sum(([uav{l}.p(x,k);uav{l}.p(y,k)]-[Xpway;Ypway]).^2)));
                        lpway_ccw = [lpway(start_node(FTUAVs(l)):Npway) lpway(1:start_node(FTUAVs(l))-1)];
                        lpway_cw = flip(lpway_ccw);
                        for j=2:Npway
                            dpway_ccw(j) = dpway_ccw(j-1)+lpway_ccw(j-1);
                            dpway_cw(j) = dpway_cw(j-1)+lpway_cw(j-1);
                        end
                        dpway_ccw = circshift(dpway_ccw,start_node(FTUAVs(l))-1);
                        dpway_cw = circshift(flip(dpway_cw),start_node(FTUAVs(l)));
                        dpway(FTUAVs(l),:) = min(dpway_cw,dpway_ccw);
                    end
    
                    total_dpway = zeros(1,Npway);
                    uav_partition = zeros(1,Npway);
                    for j=1:Npway
                        temp_dpway = Inf;
                        for l=1:numberofFTUAVs
                            if dpway(FTUAVs(l),j)<temp_dpway
                                temp_dpway = dpway(FTUAVs(l),j);
                                temp_uav = FTUAVs(l);
                            end
                        end
                        total_dpway(j) = temp_dpway;
                        uav_partition(j) = temp_uav;
                    end
    
                    msg = [];
                    for l=1:numberofFTUAVs
                        index = find(uav_partition==FTUAVs(l));
                        [~,index_max] = max(dist_way(uav_partition==FTUAVs(l)));
                        target_node(FTUAVs(l)) = index(index_max);%+kappa_s(i)-1;
    
                        [uav{FTUAVs(l)}.waypoints,uav{FTUAVs(l)}.waypointsSize,uav{FTUAVs(l)}.path] = getDesirePath(A,start_node(FTUAVs(l)),target_node(FTUAVs(l)),Xpway,Ypway,Zpway);
                        for j=1:uav{FTUAVs(l)}.waypointsSize
                            uav{FTUAVs(l)}.waypoints(3,j) = Zpway(target_node(FTUAVs(l)));
                        end
                        uav{FTUAVs(l)}.next_waypoint=1;
                        uav{FTUAVs(l)}.target(:,k) = [Xpway(target_node(FTUAVs(l)));Ypway(target_node(FTUAVs(l)))];
                        uav{FTUAVs(l)}.pgoal(:,k) = uav{FTUAVs(l)}.waypoints(:,uav{FTUAVs(l)}.next_waypoint);
                        uav{FTUAVs(l)}.mode(k)=4;
                        msg = [msg 'UAV' num2str(FTUAVs(l)) ' '];
                    end
                    fprintf(['completed,\n           ' msg 'now following the new path.\n']);
                end

% figure(4)
% clf
% hold on
% color={'r','k','g','m','b','y'};
% plot(fire.q{k}(:,x),fire.q{k}(:,y),fire.qm{k}(:,x),fire.qm{k}(:,y),fire.qest{k}(:,x),fire.qest{k}(:,y));
% plot(Xsgmt_plot,Ysgmt_plot);
% plot(Xpway,Ypway,'o');
% for l=1:sim.numberofUAVs
%     plot(uav{l}.p(x,k),uav{l}.p(y,k),'x','MarkerEdgeColor',color{l});
%     plot(uav{l}.FOV{k-1}(x,:),uav{l}.FOV{k-1}(y,:),'-','Color',color{l});
%     if uav{l}.phase(k) == 3
%         plot(Xpway(uav_partition==l),Ypway(uav_partition==l),'*','MarkerEdgeColor',color{l});
%         %plot([uav{l}.p(x,k) uav{l}.waypoints(x,:)],[uav{l}.p(y,k) uav{l}.waypoints(y,:)],'-','Color',color{l});
%     end
% end
% grid on

                if uav{i}.mode(k) == 4  % path_following
                    if (norm(uav{i}.p(:,k)-uav{i}.waypoints(:,uav{i}.next_waypoint))<=uav{i}.targetError)
                        uav{i}.next_waypoint = uav{i}.next_waypoint + 1;
                        if uav{i}.next_waypoint > uav{i}.waypointsSize
                            uav{i}.mode(k)=3;
                            uav{i}.next_waypoint = uav{i}.next_waypoint - 1;
                            fprintf([num2str(k) ' ' num2str(sim.t(k),'%.1f') ' min. - UAV' num2str(i) ' reached target waypoint in the current path!\n']);
                        end
                    end
                    uav{i}.pgoal(:,k) = uav{i}.waypoints(:,uav{i}.next_waypoint);
                end    
               
                %% fuse fire measurement and fire estimation perimeters
                % if kmonitor~=k
                temp_qm = fire.qm{k};
                temp_qest = fire.qest{k+1};
                %     kmonitor = k;
                %     fire.qm{k} = fire.qm{k-1};
                % end
                
                % fire detection throught FOV
                q_detected = inpolygon(fire.q{k}(:,x),fire.q{k}(:,y),uav{i}.FOV{k}(x,:),uav{i}.FOV{k}(y,:));
                q1_detected = inpolygon(fire.q{k+1}(:,x),fire.q{k+1}(:,y),uav{i}.FOV{k}(x,:),uav{i}.FOV{k}(y,:));
                qm_detected = inpolygon(temp_qm(:,x),temp_qm(:,y),uav{i}.FOV{k}(x,:),uav{i}.FOV{k}(y,:));
                qest_detected = inpolygon(temp_qest(:,x),temp_qest(:,y),uav{i}.FOV{k}(x,:),uav{i}.FOV{k}(y,:));
                
                % update measurements vector and estimations vector
                if sum(q_detected) > 0 && sum(qm_detected) > 0 && sum(qest_detected) > 0 && sum(q1_detected) > 0
                    dq = find(q_detected);
                    dqm = find(qm_detected);
                    [difdq,indexq_max]=max(dq(2:end)-dq(1:end-1));
                    [difdqm,indexqm_max]=max(dqm(2:end)-dqm(1:end-1));
                    
                    % update measurements vector
                    if difdqm>size(qm_detected,1)/4
                        if difdq>size(q_detected,1)/4
                            temp_qmX = [fire.q{k}(dq(indexq_max),x); temp_qm(dqm(indexqm_max)+1:dqm(indexqm_max+1)-1,x);...
                                fire.q{k}(dq(indexq_max+1):end,x); fire.q{k}(dq(1):dq(indexq_max),x)];

                            temp_qmY = [fire.q{k}(dq(indexq_max),y); temp_qm(dqm(indexqm_max)+1:dqm(indexqm_max+1)-1,y);...
                                fire.q{k}(dq(indexq_max+1):end,y); fire.q{k}(dq(1):dq(indexq_max),y)];
                        else
                            temp_qmX = [fire.q{k}(dq,x); temp_qm(dqm(indexqm_max)+1:dqm(indexqm_max+1)-1,x); fire.q{k}(dq(1),x)];
                            temp_qmY = [fire.q{k}(dq,y); temp_qm(dqm(indexqm_max)+1:dqm(indexqm_max+1)-1,y); fire.q{k}(dq(1),y)];
                        end
                    else
                        old_start = find(qm_detected,1);
                        old_finish = find(qm_detected,1,'last');
                        if difdq>size(q_detected,1)/4
                            temp_qmX = [temp_qm(1:old_start-1,x); fire.q{k}(dq(indexq_max+1:end),x);...
                                fire.q{k}(dq(2:indexq_max),x); temp_qm(old_finish+1:end,x)];
                            temp_qmY = [temp_qm(1:old_start-1,y); fire.q{k}(dq(indexq_max+1:end),y);...
                                fire.q{k}(dq(2:indexq_max),y); temp_qm(old_finish+1:end,y)];
                        else
                            new_start = find(q_detected,1);
                            new_finish = find(q_detected,1,'last');
                            temp_qmX = [temp_qm(1:old_start-1,x); fire.q{k}(q_detected,x); temp_qm(old_finish+1:end,x)];
                            temp_qmY = [temp_qm(1:old_start-1,y); fire.q{k}(q_detected,y); temp_qm(old_finish+1:end,y)];
                        end
                    end
                    temp_qm = [temp_qmX temp_qmY];
                    fire.qm{k} = temp_qm;

                    % update estimations vector
                    dq1 = find(q1_detected);
                    dqest = find(qest_detected);
                    [difdq1,indexq1_max]=max(dq1(2:end)-dq1(1:end-1));
                    [difdqest,indexqest_max]=max(dqest(2:end)-dqest(1:end-1));

                    qest = temp_qest;

                    if difdqest>size(qest_detected,1)/4
                        if difdq1>size(q1_detected,1)/4
                            temp_qestX = [fire.q{k+1}(dq1(indexq1_max),x); qest(dqest(indexqest_max)+1:dqest(indexqest_max+1)-1,x);...
                                fire.q{k+1}(dq1(indexq1_max+1):end,x); fire.q{k+1}(dq1(1):dq1(indexq1_max),x)];
                            temp_qestY = [fire.q{k+1}(dq1(indexq1_max),y); qest(dqest(indexqest_max)+1:dqest(indexqest_max+1)-1,y);...
                                fire.q{k+1}(dq1(indexq1_max+1):end,y); fire.q{k+1}(dq1(1):dq1(indexq1_max),y)];
                        else
                            temp_qestX = [fire.q{k+1}(dq1,x); qest(dqest(indexqest_max)+1:dqest(indexqest_max+1)-1,x); fire.q{k+1}(dq1(1),x)];
                            temp_qestY = [fire.q{k+1}(dq1,y); qest(dqest(indexqest_max)+1:dqest(indexqest_max+1)-1,y); fire.q{k+1}(dq1(1),y)];
                        end
                    else
                        old_start = find(qest_detected,1);
                        old_finish = find(qest_detected,1,'last');
                        if difdq1>size(q1_detected,1)/4
                            temp_qestX = [qest(1:old_start-1,x); fire.q{k+1}(dq1(indexq1_max+1:end),x);...
                                fire.q{k+1}(dq1(2:indexq1_max),x); qest(old_finish+1:end,x)];
                            temp_qestY = [qest(1:old_start-1,y); fire.q{k+1}(dq1(indexq1_max+1:end),y);...
                                fire.q{k+1}(dq1(2:indexq1_max),y); qest(old_finish+1:end,y)];
                        else
                            new_start = find(q1_detected,1);
                            new_finish = find(q1_detected,1,'last');
                            
                            temp_qestX = [qest(1:old_start-1,x); fire.q{k+1}(q1_detected,x); qest(old_finish+1:end,x)];
                            temp_qestY = [qest(1:old_start-1,y); fire.q{k+1}(q1_detected,y); qest(old_finish+1:end,y)];
                        end
                    end
                    temp_qest = [temp_qestX temp_qestY];
                    fire.qest{k+1} = temp_qest;
                    fire.Nest = size(fire.qest{k+1},1);
                    uav{i}.emptyFOV = 0;
                else
                    uav{i}.emptyFOV = uav{i}.emptyFOV + 1;
                    if uav{i}.emptyFOV >= uav{i}.emptyFOVThreshold
                        fprintf([num2str(k) ' ' num2str(sim.t(k),'%.1f') ' min. - UAV' num2str(i) ' empty FOV for a while!\n']);
                        uav{i}.mode(k)=3;
                        uav{i}.emptyFOV = 0;
                    end
                end

        end

        
        %% Measurment process
        fire.Rfus(xfov(1):xfov(2),yfov(4):yfov(1)) = fire.R(xfov(1):xfov(2),yfov(4):yfov(1));
        fire.Ufus(xfov(1):xfov(2),yfov(4):yfov(1)) = fire.U(xfov(1):xfov(2),yfov(4):yfov(1));
        fire.thetafus(xfov(1):xfov(2),yfov(4):yfov(1)) = fire.theta(xfov(1):xfov(2),yfov(4):yfov(1));

        %% Control process
        uav{i}.phase(k+1) = uav{i}.phase(k);
        uav{i}.mode(k+1) = uav{i}.mode(k);
        uav{i}.target(:,k+1) = uav{i}.target(:,k);
        [uav{i}.v(:,k),uav{i}.w(:,k)] = control(uav{i}.p(:,k),uav{i}.pgoal(:,k),uav{i}.vref);
        uav{i}.p(:,k+1) = uav{i}.p(:,k) + uav{i}.v(:,k).*sim.dt;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %plotData(plt,sim,fire,uav,k)
    sim.t(k+1) = sim.t(k) + sim.dt;
    %hold on
    %pause(0.1)
end
close(hw);
toc
kp = k;

makeanimation(sim,fire,uav,k)

figure(2)
hold on
lbl={};
i=0;
for t=0:10:sim.T
    k = t/sim.dt+1;
    plot(fire.q{k}(:,x),fire.q{k}(:,y))
%     plot(fire.q{k}(1,x),fire.q{k}(1,y),'x')
%     plot(fire.q{k}(2,x),fire.q{k}(2,y),'*')
%     plot(fire.q{k}(3,x),fire.q{k}(3,y),'*')
    i=i+1;
    lbl{i} = [num2str(t) ' min'];
end
legend(lbl)
grid on

% figure(3)
% [row,col] = find(fire.Rfus);
% hold on
% plot(row,col,'o')
% plot(uav.p(x,:),uav.p(y,:),'o-')
% grid on


figure(4)
hold on
plot(fire.q{kp}(:,x),fire.q{kp}(:,y))
plot(fire.qm{kp}(:,x),fire.qm{kp}(:,y))
plot(fire.qest{kp}(:,x),fire.qest{kp}(:,y))
plot(fire.qpre{kp}(:,x),fire.qpre{kp}(:,y))
grid on
legend('Real Fire','Fire Measurments','Fire Estimation',...
    'Fire Prediction (historic data)',Location='southeast');

% figure(3)
% hold on
% plot(fire.q{kf}(:,x),fire.q{kf}(:,y),fire.qm{kf}(:,x),fire.qm{kf}(:,y));
% plot(Xsgmt_plot,Ysgmt_plot);
% plot(Xpway,Ypway,'o');
% plot(waypoints(x,:),waypoints(y,:),'-*');
% plot(waypoints(x,1),waypoints(y,1),'g*');
% plot(waypoints(x,end),waypoints(y,end),'r*');
% plot(uav.p(x,kf),uav.p(y,kf),'rx');
% grid on
%}

%% Functions in /functions folder

function clfall
    FigList = findall(groot, 'Type', 'figure');
    for iFig = 1:numel(FigList)
        try
            clf(FigList(iFig));
        catch
            % Nothing to do
        end
    end
end
