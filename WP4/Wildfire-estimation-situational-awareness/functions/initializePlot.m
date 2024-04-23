function plt = initializePlot(sim,fire,uav)
    x=1; y=2;
    figure(1);
    %color={'r','k','g','m','b','y'};
    hold on
    n=0;
    n=n+1;
    plt{n} = fill(uav.FOVCoords.x(:,1), uav.FOVCoords.y(:,1),'yellow','FaceAlpha',0.5, 'EdgeColor', 'none');
    n=n+1;
    plt{n} = plot(fire.q{1}(:,x),fire.q{1}(:,y),LineWidth=1);
    n=n+1;
    plt{n} = scatter(0,0,5,"red","filled");
    n=n+1;
    plt{n} = plot(uav.states(1,x),uav.states(1,y),'kx',MarkerSize=10);
    n=n+1;
    plt{n} = plot(uav.FOV{1}(x,:),uav.FOV{1}(y,:),'-k');
    axis([0 sim.fire.gridsizeXaxis 0 sim.fire.gridsizeYaxis]);
%     axis([850 1150 150 400]);
    ylabel('y - [m]')
    xlabel('x - [m]')
%     title(['Time = '  num2str(0) 'min:' num2str(0)...
%         'sec,   Fire Perimeter = ' num2str(0)...
%         'm,   UAV altitude = ' num2str(0) 'm'])
    %legend('q','q_m','p_{target}','p','fov')
    grid on
%    legend('Real Fire','Fire Measurments','Fire Estimation','Fire Prediction (historic data)',...
%      'UAV 2D position','UAV FoV','Next Waypoint','Target',Location='southeast');