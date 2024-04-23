function plt = initializePlot(sim,fire,uav)
    global x y
    figure(1);
    color={'r','k','g','m','b','y'};
    hold on
    n=0;
    n=n+1;
    plt{n} = plot(fire.q{1}(:,x),fire.q{1}(:,y));
    n=n+1;
    plt{n} = plot(fire.qm{1}(:,x),fire.qm{1}(:,y));
    n=n+1;
    plt{n} = plot(fire.qest{1}(:,x),fire.qest{1}(:,y));
    n=n+1;
    plt{n} = plot(fire.qpre{1}(:,x),fire.qpre{1}(:,y));
    for i=1:sim.numberofUAVs
        n=n+1;
        plt{n} = plot(uav{i}.p(x,1),uav{i}.p(y,1),'x','MarkerEdgeColor',color{i});
        n=n+1;
        plt{n} = plot(uav{i}.FOV{1}(x,:),uav{i}.FOV{1}(y,:),color{i});
        n=n+1;
        plt{n} = plot(uav{i}.pgoal(x,1),uav{i}.pgoal(y,1),'*','MarkerEdgeColor',color{i});
        n=n+1;
        plt{n} = plot(uav{i}.target(x,1),uav{i}.target(y,1),'+','MarkerEdgeColor',color{i});
    end

    axis([0 fire.gridSize.Xaxis 0 fire.gridSize.Yaxis]);
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