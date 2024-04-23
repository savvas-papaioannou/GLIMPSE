function plotData(plt,sim,fire,uav,k)
    global x y
    n=0;
    n=n+1;
    set(plt{n},'XData',fire.q{k}(:,x),'YData',fire.q{k}(:,y));
    n=n+1;
    set(plt{n},'XData',fire.qm{k}(:,x),'YData',fire.qm{k}(:,y));
    n=n+1;
    set(plt{n},'XData',fire.qest{k}(:,x),'YData',fire.qest{k}(:,y));
    n=n+1;
    set(plt{n},'XData',fire.qpre{k}(:,x),'YData',fire.qpre{k}(:,y));
    for i=1:sim.numberofUAVs
        n=n+1;
        set(plt{n},'XData',uav{i}.p(x,k), 'YData', uav{i}.p(y,k));
        n=n+1;
        set(plt{n},'XData',uav{i}.FOV{k}(x,:), 'YData', uav{i}.FOV{k}(y,:));
        n=n+1;
        set(plt{n},'XData',uav{i}.pgoal(x,k), 'YData', uav{i}.pgoal(y,k));
        n=n+1;
        set(plt{n},'XData',uav{i}.target(x,k), 'YData', uav{i}.target(y,k));
    end


%     n=n+1;
%     set(plt{n},'XData',uav.p(x,k), 'YData', uav.p(y,k));
%     n=n+1;
%     set(plt{n},'XData',uav.FOV{k}(x,:),'YData',uav.FOV{k}(y,:));
%     n=n+1;
%     set(plt{n},'XData',uav.pgoal(x,k), 'YData', uav.pgoal(y,k));
%     n=n+1;
%     set(plt{n},'XData',uav.target(x,k), 'YData', uav.target(y,k));
%     if uav.mode(k)==4
%         uavalt = 0;
%     else
%         uavalt = uav.p(z,k);
%     end
    title(['k=' num2str(k) ' Time = '  num2str(k*sim.dt,'%.2f') ' min., Fire perimeter = ' num2str(fire.perimeter(k),'%.0f')...
        'm, Wind Speed = ' num2str(fire.U_mean,'%.0f'),...
        'm/min, Wind Direction = ' num2str(rad2deg(fire.theta_mean),'%.0f') 'deg'])
    %num2str(mod(k,60))'sec
    % wind direction arrow

%     delete(findall(gcf,'type','annotation'))
%     xa = [0.85 0.85+0.035*sin(fire.theta_mean)]; %[x_begin x_end]
%     ya = [0.85 0.85+0.035*cos(fire.theta_mean)]; %[y_begin y_end]
%     annotation('arrow',xa,ya);
%     annotation('textbox',[.85 .2 .02 .02],'String','Wind Direction','LineStyle','none');
    drawnow;
end

