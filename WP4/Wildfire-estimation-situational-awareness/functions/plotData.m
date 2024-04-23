function plotData(plt,sim,fire,uav,k,kf)
    x=1; y=2;
    n=0;
    n=n+1;
    set(plt{n},'XData',uav.FOVCoords.x(:,1:k),'YData', uav.FOVCoords.y(:,1:k));
    n=n+1;
    set(plt{n},'XData',fire.q{kf}(:,x),'YData',fire.q{kf}(:,y));
    n=n+1;
    set(plt{n},'XData',uav.waypoints{k}(:,x),'YData',uav.waypoints{k}(:,y));
    n=n+1;
    set(plt{n},'XData',uav.states(x,k),'YData',uav.states(y,k));
    n=n+1;
    set(plt{n},'XData',uav.FOV{k}(x,:),'YData',uav.FOV{k}(y,:));
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
    title(['Time: '  num2str(k*sim.dtmins,'%.2f mins')])% 'k: ' num2str(k) ' |...
        % ' min | W. Speed: ' num2str(sim.fire.U_mean(sim.fire.weather.idx(k))) ' m/s | W. Angle: ' num2str(rad2deg(sim.fire.theta_mean(sim.fire.weather.idx(k)))) ' deg. | Per. = ' num2str(fire.perimeterSize(k)/1000,'%.0f')...
        % 'km | Area: ' num2str(fire.areaSize(k)/ 10000,'%.0f') 'ha'])%, Wind Speed = ' num2str(fire.U_mean,'%.0f'),...
        %'m/min, Wind Direction = ' num2str(rad2deg(fire.theta_mean),'%.0f') 'deg'])
    %num2str(mod(k,60))'sec
    % wind direction arrow

%     delete(findall(gcf,'type','annotation'))
%     xa = [0.85 0.85+0.035*sin(fire.theta_mean)]; %[x_begin x_end]
%     ya = [0.85 0.85+0.035*cos(fire.theta_mean)]; %[y_begin y_end]
%     annotation('arrow',xa,ya);
%     annotation('textbox',[.85 .2 .02 .02],'String','Wind Direction','LineStyle','none');
    drawnow;
end

