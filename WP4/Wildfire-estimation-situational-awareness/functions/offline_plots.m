n=1;
x=1;y=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Propagation Plots
n=n+1;
figure(n)
tl = tiledlayout(5, 3);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

step=sim.updateStep/sim.fire.dt;

for kp=sim.startEnKF:step:kf-sim.Horizon+1
    nexttile
    hold on
    p1=plot(fire.q{kp}(:,x),fire.q{kp}(:,y),'-b',LineWidth=1);
    %p2=plot(firehist.q{kp}(:,x),firehist.q{kp}(:,y),'-k',LineWidth=2);
    %p3=plot(fire.yo{kp}(:,x),fire.yo{kp}(:,y),'+k',LineWidth=2);
    p2=plot(firesimulation.q_abar{kp}(:,x),firesimulation.q_abar{kp}(:,y),'-g',LineWidth=1);

    % for i=1:sim.N_e
    %     p4=plot(firesimulation.futuresave{i}.q{kp}(:,x),firesimulation.futuresave{i}.q{kp}(:,y),'--r');
    % end
    p3=plot(fire.q{kp+sim.Horizon-1}(:,x),fire.q{kp+sim.Horizon-1}(:,y),'--b',LineWidth=1);
    %p5=plot(firehist.q{kp+sim.Horizon-1}(:,x),firehist.q{kp+sim.Horizon-1}(:,y),'--k',LineWidth=2);

    p4=plot(firesimulation.q_sbar{kp}(:,x),firesimulation.q_sbar{kp}(:,y),'--g',LineWidth=1);
    legend([p1 p2 p3 p4], 'True Fire', 'EnKF Estimation',...
        ['True Fire @' num2str((kp+sim.Horizon-2)*sim.fire.dtmins) 'mins'],...
        ['Future Estimated Fire Propagation @' num2str((kp+sim.Horizon-2)*sim.fire.dtmins) 'mins']);
    hold off
    xlim([0 sim.fire.gridsizeXaxis]);
    ylim([0 sim.fire.gridsizeYaxis]);
    title(['Time ' num2str((kp-1)*sim.fire.dtmins) 'mins'])
    grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metrics Plots
n=n+1;
figure(n)
tl = tiledlayout(3, 1);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

nexttile
hold on
p1 = plot(firesimulation.metrics.time,firesimulation.metrics.sdiCurrent,'-o');
p2 = plot(firesimulation.metrics.time,firesimulation.metrics.adiCurrent,'-o');
p3 = plot(firesimulation.metrics.time,firesimulation.metrics.jcCurrent,'-o');
p4 = plot(firesimulation.metrics.time,firesimulation.metrics.ssCurrent,'-o');
legend([p1 p2 p3 p4], 'Shape Deviation Index (SDI)', 'Area Difference Index (ADI)', 'Jaccard Coefficient (JC)', 'Sorensen Similarity (SS)');
title('current vs enkf estimated')
ylim([0 1])
hold off
grid on

nexttile
hold on
p1 = plot(firesimulation.metrics.time,firesimulation.metrics.sdiFuture,'-o');
p2 = plot(firesimulation.metrics.time,firesimulation.metrics.adiFuture,'-o');
p3 = plot(firesimulation.metrics.time,firesimulation.metrics.jcFuture,'-o');
p4 = plot(firesimulation.metrics.time,firesimulation.metrics.ssFuture,'-o');
title('future vs enkf future estimated')
hold off
ylim([0 1])
legend([p1 p2 p3 p4], 'Shape Deviation Index (SDI)', 'Area Difference Index (ADI)', 'Jaccard Coefficient (JC)', 'Sorensen Similarity (SS)');
grid on

nexttile
hold on
p1 = plot(firesimulation.metrics.time,firesimulation.metrics.sdiPropagation,'-o');
p2 = plot(firesimulation.metrics.time,firesimulation.metrics.adiPropagation,'-o');
p3 = plot(firesimulation.metrics.time,firesimulation.metrics.jcPropagation,'-o');
p4 = plot(firesimulation.metrics.time,firesimulation.metrics.ssPropagation,'-o');
title('future vs propagation')
hold off
ylim([0 1])
legend([p1 p2 p3 p4], 'Shape Deviation Index (SDI)', 'Area Difference Index (ADI)', 'Jaccard Coefficient (JC)', 'Sorensen Similarity (SS)');
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Burn Propability Plots
n=n+1;
figure(n)
tl = tiledlayout(5, 3);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

step=sim.updateStep/sim.fire.dt;
for kp=sim.startEnKF:step:kf-sim.Horizon+1
    nexttile;
    % Create a mesh grid
    [x2, y2] = meshgrid(1:sim.fire.matrixsizeY, 1:sim.fire.matrixsizeX);
    % Create the surface plot
    surf(x2, y2, mapping.Oburn{kp});
    view(0, 90);
    % Enhancements for better visualization
    shading interp; % Interpolate colors across faces and edges
    c = colorbar;   % Display a color bar to indicate the scale of values
    c.Label.String = 'Burn Probability (0-1)';
    title(['Burn Probability @ \tau=' num2str((kp-1)*sim.fire.dtmins) 'mins']);
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Burn Probability');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input plots

n=n+1;
figure(n)
tl = tiledlayout(3, 4);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
nexttile;
% Create the surface plot
surf(x2, y2, fusion.R_mean);
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
c.Label.String = 'm/min';
title('Rate of spread mean');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Rate of spread m/min');

nexttile;
% Create the surface plot
surf(x2, y2, fusion.R_sigma);
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
c.Label.String = 'm/min';
title('Rate of spread uncertainty');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Rate of spread m/min');

nexttile;
% Create the surface plot
surf(x2, y2, fire.R);
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
c.Label.String = 'm/min';
title('True Rate of spread');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Rate of spread m/min');

nexttile;
% Create the surface plot
surf(x2, y2, fire.R-fusion.R_mean);
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
c.Label.String = 'm/min';
title('Rate of spread dif. from mean');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Rate of spread m/min');

nexttile;
% Create the surface plot
surf(x2, y2, fusion.U_mean{1});
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
c.Label.String = 'm/s';
title('Wind speed mean');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Wind speed m/s');

nexttile;
% Create the surface plot
surf(x2, y2, fusion.U_sigma{1});
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
c.Label.String = 'm/s';
title('Wind speed uncertainty');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Wind speed m/s');

nexttile;
% Create the surface plot
surf(x2, y2, fire.U{1});
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
c.Label.String = 'm/s';
title('True Wind speed');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Wind speed m/s');


nexttile;
% Create the surface plot
surf(x2, y2, fire.U{1}-fusion.U_mean{1});
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
c.Label.String = 'm/s';
title('Wind speed dif. from mean');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Wind speed m/s');

nexttile;
% Create the surface plot
surf(x2, y2, fusion.theta_mean{1});
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
c.Label.String = 'rad';
title('Wind angle mean');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Wind angle rad');

nexttile;
% Create the surface plot
surf(x2, y2, fusion.theta_sigma{1});
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
c.Label.String = 'rad';
title('Wind angle uncertainty');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Wind angle rad');


nexttile;
% Create the surface plot
surf(x2, y2, fire.theta{1});
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
c.Label.String = 'rad';
title('True Wind angle');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Wind angle rad');

nexttile;
% Create the surface plot
surf(x2, y2, fire.theta{1}-fusion.theta_mean{1});
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
c.Label.String = 'rad';
title('Wind angle dif. from mean');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Wind angle rad');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cover Map Plots
n=n+1;
figure(n)
tl = tiledlayout(5, 3);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

step=sim.updateStep/sim.fire.dt;
for kp=sim.startEnKF:step:kf-sim.Horizon+1
    nexttile;
    % Create a mesh grid
    [x2, y2] = meshgrid(1:sim.fire.matrixsizeY, 1:sim.fire.matrixsizeX);
    % Create the surface plot
    surf(x2, y2, mapping.CoverMap{kp});
    view(0, 90);
    % Enhancements for better visualization
    shading interp; % Interpolate colors across faces and edges
    c = colorbar;   % Display a color bar to indicate the scale of values
    c.Label.String = 'Cover Probability';
    title(['Cover Probability @ \tau=' num2str((kp-1)*sim.fire.dtmins) 'mins']);
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Cover Probability');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sim.runMethod
    % Binary Cover Map and Path Planning
    n=n+1;
    figure(n)
    tl = tiledlayout(5, 3);
    tl.TileSpacing = 'compact';
    tl.Padding = 'compact';
    
    step=sim.updateStep/sim.fire.dt;
    for kp=sim.startEnKF:step:kf-sim.Horizon+1
        nexttile;
        % Find indices of '1's in the binary matrix
        hold on
        tempmap = mapping.CoverMapBinary{kp} - (mapping.CoverMapBinary{kp} & firesimulation.Oanalysis{kp});
        [rows, cols] = find(tempmap);
        % Create a scatter plot of these indices
        scatter(cols, rows, 'filled'); % Use 'filled' for filled markers
    
        [rows, cols] = find(firesimulation.Oanalysis{kp} == 1);
        % Create a scatter plot of these indices
        scatter(cols, rows, 'filled'); % Use 'filled' for filled markers
    
        plot(planning.path{kp}(:,x),planning.path{kp}(:,y));
        
        hold off
        xlim([1, sim.fire.matrixsizeX]); % Adjust x-axis limits to match matrix size
        ylim([1, sim.fire.matrixsizeY]); % Adjust y-axis limits to match matrix size
        % Adding labels and title for clarity
        xlabel('X-axis');
        ylabel('Y-axis');
        zlabel('Value');
        title(['Path Planning @ \tau=' num2str((kp-1)*sim.fire.dtmins) 'mins']);
        legend('Uncertain Points to Cover','Fire Burning Points to Avoid','UAV Path')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    n=n+1;
    figure(n)
    tl = tiledlayout(5, 3);
    tl.TileSpacing = 'compact';
    tl.Padding = 'compact';
    
    step=sim.updateStep/sim.fire.dt;
    for kp=sim.startEnKF:step:kf-sim.Horizon+1
        nexttile;
        % Find indices of '1's in the binary matrix
        hold on
    
        p1=plot(fire.q{kp}(:,x),fire.q{kp}(:,y),'-b',LineWidth=2);
        p2=plot(firesimulation.q_sbar{kp}(:,x),firesimulation.q_sbar{kp}(:,y),'--g',LineWidth=3);
        p3=plot(planning.waypoints{kp}(:,x),planning.waypoints{kp}(:,y));
        legend([p1 p2 p3], 'True Fire',['Future Estimated Fire Propagation @'...
            num2str((kp+sim.Horizon-2)*sim.fire.dtmins) 'mins'],'UAV Path');
        hold off
        xlim([0 sim.fire.gridsizeXaxis]);
        ylim([0 sim.fire.gridsizeYaxis]);
        grid on
        xlabel('X-axis');
        ylabel('Y-axis');
        title(['Path Planning @ \tau=' num2str((kp-1)*sim.fire.dtmins) 'mins']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Propagation Plots with fill
    n=n+1;
    figure(n)
    tl = tiledlayout(5, 3);
    tl.TileSpacing = 'compact';
    tl.Padding = 'compact';
    
    step=sim.updateStep/sim.fire.dt;
    ks=1;
    xCoords = [];
    yCoords = [];
    for kp=sim.startEnKF:step:kf-sim.Horizon+1
        
        kuav = kp*sim.UAViterationsPerFire;
        fovSubset = uav.FOV(ks:kuav);
        ks=kuav;
        % Use cellfun to extract the first row (xCoords) from the subset
        xCell = cellfun(@(c) c(1,:), fovSubset, 'UniformOutput', false);

        % Use cellfun to extract the second row (yCoords) from the subset
        yCell = cellfun(@(c) c(2,:), fovSubset, 'UniformOutput', false);

        % Convert cell arrays to numeric matrices
        xCoordsCurrent = cell2mat(xCell')';
        yCoordsCurrent = cell2mat(yCell')';
        xCoords = [xCoords xCoordsCurrent];
        yCoords = [yCoords yCoordsCurrent];

        nexttile
        hold on
        fill(xCoords, yCoords,'yellow','FaceAlpha',0.5, 'EdgeColor', 'none');
        %fill(xCoordsCurrent, yCoordsCurrent,'magenta','FaceAlpha',0.5, 'EdgeColor', 'none');
        p1=plot(fire.q{kp}(:,x),fire.q{kp}(:,y),'-m',LineWidth=1);
        
        p2=plot(fire.q{kp+sim.propagationHorizon-1}(:,x),fire.q{kp+sim.propagationHorizon-1}(:,y),'-k',LineWidth=2);
        p3=plot(firesimulation.q_pbar{kp}(:,x),firesimulation.q_pbar{kp}(:,y),'-r',LineWidth=1);

        p4=plot(fire.q{kp+sim.Horizon-1}(:,x),fire.q{kp+sim.Horizon-1}(:,y),'-b',LineWidth=2);
        p5=plot(firesimulation.q_sbar{kp}(:,x),firesimulation.q_sbar{kp}(:,y),'-g',LineWidth=1);
        legend([p1 p2 p3 p4 p5], ['True Fire @' num2str((kp-1)*sim.fire.dtmins) 'mins'], ...
            ['True Fire @' num2str((kp+sim.propagationHorizon-2)*sim.fire.dtmins) 'mins'], ...
            ['Fire Propagation @' num2str((kp+sim.propagationHorizon-2)*sim.fire.dtmins) 'mins'], ...
            ['True Fire @' num2str((kp+sim.Horizon-2)*sim.fire.dtmins) 'mins'], ...
            ['Fire Propagation @' num2str((kp+sim.Horizon-2)*sim.fire.dtmins) 'mins']);
        %legend([p1 p2 p3], ['UAV coverage @' num2str((kp-1)*sim.fire.dt) 'mins'],['True Fire @' num2str((kp-1)*sim.fire.dt) 'mins'],['Future Estimated Fire Propagation @' num2str((kp+sim.Horizon-2)*sim.fire.dt) 'mins']);
        hold off
        xlim([0 sim.fire.gridsizeXaxis]);
        ylim([0 sim.fire.gridsizeYaxis]);
        title(['Time ' num2str((kp-1)*sim.fire.dtmins) 'mins'])
        grid on
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%% Plot for paper
kp = 421;
n=n+1;
figure(n)
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 18, 7]);
tl = tiledlayout(1, 2);
tl.TileSpacing = 'tight';
tl.Padding = 'tight';

nexttile;
% Create a mesh grid
[x2, y2] = meshgrid(1:sim.fire.matrixsizeY, 1:sim.fire.matrixsizeX);
% Create the surface plot
surf(x2, y2, mapping.Oburn{kp});
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
% c = colorbar;   % Display a color bar to indicate the scale of values
% c.Label.String = 'Burn Probability';
title(['(a) Burn Probability Map at t=' num2str((kp-1)*sim.fire.dtmins) 'mins']);
xlabel('Columns');
ylabel('Rows');
zlabel('Burn Probability');

%Uncertainty Map
nexttile;
% Create the surface plot
surf(x2, y2, mapping.CoverMap{kp});
view(0, 90);
% Enhancements for better visualization
shading interp; % Interpolate colors across faces and edges
c = colorbar;   % Display a color bar to indicate the scale of values
%c.Label.String = 'Uncertainty';
title(['(b) Uncertainty Map at t=' num2str((kp-1)*sim.fire.dtmins) 'mins']);
xlabel('Columns');
ylabel('Rows');
zlabel('Uncertainty');

%Cover Map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nexttile([1 2]);
% % Find indices of '1's in the binary matrix
% hold on
% tempmap = mapping.CoverMapBinary{kp} - (mapping.CoverMapBinary{kp} & firesimulation.Oanalysis{kp});
% [rows, cols] = find(tempmap);
% % Create a scatter plot of these indices
% scatter(cols, rows, 'filled'); % Use 'filled' for filled markers
% 
% [rows, cols] = find(firesimulation.Oanalysis{kp} == 1);
% % Create a scatter plot of these indices
% %scatter(cols, rows, 'filled'); % Use 'filled' for filled markers
% 
% plot(planning.path{kp}(:,x),planning.path{kp}(:,y));
% 
% % plot(firesimulation.q_abar{kp}(:,x)/sim.fire.gridResolution,firesimulation.q_abar{kp}(:,y)/sim.fire.gridResolution,'-k',LineWidth=1);
% % plot(firesimulation.q_sbar{kp}(:,x)/sim.fire.gridResolution,firesimulation.q_sbar{kp}(:,y)/sim.fire.gridResolution,'-g',LineWidth=1);
% 
% hold off
% xlim([1, sim.fire.matrixsizeX]); % Adjust x-axis limits to match matrix size
% ylim([1, sim.fire.matrixsizeY]); % Adjust y-axis limits to match matrix size
% % Adding labels and title for clarity
% xlabel('Columns');
% ylabel('Rows');
% zlabel('Value');
% title(['(c) Cover Map and Path Planning at t=' num2str((kp-1)*sim.fire.dtmins) 'mins']);
% legend('Waypoints to visit','UAV Trajectory');%,'Current Est. Fire t=70mins','Future Ets. Fire t=100mins')
% grid on



