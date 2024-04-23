function [sdi,adi,jc,ss] = getMetrics(qF,qS,sim)
    %% Inputs
    % qF - points of true fire
    % qS - points of simulated fire

    gridSizeX = sim.fire.gridsizeXaxis;
    gridSizeY = sim.fire.gridsizeYaxis;
    gridResolution = sim.fire.gridResolutionMetrics;

    occupancyGridF = getBinaryOccupancyMatrix(qF,gridSizeX,gridSizeY,gridResolution);
    occupancyGridS = getBinaryOccupancyMatrix(qS,gridSizeX,gridSizeY,gridResolution);
    
    % % Convert them to (x, y) coordinates
    % [yF, xF] = find(occupancyGridF);
    % [yS, xS] = find(occupancyGridS);
    % 
    % % Create a new figure
    % figure(2);
    % 
    % % Plot the first matrix
    % scatter(xF, yF, 'r'); % Red color for matrix1
    % hold on; % Retain the current plot when adding new plots
    % 
    % 
    % % Overlay the third matrix
    % scatter(xS, yS, 'b'); % Blue color for matrix3
    % 
    % % Overlay the second matrix
    % % scatter(x2, y2, 'g'); % Green color for matrix2
    % 
    % % Adjust the plot
    % % legend('true fire','analysis', 'forecast');
    % title('Overlay of Occupancy Matrices');
    % xlabel('X-axis');
    % ylabel('Y-axis');
    % axis equal; % To keep the x and y scale the same
    % grid on; % Optional, for a grid

    sdi = getSDI(occupancyGridF,occupancyGridS);
    adi = getADI(occupancyGridF,occupancyGridS);
    jc = getJaccardCoefficient(occupancyGridF,occupancyGridS);
    ss = getSorensenSmilarity(occupancyGridF,occupancyGridS);
end

function sdi = getSDI(F,S)
    I = and(S,F);
    OE = sum(S-I,'all');
    UE = sum(F-I,'all');
    sdi = (OE+UE)/sum(F,'all');
end

function adi = getADI(F,S)
    I = and(S,F);
    OE = sum(S-I,'all');
    UE = sum(F-I,'all');
    adi = (OE+UE)/sum(I,'all');
end

function jc = getJaccardCoefficient(F,S)
    I = sum(and(S,F),'all');
    jc = (I)/((sum(F,'all')+sum(S,'all'))-I);
end

function ss = getSorensenSmilarity(F,S)
    I = sum(and(S,F),'all');
    ss = (2*I)/(sum(F,'all')+sum(S,'all'));
end