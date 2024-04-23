function makeanimation(sim,fire,uav,N)
    clfall
    plt = initializePlot(sim,fire,uav);
    for k=1:N
        plotData(plt,sim,fire,uav,k)
        % if k>1500
        %     pause(0.1)
        % end
    end
end

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