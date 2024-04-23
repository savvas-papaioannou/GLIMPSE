function makeanimation(sim,fire,uav,N)
    plt = initializePlot(sim,fire,uav);
    kf=sim.startEnKF;
    counter = 0;

    for k=sim.start:N
        counter = counter + 1;
        if counter >= sim.UAViterationsPerFire
            kf = kf+1;
            counter=0;
        end
        plotData(plt,sim,fire,uav,k,kf)
        % if k>1500
             %pause(0.02)
        % end
    end
end