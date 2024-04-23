function fireestimation = updateFireEstimationParameters(fireestimation,fusion,sim)
    fireestimation.U = struct([]);
    fireestimation.theta = struct([]);
    for i=1:sim.fire.weather.numberofchange
        fireestimation.U{i} = abs(normrnd(fusion.U_mean{i},fusion.U_sigma{i}));
        fireestimation.theta{i} = abs(normrnd(fusion.theta_mean{i},fusion.theta_sigma{i}));
    end
    fireestimation.R = abs(normrnd(fusion.R_mean,fusion.R_sigma));
end