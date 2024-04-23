function U = getWindMagnitude(U_mean,U_sigma,sim)
    U = abs(normrnd(U_mean,U_sigma,sim.fire.matrixsizeX,sim.fire.matrixsizeY));
end