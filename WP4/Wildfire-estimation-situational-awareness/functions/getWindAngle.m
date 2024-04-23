function theta = getWindAngle(theta_mean,theta_sigma,sim)
    theta = abs(normrnd(theta_mean,theta_sigma,sim.fire.matrixsizeX,sim.fire.matrixsizeY));
    % need to address the case of rounding the values between 0 and 2pi
end