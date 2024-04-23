function fusion = makeFusion(fusion,fov,fire,sim)
    x=1;y=2;
    fov_scaled = round(fov(1:2,:)./sim.fire.gridResolution);
    %The correct ones
    indX = [min(fov_scaled(y,:)) max(fov_scaled(y,:))-1];
    indY = [min(fov_scaled(x,:)) max(fov_scaled(x,:))-1];

    fusion.R_mean(indX,indY) = fire.R(indX,indY);
    fusion.R_sigma(indX,indY) = 0;

    fusion.U_mean{1}(indX,indY) = fire.U{1}(indX,indY);
    fusion.U_sigma{1}(indX,indY) = 0;

    fusion.theta_mean{1}(indX,indY) = fire.theta{1}(indX,indY);
    fusion.theta_sigma{1}(indX,indY) = 0;

end