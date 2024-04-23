function fire = fireModel(sim,fire,k)
    x=1; y=2;

    %% Forward Model
    qD = ceil(fire.q{k-1}/sim.fire.gridResolution);                       % discretization of cordinates
    xs = [fire.q{k-1}(2,x);fire.q{k-1}(3:end,x);fire.q{k-1}(2,x)] -...        % x component differential
        [fire.q{k-1}(end-1,x);fire.q{k-1}(1:end-2,x);fire.q{k-1}(end-1,x)];
    ys = [fire.q{k-1}(2,y);fire.q{k-1}(3:end,y);fire.q{k-1}(2,y)] -...        % y component differential
        [fire.q{k-1}(end-1,y);fire.q{k-1}(1:end-2,y);fire.q{k-1}(end-1,y)];

    % xindex = sum(find(qD(:,x)>sim.fire.matrixsizeX));
    % yindex = sum(find(qD(:,y)>sim.fire.matrixsizeY));
    % if xindex>0 || yindex>0
    %     qD(:,x)
    % end

    idx = sub2ind(size(fire.R),qD(:,y),qD(:,x));                        % subscripts convertion to linear indices

    
    fire.u{k-1} = [fire.R(idx) fire.U{sim.fire.weather.idx(k-1)}(idx) fire.theta{sim.fire.weather.idx(k-1)}(idx) xs ys];        % input matrix

    fire.q{k} = fire.q{k-1} + sim.fire.dt.*(Qdif(fire.u{k-1})./(sim.minTosec));                  % firefront propagation
    fire.N(k) =  fire.N(k-1);

    %% Clipping
    % [X,Y,N] = clipping_fcn(fire.q{k}(:,x),fire.q{k}(:,y),...
    %     fire.N(k),sim.fire.theta_mean(sim.fire.weather.idx(k-1)));
    % if N~=-1
    %     fire.N(k) = N;
    %     fire.q{k} = [];
    %     fire.q{k} = [X Y];
    % end
    % 
    q_temp = clipping_fcn_Matlab(fire.q{k});
    N = size(q_temp,1);
    if N~=fire.N(k)
        fire.N(k) = N;
        fire.q{k} = [];
        fire.q{k} =  q_temp;
    end

    % Regriding
    [X,Y,N] = regriding(fire.q{k}(:,x),fire.q{k}(:,y),fire.N(k),sim.fire.T);
    if N~=fire.N(k)
        fire.N(k) = N;
        fire.q{k} = [];
        fire.q{k} = [X Y];
    end

    % Calculate distance
    distances = [vecnorm(diff(fire.q{k}), 2, 2); norm(fire.q{k}(1, :) - fire.q{k}(fire.N(k), :))];
    fire.perimeterSize(k) = sum(distances);
    
    % Calculate the area using the shoelace formula
    x = fire.q{k}(:,1);
    y = fire.q{k}(:,2);

    fire.areaSize(k) = 0.5 * abs(sum(x(1:end - 1) .* y(2:end)) + x(end) * y(1) - sum(x(2:end) .* y(1:end - 1)) - x(1) * y(end));

end

function Q = Qdif(u)
    % Evaluate differentials
    R = u(:,1);
    U = u(:,2);
    theta = u(:,3);
    xs = u(:,4);
    ys = u(:,5);
    
    LB = 0.936.*exp(0.2566*U) + 0.461.*exp(-0.1548*U) - 0.397;
    HB = (LB + (LB.^2 - 1).^0.5) ./ (LB - (LB.^2 - 1).^0.5);
    
    a = (0.5*(R+(R./HB)))./LB;
    b = (R+(R./HB))./2;
    c = b - R./HB;
    
    lowPart = sqrt(b.^2.*(xs.*cos(theta)-ys.*sin(theta)).^2 + a.^2.*(xs.*sin(theta)+ys.*cos(theta)).^2);
    xupPart = a.^2.*cos(theta).*(xs.*sin(theta)+ys.*cos(theta)) - b.^2.*sin(theta).*(xs.*cos(theta)-ys.*sin(theta));
    yupPart = -a.^2.*sin(theta).*(xs.*sin(theta)+ys.*cos(theta)) - b.^2.*cos(theta).*(xs.*cos(theta)-ys.*sin(theta));
    
    xRatio = xupPart./lowPart;
    xRatio(isinf(xRatio)) = 0;
    yRatio = yupPart./lowPart;
    yRatio(isinf(yRatio)) = 0;
    
    Xt = xRatio + c.*sin(theta);
    Yt = yRatio + c.*cos(theta);
    Q = [Xt Yt];
end