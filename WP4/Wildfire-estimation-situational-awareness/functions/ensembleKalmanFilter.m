function [x_fbar,x_abar,x_ensemble] = ensembleKalmanFilter(N_e,fireestimation,y,Sigma,k)
    % Forecast step: propagate each ensemble member
    % for i = 1:N_e 
    %     x_f(:,i) = f(x_ensemble(:,i)) + sqrt(Q) * randn(n, 1);
    %     y_f(:,i) = h(x_f(:,i));
    %     yv(:,i) = y + sqrt(R)*randn(p,1);
    % end

    % tranformation to vector of the form X = [(x1,y1),...,(xN,yN)];
    
    % enkf setup
    %x_ensemble = repmat(x_ini, 1, N_e) + sqrt(Q) * randn(n, 1);
    sz = size(fireestimation{1}.q{k});
    C = eye(sz(1)*sz(2));
    R = diag(ones(1,sz(1)*sz(2))*Sigma^2);
    
    y_res = reshape(y',[],1);

    for i = 1:N_e
        x_f(:,i) = reshape(fireestimation{i}.q{k}',[],1);
        y_f(:,i) = C*x_f(:,i);
        yv(:,i) = y_res + normrnd(0,Sigma,size(y_res));
    end
    x_fbar=mean(x_f,2);    

    X = x_f-x_fbar;
    P_f=X*X'/(N_e-1);


    K = P_f * C' / (C * P_f * C'+ R); % Kalman gain

    for i = 1:N_e
        x_ensemble(:,i)=x_f(:,i)+K*(yv(:,i)-y_f(:,i));
    end
    x_abar=mean(x_ensemble,2);
end


    % % Covariance update
    % P = X * X' / (N - 1);
    % K = P * H' / (H * P * H' + R); % Kalman gain
    % 
    % % Update each ensemble member
    % for i = 1:N
    %     innovation = z + sqrt(R) * randn(stateDim, 1) - H * ensemble(:, i);
    %     ensemble(:, i) = ensemble(:, i) + K * innovation;
    % end