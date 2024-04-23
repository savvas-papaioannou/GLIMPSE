function [Xpway,Ypway,Zpway,Xsgmt,Ysgmt,Xsgmt_plot,Ysgmt_plot,dist_way] = findPossibleWaypoints(Xm,Ym,Xe,Ye,AFOVx,s)
% Finds possible waypoints between measured and estimated fire by
% calculating middle points of the distance segments that are perpendicular
% from the measured perimeter segments and cross the estimated perimeter segments.

% Inputs:   Xm -

% Outputs:  Ym - 
    
    % Calculate line parameters: m-slope and b-constant value of y when x=0, i.e., y = m*x + b -> line equation
    Nm = length(Xm);
    mm = (Ym(1:Nm-1)-Ym(2:Nm))./(Xm(1:Nm-1)-Xm(2:Nm));  % slope
    bm = Ym(1:Nm-1)-mm.*Xm(1:Nm-1);                     % constant
    
    Ne = length(Xe);
    me = (Ye(1:Ne-1)-Ye(2:Ne))./(Xe(1:Ne-1)-Xe(2:Ne));  % slope
    be = Ye(1:Ne-1)-me.*Xe(1:Ne-1);                     % constant

    % find middle point on fire measured segments
    % [Xm(i),Ym(i)]---------[Xa(i),Ya(i)]---------[Xm(i+1),Ym(i+1)]
    Xa = (Xm(1:Nm-1) + Xm(2:Nm))/2;
    Ya = (Ym(1:Nm-1) + Ym(2:Nm))/2;
    Na = length(Xa);
    

    % find the angle of measured segments middle points
    opposite = Ym(1:Nm-1)-Ym(2:Nm);
    adjacent = Xm(1:Nm-1)-Xm(2:Nm);
    th = atan2(opposite,adjacent)-pi/2;
    %thd = atan2d(opposite,adjacent);
    
    % calculate perpendicular segments line equations
    d=5000; % length of segments [m]
    Xaa = Xa-d*cos(th);
    Yaa = Ya-d*sin(th);
    ma = (Ya-Yaa)./(Xa-Xaa);    % slope
    ba = Ya-ma.*Xa;             % constant
    
    % find intersection between measured perimeter segments 
    % and estimate perimeter segments
    Xsgmt = []; Ysgmt = [];
    Xsgmt_plot = []; Ysgmt_plot = [];
    Xpway =[]; Ypway = [];
    dist_way = [];
    
    for i=1:Na
        Xsgmt_temp = []; Ysgmt_temp = [];
        rmv = false;
        dist = Inf;
        for j=1:Ne-1
            Xi = (ba(i)-be(j))/(me(j)-ma(i));
            Yi = me(j)*Xi+be(j);
            if checkIntersection(Xi,Yi,...
                    min(Xe(j),Xe(j+1)),max(Xe(j),Xe(j+1)),...
                    min(Ye(j),Ye(j+1)),max(Ye(j),Ye(j+1)),...
                    min(Xa(i),Xaa(i)),max(Xa(i),Xaa(i)),...
                    min(Ya(i),Yaa(i)),max(Ya(i),Yaa(i)))
                dist_temp = norm([Xa(i) Ya(i)]-[Xi Yi]);
    
                if dist>dist_temp
                    dist = dist_temp;
                    Xsgmt_temp = [Xa(i) Xi];
                    Ysgmt_temp = [Ya(i) Yi];
                end
            end
        end
        if ~isempty(Xsgmt_temp)
            % Xsgmt = [Xsgmt ; Xsgmt_temp];
            % Ysgmt = [Ysgmt ; Ysgmt_temp];
            % Xsgmt_plot = [Xsgmt_plot Xsgmt_temp NaN];
            % Ysgmt_plot = [Ysgmt_plot Ysgmt_temp NaN];
            % dist_way = [dist_way dist];
            for j=1:Nm-1
                Xchk = (ba(i)-bm(j))/(mm(j)-ma(i));
                Ychk = mm(j)*Xchk+bm(j);
                if checkIntersection2(Xchk,Ychk,...
                        min(Xm(j),Xm(j+1)),max(Xm(j),Xm(j+1)),...
                        min(Ym(j),Ym(j+1)),max(Ym(j),Ym(j+1)),...
                        min(Xa(i),Xsgmt_temp(2)),max(Xa(i),Xsgmt_temp(2)),...
                        min(Ya(i),Ysgmt_temp(2)),max(Ya(i),Ysgmt_temp(2)))
                    rmv=true;
                    break
                end
            end

            if ~rmv
                Xsgmt = [Xsgmt ; Xsgmt_temp];
                Ysgmt = [Ysgmt ; Ysgmt_temp];
                Xsgmt_plot = [Xsgmt_plot Xsgmt_temp NaN];
                Ysgmt_plot = [Ysgmt_plot Ysgmt_temp NaN];
                dist_way = [dist_way dist];
            else
                rmv=false;
            end
        end
    end

    % figure(5)
    % clf
    % hold on
    % plot(Xm,Ym,Xe,Ye)
    % plot(Xa,Ya,'o');%,Xaa,Yaa,'+')
    % plot(Xsgmt_plot,Ysgmt_plot,'-*')
    % grid on
    % hold off

    % Remove intersect segments to have possible waypoints in order
    sgmtInt=1;
    while sgmtInt > 0
        N=length(Xsgmt);
        % y = m*x + b -> line equation
        m = (Ysgmt(:,1)-Ysgmt(:,2))./(Xsgmt(:,1)-Xsgmt(:,2));
        b = Ysgmt(:,1)-m.*Xsgmt(:,1);
        Xi = zeros(N,N);
        Yi = zeros(N,N);
        intAdjacencyMatrix = zeros(N,N);
        for i=1:N-1
            for j=i+1:N
                Xi_temp = (b(i)-b(j))/(m(j)-m(i));
                Yi_temp = m(j)*Xi_temp+b(j);
                if checkIntersection(Xi_temp,Yi_temp,...
                        min(Xsgmt(i,1),Xsgmt(i,2)),max(Xsgmt(i,1),Xsgmt(i,2)),...
                        min(Ysgmt(i,1),Ysgmt(i,2)),max(Ysgmt(i,1),Ysgmt(i,2)),...
                        min(Xsgmt(j,1),Xsgmt(j,2)),max(Xsgmt(j,1),Xsgmt(j,2)),...
                        min(Ysgmt(j,1),Ysgmt(j,2)),max(Ysgmt(j,1),Ysgmt(j,2)))
                    Xi(i,j) = Xi_temp;
                    Yi(i,j) = Yi_temp;
                    intAdjacencyMatrix(i,j) = 1;
                    intAdjacencyMatrix(j,i) = 1;
                end
            end
        end
        intPerSegment = sum(intAdjacencyMatrix,2);
        [sgmtInt,indexrmv] = max(intPerSegment);
        if sgmtInt>0
            Xsgmt(indexrmv,:) = [];
            Ysgmt(indexrmv,:) = [];
            dist_way(indexrmv) = [];
        end
    end

    % Calculate possible waypoints 3D coordinates
    N = length(Xsgmt);
    Xsgmt_plot = [];
    Ysgmt_plot = [];
    Xpway = [];
    Ypway = [];
    Zpway = [];
    for i=1:N
        Xsgmt_plot = [Xsgmt_plot Xsgmt(i,:) NaN];
        Ysgmt_plot = [Ysgmt_plot Ysgmt(i,:) NaN];
        Xpway = [Xpway Xsgmt(i,2)];
        Ypway = [Ypway Ysgmt(i,2)];
        FOVx = 2*norm([Xsgmt(i,1) Ysgmt(i,1)]-[Xsgmt(i,2) Ysgmt(i,2)]);
        z_temp = FOVx*s/(2*tan(AFOVx/2));
        Zpway = [Zpway z_temp];
    end


    % figure(6)
    % clf
    % hold on
    % plot(Xm,Ym,Xe,Ye)
    % plot(Xa,Ya,'o');%,Xaa,Yaa,'+')
    % plot(Xsgmt_plot,Ysgmt_plot,'-*')
    % grid on
    % hold off

end
