function fire = fireEvolution(sim,fire,k)
    global x y
    %% Ground truth fire
    for i=1:fire.N
        qX = round(fire.q{k}(i,x));
        qY = round(fire.q{k}(i,y));
        if i==1
            xs = fire.q{k}(i+1,x) - fire.q{k}(fire.N-1,x);
            ys = fire.q{k}(i+1,y) - fire.q{k}(fire.N-1,y);
        elseif i==fire.N
            xs = fire.q{k}(2,x) - fire.q{k}(i-1,x);
            ys = fire.q{k}(2,y) - fire.q{k}(i-1,y);
        else
            xs = fire.q{k}(i+1,x) - fire.q{k}(i-1,x);
            ys = fire.q{k}(i+1,y) - fire.q{k}(i-1,y);
        end
        fire.u{k}(i,:) = [fire.R(qX,qY); fire.U(qX,qY); fire.theta(qX,qY); xs; ys; sim.dt];
        fire.q{k+1}(i,:) = fireStateFcn(fire.q{k}(i,:),fire.u{k}(i,:));
    end

    %% Predicted from historical data
    for i=1:fire.Npre
        qXpre = round(fire.qpre{k}(i,x));
        qYpre = round(fire.qpre{k}(i,y));
        if i==1
            xspre = fire.qpre{k}(i+1,x) - fire.qpre{k}(fire.Npre-1,x);
            yspre = fire.qpre{k}(i+1,y) - fire.qpre{k}(fire.Npre-1,y);
        elseif i==fire.Npre
            xspre = fire.qpre{k}(2,x) - fire.qpre{k}(i-1,x);
            yspre = fire.qpre{k}(2,y) - fire.qpre{k}(i-1,y);
        else
            xspre = fire.qpre{k}(i+1,x) - fire.qpre{k}(i-1,x);
            yspre = fire.qpre{k}(i+1,y) - fire.qpre{k}(i-1,y);
        end
        fire.upre{k}(i,:) = [fire.Rhis(qXpre,qYpre); fire.Uhis(qXpre,qYpre);...
            fire.thetahis(qXpre,qYpre); xspre; yspre; sim.dt];
        fire.qpre{k+1}(i,:) = fireStateFcn(fire.qpre{k}(i,:),fire.upre{k}(i,:));
    end

    %% Estimated from fused data
    if fire.est
        for i=1:fire.Nest
            qXest = round(fire.qest{k}(i,x));
            qYest = round(fire.qest{k}(i,y));
            if i==1
                xsest = fire.qest{k}(i+1,x) - fire.qest{k}(fire.Nest-1,x);
                ysest = fire.qest{k}(i+1,y) - fire.qest{k}(fire.Nest-1,y);
            elseif i==fire.Nest
                xsest = fire.qest{k}(2,x) - fire.qest{k}(i-1,x);
                ysest = fire.qest{k}(2,y) - fire.qest{k}(i-1,y);
            else
                xsest = fire.qest{k}(i+1,x) - fire.qest{k}(i-1,x);
                ysest = fire.qest{k}(i+1,y) - fire.qest{k}(i-1,y);
            end
            switch fire.mode
                case 1
                   fire.uest{k}(i,:) = [fire.Rhis(qXest,qYest); fire.Uhis(qXest,qYest);...
                        fire.thetahis(qXest,qYest); xsest; ysest; sim.dt];
                case 2
                   fire.uest{k}(i,:) = [fire.Rfus(qXest,qYest); fire.Uhis(qXest,qYest);...
                        fire.thetahis(qXest,qYest); xsest; ysest; sim.dt];
                case 3
                   fire.uest{k}(i,:) = [fire.Rhis(qXest,qYest); fire.Ufus(qXest,qYest);...
                        fire.thetafus(qXest,qYest); xsest; ysest; sim.dt];
                case 4
                    fire.uest{k}(i,:) = [fire.Rfus(qXest,qYest); fire.Ufus(qXest,qYest);...
                        fire.thetafus(qXest,qYest); xsest; ysest; sim.dt];
            end
            fire.qest{k+1}(i,:) = fireStateFcn(fire.qest{k}(i,:),fire.uest{k}(i,:));
        end
    end
end