function sim = updateParams(sim, k)
    if k>1 && k > sim.fire.weather.Change(sim.fire.weather.idx(k-1) + 1)
        sim.fire.weather.idx(k) = sim.fire.weather.idx(k-1) + 1;
    else
        sim.fire.weather.idx(k) = sim.fire.weather.idx(k-1);
    end
end