function [waypoints,waypointsSize,path] = getDesirePath(A,start_node,target_node,Xpway,Ypway,Zpway)
    [~,path] = dijkstra(A,start_node,target_node);
    if length(path(path==start_node)) > 1
        ind = find(path==start_node);
        path(ind(1))=1;
    end
    path = flip(path);
    waypoints = [Xpway(path); Ypway(path); Zpway(path)];
    waypointsSize = size(waypoints,2);
end