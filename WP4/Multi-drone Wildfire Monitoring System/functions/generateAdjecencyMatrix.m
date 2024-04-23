function A = generateAdjecencyMatrix(Xpway,Ypway)
    % Generate Adjacency Matrix from the possible waypoints
    N=length(Xpway);
    A = zeros(N,N);
    for i=1:N-1
        A(i,i+1) = norm([Xpway(i) Ypway(i)]-[Xpway(i+1) Ypway(i+1)]);
        A(i+1,i) = A(i,i+1);
    end
    A(1,N) = norm([Xpway(1) Ypway(1)]-[Xpway(N) Ypway(N)]);
    A(N,1) = A(1,N);
    
    % G = digraph(A);    
end