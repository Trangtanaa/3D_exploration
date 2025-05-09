function G=Graph(pop,rc)
N=size(pop,1);
adj_pop=zeros(N);
for i=1:N
    for j=1:N
        dist = sqrt((pop(i,1)-pop(j,1))^2+(pop(i,2)-pop(j,2))^2+(pop(i,3)-pop(j,3))^2);
        if dist <= rc && dist ~=0
            adj_pop(i,j)=dist;
        end
    end
end
G=graph(adj_pop);
