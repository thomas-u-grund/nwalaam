import delimited "dev/karate/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "dev/karate/karate_nodes.csv", clear
_nwsyntax karate
mata:
adj = *(`netobj'->get_matrix())
y = st_data(., "faction")
AlaamChangeStats(adj, 0, y, (1,2,6), (0,0,0), (0,0,0), J(0,34,.), D=., z=.)
X = D
b = J(3,1,0)
p = invlogit(X*b)
g = X' * (y - p)
H = cross(X, p :* (1 :- p), X)
rank(H)
H
d = invsym(H)*g
d
bn = b + d
sum(y :* ln(invlogit(X * bn)) + (1 :- y) :* ln(invlogit(-X * bn)))
end
