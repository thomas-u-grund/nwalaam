import delimited "dev/karate/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "dev/karate/karate_nodes.csv", clear
_nwsyntax karate
mata:
adj = *(`netobj'->get_matrix())
y = st_data(., "faction")
AlaamChangeStats(adj, 0, y, (1,2,6), (0,0,0), (0,0,0), J(0,34,.), D=., z=.)
X = D
	b = J(cols(X), 1, 0)
	ll = -rows(X) * ln(2)
	for (it = 1; it <= 8; it++) {
		p = invlogit(X * b)
		g = X' * (y - p)
		H = cross(X, p :* (1 :- p), X)
		d = invsym(H) * g
		step = 1
		for (k = 1; k <= 40; k++) {
			bn = b + step * d
			lln = sum(y :* ln(invlogit(X * bn)) + (1 :- y) :* ln(invlogit(-X * bn)))
			if (lln >= ll - 1e-12) break
			step = step / 2
		}
		b = bn
		ll = lln
		it, step, ll, b'
	}
end
