import delimited "dev/karate/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "dev/karate/karate_nodes.csv", clear
nwalaam faction karate, effects(density activity contagion) runs(5) iterations(20000) seed(1)
matrix list e(obs_stats)
matrix list e(b_start)
