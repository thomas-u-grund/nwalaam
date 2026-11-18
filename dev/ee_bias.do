* Is the difference between the EE estimate and the exact maximum likelihood estimate (karate club; exact: Density 0.535, Activity -0.832,
* Contagion 1.461) a matter of the learning rate or of the number of runs? EE with three learning rates, 30 runs each.
import delimited "dev/karate/karate_edges.csv", clear
nwfromedge from to, undirected name(karate)
import delimited "dev/karate/karate_nodes.csv", clear
foreach lr in 0.05 0.01 0.002 {
	local it = cond(`lr' == 0.002, 250000, 100000)
	nwalaam faction karate, effects(density activity contagion) method(ee) runs(30) iterations(`it') lrate(`lr') seed(5) noci
	matrix b = e(b)
	di "RESULT lrate `lr': " %7.3f b[1,1] " " %7.3f b[1,2] " " %7.3f b[1,3]
}
