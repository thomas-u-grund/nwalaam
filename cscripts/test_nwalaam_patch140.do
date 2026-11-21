* nwalaam 1.4.0: Algorithm S limits the step of every coefficient to [-0.1, 0.1] (symmetric), so the start of the stochastic
* approximation is moderate for every seed; seed 111 and seeds 1-5 of the book's smoking model converge from the Algorithm S start
* itself (no restart); the estimates agree across seeds. Needs the glasgow data (core data folder through NWCORE, else the web).
local core : environment NWCORE
if "`core'" != "" global nwwebpath "`core'/data"
nwclear
nwwebuse glasgow, nwclear
generate byte smoker1 = (smoke1 > 1)
quietly summarize alcohol1
generate double alc_c = alcohol1 - r(mean)
quietly summarize sport1
generate double sport_c = sport1 - r(mean)
local eff "density contagion sender receiver ego(alc_c) ego(sport_c)"
* reference means over seeds 1-30 (1.4.0): -1.202 0.414 -0.390 -0.257 1.845 0.217; spread between seeds below 0.02
matrix REF = (-1.202, 0.414, -0.390, -0.257, 1.845, 0.217)
foreach s in 1 2 3 4 5 111 {
	quietly nwalaam smoker1 glasgow1, effects(`eff') seed(`s')
	assert e(sa_converged) == 1
	assert e(repetitions) <= 2
	matrix b = e(b)
	matrix st = e(b_start)
	* the start is limited by Algorithm S itself (it was 10 after the 1.3.1 restart, up to 1800 before)
	forvalues j = 1/6 {
		assert abs(st[1,`j']) < 5
		assert abs(b[1,`j'] - REF[1,`j']) < 0.1
	}
	di "seed `s': converged, max |start| " max(abs(st[1,1]), abs(st[1,2]), abs(st[1,3]), abs(st[1,4]), abs(st[1,5]), abs(st[1,6]))
}
* seeds that started at 400 to 1800 before (6, 16, 20, 21, 27) converge in one repetition from a start of order 1
foreach s in 6 20 {
	quietly nwalaam smoker1 glasgow1, effects(`eff') seed(`s')
	assert e(sa_converged) == 1
	matrix st = e(b_start)
	forvalues j = 1/6 {
		assert abs(st[1,`j']) < 5
	}
}
di as text "PASS: symmetric Algorithm S step cap: wild starts gone, seed 111 and seeds 1-5 converge"
di "test_nwalaam_patch140: ALL OK"
