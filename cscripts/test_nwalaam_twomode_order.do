* A two-mode network whose nodes are stored in label order (the levels interleaved, as the Davis southern women network is loaded) is read with
* the right actors: the outcome level is taken from the modes of the nodes, not from the position of the nodes.
nwwebuse davis, nwclear
generate byte mrs = (title == "Mrs") if _nwmode == "1"
nwalaam mrs davis, effects(density activity) mode(1) seed(1)
assert e(N) == 18 & e(n_ones) == 7 & e(nodes) == 32
di "PASS: two-mode outcome level taken from the node modes"
