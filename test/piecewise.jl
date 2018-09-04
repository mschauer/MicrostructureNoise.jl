using Test, MicrostructureNoise

@test MicrostructureNoise.piecewise([0,1,2],[10,11]) == ([0, 1, 1, 2], [10, 10, 11, 11])