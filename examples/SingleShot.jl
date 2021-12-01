#=
SingleShot: Performs and visualizes one MC simulation at a specific
temperature.
- Julia version: 1.4.1
- Author: Lenz Fiedler
- Date: 2021-12-01
=#
include("../src/IsingMC.jl")
using .IsingMC

newSim = IsingMC.MCSimulation(2.0/IsingMC.boltzmannConstant, 20, 10000, 4.0)
IsingMC.initialize(newSim, "positive")
IsingMC.timeEvolve(newSim, true)
IsingMC.visualizeSim(newSim)
