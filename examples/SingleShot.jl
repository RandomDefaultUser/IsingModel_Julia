#=
SingleShot: Performs and visualizes one MC simulation at a specific
temperature.
- Julia version: 1.4.1
- Author: Lenz Fiedler
- Date: 2021-12-01
=#
include("../src/IsingMC.jl")
using .IsingMC

newSim = IsingMC.MCSimulation(7.0/IsingMC.boltzmannConstant, 20, 4.0)
IsingMC.initialize(newSim, "negative")
IsingMC.timeEvolve(newSim, 5000, true)
IsingMC.performSimulation(newSim)
