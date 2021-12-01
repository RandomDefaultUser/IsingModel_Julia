#=
TemperatureComparison: Runs a MC simulation for multiple temperatures and gives
the observables (currently only energy) for them.
- Julia version: 1.4.1
- Author: Lenz Fiedler
- Date: 2021-12-01
=#
include("../src/IsingMC.jl")
using .IsingMC

temperatures = range(0.5, 10.0, step=0.5)
for temp in temperatures
    energies = []
    newSim = IsingMC.MCSimulation(temp/IsingMC.boltzmannConstant, 20, 4.0)
    IsingMC.initialize(newSim, "negative")
    energy = IsingMC.performSimulation(newSim, 5000)
    append!(energies, energy)
    println(temp, " ", energy)
end
