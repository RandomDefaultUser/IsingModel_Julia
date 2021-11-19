"""
# module IsingMC.jl

- Julia version: 1.4.1
- Author: Lenz Fiedler
- Date: 2021-11-19

# Examples

```jldoctest
test
julia>
```
"""


module IsingMC
    using Plots

    ####################
    # Definitions
    ####################
    # Holds a single MC simulation.
    # I realize latticeSize should never be negative, but setting it as UInt prints it as hexad
    mutable struct MCSimulation
       temperatureK::Float64
       latticeSize::UInt64
       lattice
    end
    MCSimulation(temperatureK, latticeSize) = MCSimulation(temperatureK, latticeSize, zeros(Int64, latticeSize, latticeSize))

    # Initialize the lattice of a simulation.
    # I know this is inefficient still.
    function initialize(sim::MCSimulation)
        for i in 1:sim.latticeSize
            for j in 1:sim.latticeSize
                randomNumber = rand(0:1)
                if randomNumber == 0
                   randomNumber -= 1
                end
                sim.lattice[i,j] = randomNumber
            end
        end
    end

    # Visualize the current state of a simulation.
    # This is probably not the best way to do this.
    function visualizeSim(sim::MCSimulation)
        plusx = []
        plusy = []
        minusx = []
        minusy = []
        for i in 1:sim.latticeSize
            for j in 1:sim.latticeSize
                if sim.lattice[i,j] == 1
                    append!(plusx, i)
                    append!(plusy, j)
                else
                    append!(minusx, i)
                    append!(minusy, j)
                end

            end
        end
        plot(plusx, plusy, seriestype = :scatter, label = "plus", markershape = :square,
             markersize = 7, markercolor = :blue, markerstrokecolor=:blue)
        display(plot!(minusx, minusy, seriestype = :scatter, label = "minus",
                markershape = :square, markersize = 7, markercolor = :red,
                markerstrokecolor=:red))
    end


    ####################
    # Tests
    ####################

    newSim = MCSimulation(1, 20)
    initialize(newSim)
    visualizeSim(newSim)
    println(newSim)
end
