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
    # I don't want to deal with unitful constants right now
    boltzmannConstant = 1.380649e-23
    using Plots
        ####################
    # Definitions
    ####################
    # Holds a single MC simulation.
    # I realize latticeSize should never be negative, but setting it as UInt prints it as hexad
    mutable struct MCSimulation
       temperatureK::Float64
       latticeSize::UInt64
       stepsToEvolve::UInt64
       interactionStrength::Float64
       lattice
    end
    MCSimulation(temperatureK, latticeSize, stepsToEvolve, interactionStrength) =
            MCSimulation(temperatureK, latticeSize, stepsToEvolve, interactionStrength,
                         zeros(Int64, latticeSize, latticeSize))

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

    # Time evolve one lattice.
    function timeEvolve(sim::MCSimulation)
        lastEnergy = energyFromLattice(sim)
        for step in 1:sim.stepsToEvolve
            # Determine which point in the lattice may be flipped.
            posToFlip = rand(0:(sim.latticeSize*sim.latticeSize))
            xToFlip = posToFlip ÷ sim.latticeSize+1
            yToFlip = posToFlip % sim.latticeSize+1
            sim.lattice[xToFlip, yToFlip] *= -1
            deltaE = energyFromLattice(sim) - lastEnergy
#             println(xToFlip, " ", yToFlip, " ", deltaE, " ")
            if deltaE > 0.0
                randomNumber = rand()
                probability = exp(Float64(deltaE/(boltzmannConstant*sim.temperatureK)))
                if probability >= randomNumber
                    sim.lattice[xToFlip, yToFlip] *= -1
                end
            end
#             println(sim.lattice)
        end
    end

    function energyFromLattice(sim::MCSimulation)
        energy = 0.0
        for i in 1:sim.latticeSize
            for j in 1:sim.latticeSize
                # We assume periodic boundary conditions.
                localHamiltonian = 0.0
                thisPoint = sim.lattice[i,j]
                if i == 1
                    localHamiltonian += sim.lattice[sim.latticeSize, j]*thisPoint
                else
                    localHamiltonian += sim.lattice[i-1, j]*thisPoint
                end
                if j == 1
                    localHamiltonian += sim.lattice[i, sim.latticeSize]*thisPoint
                else
                    localHamiltonian += sim.lattice[i, j-1]*thisPoint
                end
                if i == sim.latticeSize
                    localHamiltonian += sim.lattice[1, j]*thisPoint
                else
                    localHamiltonian += sim.lattice[i+1, j]*thisPoint
                end
                if j == sim.latticeSize
                    localHamiltonian += sim.lattice[i, 1]*thisPoint
                else
                    localHamiltonian += sim.lattice[i, j+1]*thisPoint
                end
                energy += localHamiltonian
            end
        end
        return -0.5*energy*sim.interactionStrength
    end


    ####################
    # Tests
    ####################

    newSim = MCSimulation(0.1, 20, 100000, 1000.0)
    initialize(newSim)
    println(newSim)
    timeEvolve(newSim)
    println(newSim)
    visualizeSim(newSim)
end
