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
    boltzmannConstant =  8.617333262e10-5
    # Also, why not set it to one?
#     boltzmannConstant = 1.0

    using Statistics
#     using Plots
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

    function getLocalHamiltonian(sim::MCSimulation, pointX, pointY)
       thisPoint = sim.lattice[pointX,pointY]
       localHamiltonian = 0.0
        if pointX == 1
            localHamiltonian += sim.lattice[sim.latticeSize, pointY]*thisPoint
        else
            localHamiltonian += sim.lattice[pointX-1, pointY]*thisPoint
        end
        if pointY == 1
            localHamiltonian += sim.lattice[pointX, sim.latticeSize]*thisPoint
        else
            localHamiltonian += sim.lattice[pointX, pointY-1]*thisPoint
        end
        if pointX == sim.latticeSize
            localHamiltonian += sim.lattice[1, pointY]*thisPoint
        else
            localHamiltonian += sim.lattice[pointX+1, pointY]*thisPoint
        end
        if pointY == sim.latticeSize
            localHamiltonian += sim.lattice[pointX, 1]*thisPoint
        else
            localHamiltonian += sim.lattice[pointX, pointY+1]*thisPoint
        end
        return localHamiltonian
    end

    # Time evolve one lattice.
    function timeEvolve(sim::MCSimulation)
        averagedEnergy = 0.0
        energy = energyFromLattice(sim)
        acceptedSteps = 1
        for step in 1:sim.stepsToEvolve
            # Determine which point in the lattice may be flipped.
            posToFlip = rand(0:((sim.latticeSize*sim.latticeSize)-1))
            xToFlip = posToFlip ÷ sim.latticeSize+1
            yToFlip = posToFlip % sim.latticeSize+1
            sim.lattice[xToFlip, yToFlip] *= -1

            # Calculate change due to flip.
            deltaE = -0.5*sim.interactionStrength*getLocalHamiltonian(sim, xToFlip, yToFlip)

            # Determine whether the change should be accepted.
            accepted = true
            if deltaE > 0.0
                randomNumber = rand()
                probability = exp(Float64(-1.0*deltaE/(boltzmannConstant*sim.temperatureK)))
                if probability < randomNumber
                    accepted = false
                end
            end

            # If accepted, update the energy and
            if accepted == true
                energy = energyFromLattice(sim)
                acceptedSteps += 1
                averagedEnergy = ((averagedEnergy*(acceptedSteps-1))+energy)/acceptedSteps
#                 println("Accepted step, energy is now: ", energy)
            else
               sim.lattice[xToFlip, yToFlip] *= -1
            end
        end
        return averagedEnergy
    end

    function energyFromLattice(sim::MCSimulation)
        energy = 0.0
        for i in 1:sim.latticeSize
            for j in 1:sim.latticeSize
                # We assume periodic boundary conditions.
                energy += getLocalHamiltonian(sim, i, j)
            end
        end
        return -0.5*energy*sim.interactionStrength
    end


    ####################
    # Tests
    ####################

#     newSim = MCSimulation(2.0/boltzmannConstant, 20, 50000, 4.0)
#     initialize(newSim)
#     println(newSim)
#     timeEvolve(newSim)
#     println(newSim)
#     visualizeSim(newSim)
    temperatures = range(0.5, 10.0, step=0.5)
    for temp in temperatures
        energies = []
        newSim = MCSimulation(temp/boltzmannConstant, 20, 50000, 4.0)
        initialize(newSim)
        energy = timeEvolve(newSim)
        append!(energies, energy)
        println(temp, " ", energy)
    end
#     println(energies)
end
