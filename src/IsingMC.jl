"""
# module IsingMC.jl

- Julia version: 1.4.1
- Author: Lenz Fiedler
- Date: 2021-11-19

# Examples

See examples subfolder.

"""


module IsingMC
    using Statistics
    using Plots

    # Boltzmann constant in atomic units.
    boltzmannConstant =  8.617333262e10-5

    mutable struct MCSimulation
       temperatureK::Float64
       latticeSize::UInt64
       interactionStrength::Float64
       lattice
    end
    """
    Create a new MC simulation. The lattice array is created automatically.

    - `temperatureK::Float64`: Temperature of the simulation.
    - `latticeSize::UInt64`: Size of the lattice (in either direction, a quadratic lattice is
                             assumed.
    - `interactionStrength::Float64`: Strength of the spin interaction (for energy calculation)
    """
    MCSimulation(temperatureK::Float64, latticeSize::UInt64, interactionStrength::Float64) =
            MCSimulation(temperatureK, latticeSize, interactionStrength,
                         zeros(Int64, latticeSize, latticeSize))

    """
    Visualize the current state of a simulation. This is probably not the best way to do this,
    but it works for now.

    - `sim::MCSimulation`: The simulation to be visualized.
    """
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


    """
    Initialize a simulation.

    - `sim::MCSimulation`: The simulation to be initialized.
    - `initType::string`: Type of initialization to be performed. Default "random", assgning
                           spins at random. "positive" or "negative" initializes the lattice
                           entirely with positive or negative spins, respectively.
    """
    function initialize(sim::MCSimulation, initType::string="random")
        for i in 1:sim.latticeSize
            for j in 1:sim.latticeSize
                if initType == "random"
                    randomNumber = rand(0:1)
                    if randomNumber == 0
                       randomNumber -= 1
                    end
                    sim.lattice[i,j] = randomNumber
                elseif initType == "negative"
                   sim.lattice[i,j] = -1
                elseif initType == "positive"
                    sim.lattice[i,j] = 1
                else
                    throw(DomainError(initType, "Invalid initilization type"))
                end
            end
        end
    end

    """
    Calculate the local Hamiltonian around a lattice point.

    - `sim::MCSimulation`: The simulation which provides the lattice.
    - `pointX::Int`: Lattice point X coordinate.
    - `pointY::Int`: Lattice point Y coordinate.
    """
    function getLocalHamiltonian(sim::MCSimulation, pointX::Int, pointY::Int)
       thisPoint = sim.lattice[pointX,pointY]
       localHamiltonian = 0.0
       # We assume periodic boundary conditions.
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

    """
    For a given simulation, perform the "time" evolution.

    - `sim::MCSimulation`: The simulation for which the time evolution will be performed.
    - `stepsToEvolve::UInt644`: Time steps this simulation is being evolved for.
    - `printEnergies::bool`: If true, energies are printed during the simulation.
    """
    function timeEvolve(sim::MCSimulation, stepsToEvolve, printEnergies::bool=false)
        averagedEnergy = 0.0
        energy = energyFromLattice(sim)
        averagedEnergy = energy
        acceptedSteps = 1
        for step in 1:stepsToEvolve
            # Determine which point in the lattice may be flipped.
            posToFlip = rand(0:((sim.latticeSize*sim.latticeSize)-1))
            xToFlip = posToFlip ÷ sim.latticeSize+1
            yToFlip = posToFlip % sim.latticeSize+1
            sim.lattice[xToFlip, yToFlip] *= -1

            # Calculate change due to flip.
            newEnergy = energyFromLattice(sim)
            deltaE = newEnergy-energy

            # Determine whether the change should be accepted.
            accepted = true
            if deltaE > 0.0
                randomNumber = rand()
                probability = exp(Float64(-1.0*deltaE/(boltzmannConstant*sim.temperatureK)))
                if probability < randomNumber
                    accepted = false
                end
            end

            # If accepted, update the energy.
            if accepted == true
                energy = newEnergy
                acceptedSteps += 1
                averagedEnergy = ((averagedEnergy*(acceptedSteps-1))+energy)/acceptedSteps
                if printEnergies == true
                    println("Accepted step, energy is now: ", energy)
                end
            else
               sim.lattice[xToFlip, yToFlip] *= -1
            end
        end
        return averagedEnergy
    end

    """
    Calculate the current total energy for a given simulation.

    # Arguments

    - `sim::MCSimulation`: Simulation for which to calculate the current energy.

    """
    function energyFromLattice(sim::MCSimulation)
        energy = 0.0
        for i in 1:sim.latticeSize
            for j in 1:sim.latticeSize
                energy += getLocalHamiltonian(sim, i, j)
            end
        end
        return -0.5*energy*sim.interactionStrength
    end
end
