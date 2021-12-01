# IsingModel_Julia
Small Julia+MC test on the Ising model. Includes a IsingMC.jl 
with all the necessary routines and examples on how to use them.
The main purpose of this code is to get familiar with the 
Monte-Carlo method and Julia.

## Known limitations

- Currently, only the energy is calculated
- No steps are cut from the evaluation of the energy (i.e., 
  during the equilibration phase, leading to larger required 
  simulation "times")
- no parallelization / speed optimization whatsoever


