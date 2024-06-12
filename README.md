# Dynamics of Kibble-Zurek mechanism using ITensors.jl
Kibble-Zurek mechanism and errors of gapped phases

This repo has the basic code for the paper: https://arxiv.org/abs/2401.13625

### Installing the required packages
1. Get julia 1.6.6
2. Use the `Project.toml` and `Manifest.toml` to generate the suitable work environments, using the link in `reproducing_environment`

### Dynamics of the different models
General structure for 
1. `ising` -> Ising
2. `ssh` -> SSH model
3. `zxz` -> Cluster state model
4. `issh_qc` -> Extended SSH model dynamics on a quantum computer

Each of the above repo has the following structure:
`gen` -> DMRG code (might also include the dynamics)
`mps` -> MPS sampling (might also include the dynamics)
`data_mps` -> sample run for a given system size of N=500 on a HPC node (can be extended to other system sizes as well)



