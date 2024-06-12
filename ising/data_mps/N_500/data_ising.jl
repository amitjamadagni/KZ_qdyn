empty!(DEPOT_PATH)
push!(DEPOT_PATH, "jlib_1pt6")

include("ising_gen_mps_new.jl")
include("ising_e_mps.jl")


tq_arr = [2^i for i=collect(2:14)]

N = 500
# n_traj = 20000

_a = 5.
sites = siteinds("S=1/2", N)
ener, psi_init = run_ising_xpert(1., sites, N, _a)


for param in [parse(Int, ENV["SLURM_ARRAY_TASK_ID"])]
    _tq = tq_arr[param]
    println(_tq)

    tr_ = collect(-_a*_tq:0.05:0.)

    # den = run_N_den_final_time(deepcopy(psi_init), _tq, tr_, n_traj)
    den = run_N_den_exp_final_time(deepcopy(psi_init), _tq, tr_)
    println(" den : ", den)
end
