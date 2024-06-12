empty!(DEPOT_PATH)
push!(DEPOT_PATH, "jlib_1pt6")

include("zxz_mps_gen_exp.jl")
include("zxz_e_mps.jl")

tq_arr = [2^i for i=collect(2:14)]

N = 500
# n_traj = 20000

hx_init = 5.
sites = siteinds("S=1/2", N)
ener, final_psi = run_zxz_pert(sites, N, 1., hx_init)

for param in [parse(Int, ENV["SLURM_ARRAY_TASK_ID"])]
    _tq = tq_arr[param]
    # println(_tq)

    tr_ = collect(-hx_init*_tq:0.05:0.)

    den = run_N_den_exp_final_time(deepcopy(final_psi), _tq, tr_)
    println(" den : ", den)
end
