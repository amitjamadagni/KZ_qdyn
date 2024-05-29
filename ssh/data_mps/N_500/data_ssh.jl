empty!(DEPOT_PATH)
push!(DEPOT_PATH, "/dss/dsshome1/07/di93hog/jlib_1pt6")

include("/dss/dsshome1/07/di93hog/quench_dyn/ssh/gen/ssh_mps_gen_new_minus.jl")
include("/dss/dsshome1/07/di93hog/quench_dyn/ssh/mps/ssh_e_mps.jl")

tq_arr = [2^i for i=collect(2:14)]

N = 500
# n_traj = 20000

_a = 5.
ze = QN("Pf",0,-2) => 1
oe = QN("Pf",1,-2) => 1
sites = [Index(ze,oe;tags="Site,S=1/2,n=$n") for n=1:N]
ener, final_psi = run_ssh(sites, N, _a, 1.)


for param in [parse(Int, ENV["SLURM_ARRAY_TASK_ID"])]
    _tq = tq_arr[param]
    println(_tq)

    tr_ = collect(-_a*_tq:0.05:0.)

    # den = run_N_den_exp_minus(deepcopy(final_psi), _tq, tr_)#, n_traj)
    den = run_N_den_exp_minus_final_time(deepcopy(final_psi), _tq, tr_)
    println(" den : ", den)
end
