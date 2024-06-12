empty!(DEPOT_PATH)
push!(DEPOT_PATH, "jlib_1pt6")

include("issh_mps_gen_mc_zo.jl")

tq_arr = [2^i for i=collect(2:14)]

# N = 500
n_traj = 20000

# _a = 2.
# ze = QN("Pf",0,-2) => 1
# oe = QN("Pf",1,-2) => 1
# sites = [Index(ze,oe;tags="Site,S=1/2,n=$n") for n=1:N]
# ener, final_psi = run_xy_zpert(sites, N, _a, 0., 3.)

for param in [parse(Int, ENV["SLURM_ARRAY_TASK_ID"])]
    _tq = tq_arr[param]
    println(_tq)

    f = h5open("quench_dyn/issh_qc_new/data_mps/N_500/wf_data"*"/wf_tq_"*"$(_tq)"*".h5","r")
    final_psi = read(f, "psi", MPS)
    close(f)

    _time = @elapsed mc_den_zo = run_N_den_ssh_ising_final_time(final_psi, n_traj)

    println("den : ", mc_den_zo)
    println("time : ", _time)

end
