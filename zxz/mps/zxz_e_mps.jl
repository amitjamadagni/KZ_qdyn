function run_N_den_exp_final_time(psi_init, t_q, time_range)

    N = length(psi_init)
    dt = time_range[2] - time_range[1]

    # psi = deepcopy(psi_init)
    psi = convert_sp_dense(psi_init)

    n_steps = length(time_range)-1

    # println("begin evolution")
    for step=1:n_steps
        evol_gates = evolution_tebd_gates_time_dependent_it(psi, dt, time_range[step], t_q)
        psi = apply(evol_gates, psi; cutoff=1E-8)
        normalize!(psi)
    end
    # println("end evolution")

    return expect_gs(psi)
end

function run_N_den_site_exp_final_time(psi_init, t_q, time_range)

    N = length(psi_init)
    dt = time_range[2] - time_range[1]

    # psi = deepcopy(psi_init)
    psi = convert_sp_dense(psi_init)

    n_steps = length(time_range)-1

    # println("begin evolution")
    for step=1:n_steps
        evol_gates = evolution_tebd_gates_time_dependent_it(psi, dt, time_range[step], t_q)
        psi = apply(evol_gates, psi; cutoff=1E-8)
        normalize!(psi)
    end
    # println("end evolution")

    return site_expect_gs(psi)
end

function run_N_den_exp_time_steps(psi_init, t_q, time_range)#, ntraj)

    N = length(psi_init)
    dt = time_range[2] - time_range[1]

    bd = []
    push!(bd, maxlinkdim(psi_init))

    psi = convert_sp_dense(psi_init)
    d_ts = []
    push!(d_ts, expect_gs(psi))

    n_steps = length(time_range)-1

    # println("begin evolution")
    for step=1:n_steps
        evol_gates = evolution_tebd_gates_time_dependent_it(psi, dt, time_range[step], t_q)
        psi = apply(evol_gates, psi; cutoff=1E-8)
        normalize!(psi)
        push!(bd, maxlinkdim(psi))
        push!(d_ts, expect_gs(psi))
    end
    # println("end evolution")

    return bd, d_ts
end
