function mean_var(data)
    return sum(data)/float(length(data))
end

function variance_y(data)
    mv = mean_var(data)
    s = 0
    for elem in data
        s = s + (elem - mv)^2/float(length(data))
    end
    return s
end

function stand_dev(var_data)
    return sqrt(var_data)
end

function run_N_den(psi_init, t_q, time_range, ntraj)

    # psi = deepcopy(psi_init)
    N = length(psi_init)
    dt = time_range[2] - time_range[1]

    n_steps = length(time_range)

    # psi = deepcopy(psi_init)
    psi = convert_sp_dense(psi_init)

    den_err = zeros(length(time_range))

    for step=1:n_steps
        for traj=1:ntraj
            # println(traj)
            psi_dc = deepcopy(psi)
            zi_N, oi_N, pi_N = convert_bs_err_corr(mc_syn_gen_single_traj(psi_dc))
            # println(zi_N, oi_N, pi_N)
            # zi_N, oi_N, pi_N = convert_bs_err_corr_new(mc_syn_gen_single_traj_new(psi_dc, N))
            # den_err[step] = den_err[step] + length(zi_N)/N
            # den_err[step] = den_err[step] + length(oi_N)/N
            # den_err[step] = den_err[step] + length(pi_N)/N
            den_err[step] = den_err[step] + length(zi_N)/N + length(oi_N)/N + length(pi_N)/N
        end
        den_err[step] = den_err[step]/ntraj
        # println("begin evolution")
        evol_gates = evolution_tebd_gates_time_dependent_it(psi, dt, time_range[step], t_q)
        psi = apply(evol_gates, psi; cutoff=1E-8)
        normalize!(psi)
        # println("end evolution")
    end
    return den_err
end

function run_N_den_final_time(psi_init, t_q, time_range, ntraj)

    N = length(psi_init)
    dt = time_range[2] - time_range[1]

    # psi = deepcopy(psi_init)
    psi = convert_sp_dense(psi_init)

    n_steps = length(time_range)

    # println("begin evolution")
    for step=1:n_steps
        evol_gates = evolution_tebd_gates_time_dependent_it(psi, dt, time_range[step], t_q)
        psi = apply(evol_gates, psi; cutoff=1E-8)
        normalize!(psi)
    end
    # println("end evolution")

    den_err = 0.
    for traj=1:ntraj
        # println(traj)
        psi_dc = deepcopy(psi)
        zi_N, oi_N, pi_N = convert_bs_err_corr(mc_syn_gen_single_traj(psi_dc))
        # println(zi_N, oi_N, pi_N)
        den_err = den_err + length(zi_N)/N + length(oi_N)/N + length(pi_N)/N
    end
    return den_err/ntraj
end

function run_N_den_exp_minus(psi_init, t_q, time_range)#, ntraj)

    # psi = deepcopy(psi_init)
    N = length(psi_init)
    dt = time_range[2] - time_range[1]

    n_steps = length(time_range)

    psi = convert_sp_dense(psi_init)

    den_err = zeros(length(time_range))

    for step=1:n_steps
        den_err[step] = expect_minus(psi)

        # println("begin evolution")
        evol_gates = evolution_tebd_gates_time_dependent_it(psi, dt, time_range[step], t_q)
        # evol_gates = evolution_tebd_gates(psi, dt, time_range[step], t_q)
        psi = apply(evol_gates, psi; cutoff=1E-8)
        normalize!(psi)

        # println("end evolution")
    end
    return den_err
end

function run_N_den_exp_minus_final_time(psi_init, t_q, time_range)#, ntraj)

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

    return expect_minus(psi)
end

function run_N_den_site_exp_minus_final_time(psi_init, t_q, time_range)#, ntraj)

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

    return site_expect_minus(psi)
end

function run_N_den_exp_time_steps(psi_init, t_q, time_range)#, ntraj)

    N = length(psi_init)
    dt = time_range[2] - time_range[1]

    bd = []
    push!(bd, maxlinkdim(psi_init))

    psi = convert_sp_dense(psi_init)
    d_ts = []
    push!(d_ts, expect_minus(psi))

    n_steps = length(time_range)-1

    # println("begin evolution")
    for step=1:n_steps
        evol_gates = evolution_tebd_gates_time_dependent_it(psi, dt, time_range[step], t_q)
        psi = apply(evol_gates, psi; cutoff=1E-8)
        normalize!(psi)
        push!(bd, maxlinkdim(psi))
        push!(d_ts, expect_minus(psi))
    end
    # println("end evolution")

    return bd, d_ts
end
