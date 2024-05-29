# https://tensornetwork.org/mps/algorithms/timeevo/tebd.html
function evolution_tebd_gates_time_dependent(psi_init, dt, t, t_q)

    N = length(psi_init)
    s = siteinds(psi_init)

    gates = ITensor[]

    gates_e = ITensor[]
    for k=collect(2:2:N-2)
        s1 = s[k]
        s2 = s[k+1]

        hj = -4*op("Sz", s1)*op("Sz", s2) + (t/t_q)*op("Sx", s1)*op("Id", s2) + (t/t_q)*op("Id", s1)*op("Sx", s2)

        Gj = exp(-0.5im*dt*hj)
        push!(gates_e, Gj)
    end

    gates_o = ITensor[]

    hj1 = -4*op("Sz", s[1])*op("Sz", s[2]) + 2*(t/t_q)*op("Sx", s[1])*op("Id", s[2]) + (t/t_q)*op("Id", s[1])*op("Sx", s[2])

    Gj1 = exp(-1.0im*dt*hj1)
    push!(gates_o, Gj1)

    for j=collect(3:2:N-3)
        s1 = s[j]
        s2 = s[j+1]

        hj = -4*op("Sz", s1)*op("Sz", s2) + (t/t_q)*op("Sx", s1)*op("Id", s2) + (t/t_q)*op("Id", s1)*op("Sx", s2)

        Gj = exp(-1.0im*dt*hj)
        push!(gates_o, Gj)
    end

    hj2 = -4*op("Sz", s[end-1])*op("Sz", s[end]) + (t/t_q)*op("Sx", s[end-1])*op("Id", s[end]) + 2*(t/t_q)*op("Id", s[end-1])*op("Sx", s[end])

    Gj2 = exp(-1.0im*dt*hj2)
    push!(gates_o, Gj2)

    append!(gates, gates_e)
    append!(gates, gates_o)
    append!(gates, gates_e)
    return gates
end

function evolution_tebd_gates_time_dependent_it(psi_init, dt, t, t_q)

    N = length(psi_init)
    s = siteinds(psi_init)

    gates = ITensor[]

    hj1 = -4*op("Sz", s[1])*op("Sz", s[2]) + 2*(t/t_q)*op("Sx", s[1])*op("Id", s[2]) + (t/t_q)*op("Id", s[1])*op("Sx", s[2])
    Gj1 = exp(-0.5im*dt*hj1)
    push!(gates, Gj1)

    for k=collect(2:1:N-2)
        s1 = s[k]
        s2 = s[k+1]

        hj = -4*op("Sz", s1)*op("Sz", s2) + (t/t_q)*op("Sx", s1)*op("Id", s2) + (t/t_q)*op("Id", s1)*op("Sx", s2)

        Gj = exp(-0.5im*dt*hj)
        push!(gates, Gj)
    end

    hj2 = -4*op("Sz", s[end-1])*op("Sz", s[end]) + (t/t_q)*op("Sx", s[end-1])*op("Id", s[end]) + 2*(t/t_q)*op("Id", s[end-1])*op("Sx", s[end])
    Gj2 = exp(-0.5im*dt*hj2)
    push!(gates, Gj2)

    append!(gates, reverse(gates))
    return gates
end


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

    N = length(psi_init)
    dt = time_range[2] - time_range[1]

    n_steps = length(time_range)

    err_time_t_arr = zeros(n_steps, ntraj)

    # psi = deepcopy(psi_init)
    psi = convert_sp_dense(psi_init)

    den_err = zeros(length(time_range))

    for step=1:n_steps
        for traj=1:ntraj
            # println(traj)
            psi_dc = deepcopy(psi)
            bs_straj = mc_syn_gen_single_traj(psi_dc)
            oi_N = convert_bs_err_corr(bs_straj)
            # println(zi_N, oi_N, pi_N)
            den_err[step] = den_err[step] + length(oi_N)/N
        end
        den_err[step] = den_err[step]/ntraj
        println("begin evolution")
        evol_gates = evolution_tebd_gates_time_dependent_it(psi, dt, time_range[step], t_q)
        psi = apply(evol_gates, psi; cutoff=1E-8)
        # psi = apply(evol_gates_e, psi; cutoff=1E-8)
        # psi = apply(evol_gates_o, psi; cutoff=1E-8)
        println("end evolution")
    end
    return den_err
end

function run_N_den_final_time(psi_init, t_q, time_range, ntraj)

    N = length(psi_init)
    dt = time_range[2] - time_range[1]

    # psi = deepcopy(psi_init)
    psi = convert_sp_dense(psi_init)

    n_steps = length(time_range)-1

    println("begin evolution")
    for step=1:n_steps
        evol_gates = evolution_tebd_gates_time_dependent_it(psi, dt, time_range[step], t_q)
        psi = apply(evol_gates, psi; cutoff=1E-8)
        normalize!(psi)
    end
    println("end evolution")

    den_err = 0.
    for traj=1:ntraj
        # println(traj)
        psi_dc = deepcopy(psi)
        bs_straj = mc_syn_gen_single_traj(psi_dc)
        oi_N = convert_bs_err_corr(bs_straj)
        # println(zi_N, oi_N, pi_N)
        den_err = den_err + length(oi_N)/N
    end
    return den_err/ntraj
end

function run_N_den_exp_final_time(psi_init, t_q, time_range)#, ntraj)

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

function run_N_den_site_exp_final_time(psi_init, t_q, time_range)#, ntraj)

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
    push!(d_ts, expect_gs(psi_init))

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
