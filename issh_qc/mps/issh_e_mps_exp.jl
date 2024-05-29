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

function run_N_den_ssh_ising_final_time(psi_init, t_q, time_range, delta)

    N = length(psi_init)
    dt = time_range[2] - time_range[1]

    # psi = deepcopy(psi_init)
    psi = convert_sp_dense(psi_init)

    n_steps = length(time_range)-1

    println("begin evolution")
    for step=1:n_steps
        evol_gates = evolution_tebd_gates_time_dependent_it(psi, dt, time_range[step], t_q, delta)
        psi = apply(evol_gates, psi; cutoff=1E-8)
        normalize!(psi)
    end
    println("end evolution")

    return expect_err_den(psi)

end
