using Printf, HDF5, ITensors
using JLD: save
using SparseArrays

function run_xy_zpert(sites, N, J, Jp, dta)
    let
       	ampo = AutoMPO()

        for j=collect(1:2:N-1)
            ampo += (2*J, "Sx", j, "Sx", j+1)
            ampo += (-2*J, "iSy", j, "iSy", j+1)
            ampo += (2*J*dta, "Sz", j, "Sz", j+1)
        end
	for k=collect(2:2:N-2)
            ampo += (2*Jp, "Sx", k, "Sx", k+1)
            ampo += (-2*Jp, "iSy", k, "iSy", k+1)
            ampo += (2*Jp*dta, "Sz", k, "Sz", k+1)
        end

	H = MPO(ampo,sites)

        println("check pass ")

        # psi0 = randomMPS(sites)
        state = [isodd(n) ? "Up" : "Dn" for n in 1:N]
        psi0 = randomMPS(sites, state, N)
        # psi0 = productMPS(sites, state)

        sweeps = Sweeps(10)
        # Set maximum MPS bond dimensions for each sweep
        sw_arr = []
        maxdim!(sweeps, 10,20,100,200)
        # Set maximum truncation error allowed when adapting bond dimensions
        cutoff!(sweeps, 1E-15)
        @show sweeps

        # Run the DMRG algorithm, returning energy and optimized MPS
        energy, psi = dmrg(H,psi0, sweeps)
        @printf("Final energy = %.12f\n",energy)
        # for p in psi
        #     println(p)
        # end
	return energy, psi
    end
end

function convert_sp_dense(psi)
    psi_dc = deepcopy(psi)
    for i=1:length(psi_dc)
        psi_dc[i] = dense(psi_dc[i])
    end
    return psi_dc
end

function proj_op(psi, i, j)


    sind_i = siteind(psi, i)
    sind_j = siteind(psi, j)
    # construct the operator
    # two_qubit_gate = ITensor(sind_j', sind_i', sind_j, sind_j)

    return ITensor([0 0 0 0; 0 0.5 -0.5 0; 0 -0.5 0.5 0; 0 0 0 0], sind_j', sind_i', sind_j, sind_i)

end

function expectation_proj_new(psi, i, j, op_proj)
    orthogonalize!(psi, i)
    psi_n = apply(ITensor[op_proj], psi)
    # return inner(psi, dag(psi_n))
    # println(inner(psi, dag(psi_n)))
    return real(inner(psi, psi_n))
end

function expect_minus(psi)
    N = length(psi)

    ind_arr = [k for k=collect(2:2:N-2)]

    expect_m = []

    for (bs_ind, j) in enumerate(ind_arr)
        proj_op_minus = proj_op(psi, j, j+1)
        push!(expect_m, 1-expectation_proj_new(psi, j, j+1, proj_op_minus))
    end

    push!(expect_m, 1-expectation_proj_new(psi, 1, N, proj_op(psi, 1, N)))

    return sum(expect_m)/float(N)
end

#########

function gate_decomp(q1, q2, alpha, beta, gamma)
    theta = 0.5*pi - 2*gamma
    phi = 2*alpha - 0.5*pi
    lamb = 0.5*pi - 2*beta
    gates = [op("Rz", q2, ϕ=-0.5*pi),
         op("CX", q2, q1),
         op("Rz", q1, ϕ=theta),
         op("Ry", q2, θ=phi),
         op("CX", q1, q2),
         op("Ry", q2, θ=lamb),
         op("CX", q2, q1),
         op("Rz", q1, ϕ=0.5*pi),
         ]
    return gates
end

function evolution_tebd_gates_time_dependent_it(psi_init, dt, t, t_q, delta)

    N = length(psi_init)
    s = siteinds(psi_init)

    gates = ITensor[]

    for k=collect(1:1:N-1)
        s1 = s[k]
        s2 = s[k+1]
        if mod(k, 2) == 1
            Gj = gate_decomp(s1, s2, -0.25*dt*(1-(t/t_q)), -0.25*dt*(1-(t/t_q)), -0.25*dt*delta*(1-(t/t_q)))
            append!(gates, Gj)
        else
            Gj = gate_decomp(s1, s2, -0.25*dt*(1+(t/t_q)), -0.25*dt*(1+(t/t_q)), -0.25*dt*delta*(1+(t/t_q)))
            append!(gates, Gj)
        end
    end
    for k2=collect(N:-1:2)
        s1 = s[k2]
        s2 = s[k2-1]
        if mod(k2, 2) == 0
            Gj = gate_decomp(s1, s2, -0.25*dt*(1-(t/t_q)), -0.25*dt*(1-(t/t_q)), -0.25*dt*delta*(1-(t/t_q)))
            append!(gates, Gj)
        else
            Gj = gate_decomp(s1, s2, -0.25*dt*(1+(t/t_q)), -0.25*dt*(1+(t/t_q)), -0.25*dt*delta*(1+(t/t_q)))
            append!(gates, Gj)
        end
    end
    return gates
end
