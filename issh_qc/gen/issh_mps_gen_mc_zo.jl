using Printf, HDF5, ITensors #, PastaQ
using JLD: save
# include("../exact/graph_primer.jl")
# include("../exact/find_rescue_1d_periodic.jl")
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
            Gj = gate_decomp(s2, s1, -0.25*dt*(1-(t/t_q)), -0.25*dt*(1-(t/t_q)), -0.25*dt*delta*(1-(t/t_q)))
            append!(gates, Gj)
        else
            Gj = gate_decomp(s2, s1, -0.25*dt*(1+(t/t_q)), -0.25*dt*(1+(t/t_q)), -0.25*dt*delta*(1+(t/t_q)))
            append!(gates, Gj)
        end
    end

    # append!(gates, reverse(gates))

    return gates
end

############# Ising like single site measurements

function mc_syn_gen_single_traj_ssh_ising(psi)
    N = length(psi)
    bs_arr = Array{Union{Nothing, Int64}}(nothing, N)
    for ind=1:N
        r = rand()
        proj_z = op(siteind(psi, ind), "projUp")
        proj_o = op(siteind(psi, ind), "projDn")
        orthogonalize!(psi, ind)
        pz = real(scalar(psi[ind]*proj_z*dag(prime(psi[ind], "Site"))))
        po = real(scalar(psi[ind]*proj_o*dag(prime(psi[ind], "Site"))))
        if r < pz
            bs_arr[ind] = 0
            # apply the projector
            pz_psi = proj_z*psi[ind]
            noprime!(pz_psi)
            psi[ind] = pz_psi
            orthogonalize!(psi, ind)
            psi[ind] ./= norm(psi[ind])
        elseif r < pz + po
            bs_arr[ind] = 1
            po_psi = proj_o*psi[ind]
            noprime!(po_psi)
            psi[ind] = po_psi
            orthogonalize!(psi, ind)
            psi[ind] ./= norm(psi[ind])
        end
    end
    return bs_arr
end

function convert_bs_err_corr_ssh_ising(bs)
    zi_n = 0
    oi_n = 0
    for i=2:2:length(bs)-1
        if bs[i] + bs[i+1] == 2
            oi_n = oi_n + 1
        elseif bs[i] + bs[i+1] == 0
            zi_n = zi_n + 1
        end
    end
    fl = bs[1] + bs[end]
    if fl == 2
        oi_n = oi_n + 1
    elseif fl == 0
        zi_n = zi_n + 1
    end
    return zi_n, oi_n
end

function run_N_den_ssh_ising_final_time(psi_init, ntraj)

    N = length(psi_init)
    # dt = time_range[2] - time_range[1]

    # psi = deepcopy(psi_init)
    psi = convert_sp_dense(psi_init)

    den_err = 0.

    for traj=1:ntraj
        # println(traj)
        psi_dc = deepcopy(psi)
        zi_N, oi_N = convert_bs_err_corr_ssh_ising(mc_syn_gen_single_traj_ssh_ising(psi_dc))#, N))
        den_err = den_err + zi_N/N
        den_err = den_err + oi_N/N
    end
    return den_err/ntraj
end
