using Printf, HDF5, ITensors
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

function fermion_no(psi)
    expect_arr = []
    for i=1:length(psi)
        psi_dc = deepcopy(psi)
        orthogonalize!(psi_dc, i)
        push!(expect_arr, scalar(psi_dc[i]*op(siteind(psi_dc, i), "projDn")*dag(prime(psi_dc[i], "n=$i"))))
    end
    return expect_arr
end

function proj_op(psi, i, j, op)
    sind_i = siteind(psi, i)
    sind_j = siteind(psi, j)
    # construct the operator
    # two_qubit_gate = ITensor(sind_j', sind_i', sind_j, sind_j)

    if op == "|00><00|"
        return ITensor([1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0], sind_j', sind_i', sind_j, sind_i)

    elseif op == "|+><+|"
        return ITensor([0 0 0 0; 0 0.5 0.5 0; 0 0.5 0.5 0; 0 0 0 0], sind_j', sind_i', sind_j, sind_i)

    elseif op == "|-><-|"
        return ITensor([0 0 0 0; 0 0.5 -0.5 0; 0 -0.5 0.5 0; 0 0 0 0], sind_j', sind_i', sind_j, sind_i)

    elseif op == "|11><11|"
        return ITensor([0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1], sind_j', sind_i', sind_j, sind_i)
    end

end

function expectation_proj_arr(psi, i, j, op_proj_arr)
    orthogonalize!(psi, i)
    psi_ij = psi[i]*psi[j]
    psi_ij_p = prime(psi_ij, "Site")
    expect_arr = []
    for op_proj in op_proj_arr
        push!(expect_arr, real(scalar(psi_ij*op_proj*dag(psi_ij_p))))
    end
    return expect_arr
end

function expectation_proj_arr_new(psi, i, j, op_proj_arr)
    orthogonalize!(psi, i)

    expect_arr = []

    for op_proj in op_proj_arr
        psi_n = apply(ITensor[op_proj], psi)
        # push!(expect_arr, scalar(psi_ij*op_proj*dag(psi_ij_p)))
        push!(expect_arr, inner(psi, dag(psi_n)))
    end
    return expect_arr
end

function mc_syn_gen_single_traj(psi)
    N = length(psi)
    bs_arr = Array{Union{Nothing, Int64}}(nothing, Int(N/2)+1)
    # bs_arr = sparsevec(zeros(Int(N/2)+2))
    r = rand()
    proj_z = op(siteind(psi, 1), "projUp")
    proj_o = op(siteind(psi, 1), "projDn")
    orthogonalize!(psi, 1)
    # compute expectation value
    # pz = real(scalar(psi[1]*proj_z*dag(prime(psi[1], "n=1"))))
    # po = real(scalar(psi[1]*proj_o*dag(prime(psi[1], "n=1"))))
    pz = real(scalar(psi[1]*proj_z*dag(prime(psi[1], "Site"))))
    po = real(scalar(psi[1]*proj_o*dag(prime(psi[1], "Site"))))
    if r < pz
        bs_arr[1] = 0
        # apply the projector
        pz_psi = proj_z*psi[1]
        noprime!(pz_psi)
        psi[1] = pz_psi
        orthogonalize!(psi, 1)
        psi[1] ./= norm(psi[1])
    elseif r < pz + po
        bs_arr[1] = 1
        po_psi = proj_o*psi[1]
        noprime!(po_psi)
        psi[1] = po_psi
        orthogonalize!(psi, 1)
        psi[1] ./= norm(psi[1])
    end
    ind_arr = [k for k=collect(2:2:N-2)]
    for (bs_ind,j) in enumerate(ind_arr)
        r = rand()
        proj_op_arr = [proj_op(psi, j, j+1, "|00><00|"), proj_op(psi, j, j+1, "|+><+|"), proj_op(psi, j, j+1, "|-><-|"), proj_op(psi, j, j+1, "|11><11|")]
        p_zero, p_plus, p_minus, p_one = expectation_proj_arr(psi, j, j+1, proj_op_arr)
        # println("###")
        # println(p_zero, " ", p_plus, " ", p_minus, " ", p_one)
        # println("###")
        if r < p_zero
            # println(" r : ", r, " p0 : ", p_zero)
            bs_arr[bs_ind+1] = 0
            op_psi_j_jp1 = proj_op_arr[1]*psi[j]*psi[j+1]
            noprime!(op_psi_j_jp1)
            # U,S,V = svd(op_psi_j_jp1, siteind(psi, j))

            inds3 = uniqueinds(psi[j],psi[j+1])
            U,S,V = svd(op_psi_j_jp1, inds3, cutoff=1E-8)
            psi[j] = U
            psi[j+1] = S*V
            # renormalize
            orthogonalize!(psi, j)
            psi[j] ./= norm(psi[j])

        elseif r < p_plus + p_zero
            # println(" r : ", r, " p+0 : ", p_plus + p_zero)
            bs_arr[bs_ind+1] = 2
            op_psi_j_jp1 = proj_op_arr[2]*psi[j]*psi[j+1]
            noprime!(op_psi_j_jp1)
            # U,S,V = svd(op_psi_j_jp1, siteind(psi, j))

            inds3 = uniqueinds(psi[j],psi[j+1])
            U,S,V = svd(op_psi_j_jp1, inds3, cutoff=1E-8)
            psi[j] = U
            psi[j+1] = S*V
            # renormalize
            orthogonalize!(psi, j)
            psi[j] ./= norm(psi[j])

        elseif r < p_minus + p_plus + p_zero
            # println(" r : ", r, " p-+0 : ", p_minus + p_plus + p_zero)
            # println(" diff : ", r-(p_minus + p_plus + p_zero))
            bs_arr[bs_ind+1] = -1
            op_psi_j_jp1 = proj_op_arr[3]*psi[j]*psi[j+1]
            noprime!(op_psi_j_jp1)
            # U,S,V = svd(op_psi_j_jp1, siteind(psi, j))

            inds3 = uniqueinds(psi[j],psi[j+1])
            U,S,V = svd(op_psi_j_jp1, inds3, cutoff=1E-8)

            psi[j] = U
            psi[j+1] = S*V
            # renormalize
            orthogonalize!(psi, j)
            psi[j] ./= norm(psi[j])

        elseif r < p_one + p_minus + p_plus + p_zero
            # println(" r : ", r, " p1-+0 : ", p_one + p_minus + p_plus + p_zero)
            # println(" diff : ", r-(p_minus + p_plus + p_zero))
            bs_arr[bs_ind+1] = 1
            op_psi_j_jp1 = proj_op_arr[4]*psi[j]*psi[j+1]
            noprime!(op_psi_j_jp1)
            # U,S,V = svd(op_psi_j_jp1, siteind(psi, j))

            inds3 = uniqueinds(psi[j],psi[j+1])
            U,S,V = svd(op_psi_j_jp1, inds3, cutoff=1E-8)

            psi[j] = U
            psi[j+1] = S*V
            # renormalize
            orthogonalize!(psi, j)
            psi[j] ./= norm(psi[j])
        end
    end
    proj_z = op(siteind(psi, N), "projUp")
    proj_o = op(siteind(psi, N), "projDn")
    orthogonalize!(psi, N)
    # pz = real(scalar(psi[N]*proj_z*dag(prime(psi[N], "n=$N"))))
    # po = real(scalar(psi[N]*proj_o*dag(prime(psi[N], "n=$N"))))
    pz = real(scalar(psi[N]*proj_z*dag(prime(psi[N], "Site"))))
    po = real(scalar(psi[N]*proj_o*dag(prime(psi[N], "Site"))))
    # println(pz+po)
    r = rand()
    if r < pz
        bs_arr[end] = 0
    elseif r < pz + po
        bs_arr[end] = 1
    end
    return bs_arr
end

function convert_bs_err_corr(bs)
    bs_reconst = []
    if bs[1] == 0 && bs[end] == 0
        push!(bs_reconst, 0)
    elseif bs[1] == 1 && bs[end] == 1
        push!(bs_reconst, 1)
    elseif bs[1] == 0 && bs[end] == 1
        push!(bs_reconst, 2)
    else
        push!(bs_reconst, -1)
    end
    for j=2:length(bs)-1
        push!(bs_reconst, bs[j])
    end
    zi = []
    oi = []
    pi_N = []
    for (i, elem) in enumerate(bs_reconst)
        if elem == 0
            push!(zi, i)
        elseif elem == 1
            push!(oi, i)
        elseif elem == 2 && i!=1
            push!(pi_N, i)
        end
    end
    return zi, oi, pi_N
end

function mc_syn_gen_single_traj_new(psi)#, N)
    N = length(psi)
    # println(psi)
    bs_arr = Array{Union{Nothing, Int64}}(nothing, Int(N/2))
    ind_arr = [k for k=collect(2:2:N-2)]
    for (bs_ind,j) in enumerate(ind_arr)
        r = rand()
        proj_op_arr = [proj_op(psi, j, j+1, "|00><00|"), proj_op(psi, j, j+1, "|+><+|"), proj_op(psi, j, j+1, "|-><-|"), proj_op(psi, j, j+1, "|11><11|")]
        p_zero, p_plus, p_minus, p_one = expectation_proj_arr(psi, j, j+1, proj_op_arr)

        if r < p_zero
            bs_arr[bs_ind+1] = 0
            op_psi_j_jp1 = proj_op_arr[1]*psi[j]*psi[j+1]
            noprime!(op_psi_j_jp1)
            jnds3 = uniqueinds(psi[j],psi[j+1])
            # println(jnds3)
            # U,S,V = svd(op_psi_j_jp1, siteind(psi, j))
            U,S,V = svd(op_psi_j_jp1, jnds3)
            psi[j] = U
            psi[j+1] = S*V
            # renormalize
            orthogonalize!(psi, j)
            psi[j] ./= norm(psi[j])

        elseif r < p_plus + p_zero
            bs_arr[bs_ind+1] = 2
            op_psi_j_jp1 = proj_op_arr[2]*psi[j]*psi[j+1]
            noprime!(op_psi_j_jp1)
            jnds3 = uniqueinds(psi[j],psi[j+1])
            # U,S,V = svd(op_psi_j_jp1, siteind(psi, j))
            U,S,V = svd(op_psi_j_jp1, jnds3)
            psi[j] = U
            psi[j+1] = S*V
            # renormalize
            orthogonalize!(psi, j)
            psi[j] ./= norm(psi[j])

        elseif r < p_minus + p_plus + p_zero
            bs_arr[bs_ind+1] = -1
            op_psi_j_jp1 = proj_op_arr[3]*psi[j]*psi[j+1]
            noprime!(op_psi_j_jp1)
            jnds3 = uniqueinds(psi[j],psi[j+1])
            # U,S,V = svd(op_psi_j_jp1, siteind(psi, j))
            U,S,V = svd(op_psi_j_jp1, jnds3)
            psi[j] = U
            psi[j+1] = S*V
            # renormalize
            orthogonalize!(psi, j)
            psi[j] ./= norm(psi[j])

        elseif r < p_one + p_minus + p_plus + p_zero
            bs_arr[bs_ind+1] = 1
            op_psi_j_jp1 = proj_op_arr[4]*psi[j]*psi[j+1]
            noprime!(op_psi_j_jp1)
            jnds3 = uniqueinds(psi[j],psi[j+1])
            # U,S,V = svd(op_psi_j_jp1, siteind(psi, j))
            U,S,V = svd(op_psi_j_jp1, jnds3)
            psi[j] = U
            psi[j+1] = S*V
            # renormalize
            orthogonalize!(psi, j)
            psi[j] ./= norm(psi[j])
        end
    end

    println(norm(psi))
    proj_op_arr = [proj_op(psi, 1, N, "|00><00|"), proj_op(psi, 1, N, "|+><+|"), proj_op(psi, 1, N, "|-><-|"), proj_op(psi, 1, N, "|11><11|")]
    # p_zero, p_plus, p_minus, p_one = expectation_proj_arr(psi, 1, N, proj_op_arr)
    p_zero, p_plus, p_minus, p_one = expectation_proj_arr_new(psi, 1, N, proj_op_arr)

    println(p_zero + p_plus + p_minus + p_one)

    r = rand()
    if r < real(p_zero)
        bs_arr[1] = 0
    elseif r < real(p_zero) + real(p_plus)
        bs_arr[1] = 2
    elseif r < real(p_zero) + real(p_plus) + real(p_minus)
        bs_arr[1] = -1
    else
        bs_arr[1] = 1
    end

    return bs_arr
end

function convert_bs_err_corr_new(bs)
    zi = []
    oi = []
    pi_N = []
    for (i, elem) in enumerate(bs)
        if elem == 0
            push!(zi, i)
        elseif elem == 1
            push!(oi, i)
        elseif elem == 2
            push!(pi_N, i)
        end
    end
    return zi, oi, pi_N
end

function run_N_mc_final_time(psi_init, ntraj)

    N = length(psi_init)

    # psi = deepcopy(psi_init)
    psi = convert_sp_dense(psi_init)

    dz_N = 0.
    do_N = 0.
    dp_N = 0.

    den_err = 0.
    for traj=1:ntraj
        # println(traj)
        psi_dc = deepcopy(psi)
        zi_N, oi_N, pi_N = convert_bs_err_corr(mc_syn_gen_single_traj(psi_dc))
        dz_N = dz_N + length(zi_N)/N
        do_N = do_N + length(oi_N)/N
        dp_N = dp_N + length(pi_N)/N
        # println(zi_N, oi_N, pi_N)
        den_err = den_err + length(zi_N)/N + length(oi_N)/N + length(pi_N)/N
    end
    println(dz_N/ntraj)
    println(do_N/ntraj)
    println(dp_N/ntraj)
    return den_err/ntraj
end
