using Printf, HDF5, ITensors
using JLD: save
using SparseArrays

function run_ising_xpert(J, sites, N, hx)
    let
        ampo = AutoMPO()

        for j=collect(1:N-1)
            ampo += (-4*J, "Sz", j, "Sz", j+1)
        end
        for k=collect(1:N)
            ampo += (-2*hx, "Sx", k)
        end

        H = MPO(ampo,sites)

        psi0 = randomMPS(ComplexF64, sites)
        # state = [isodd(n) ? "Up" : "Dn" for n in 1:N]
        # psi0 = randomMPS(ComplexF64, sites, state)
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


function mc_syn_gen_single_traj(psi)
    N = length(psi)
    bs_arr = Array{Union{Nothing, Int64}}(nothing, N)
    for ind=1:N
        r = rand()
        proj_z = op(siteind(psi, ind), "projUp")
        proj_o = op(siteind(psi, ind), "projDn")
        orthogonalize!(psi, ind)
        # compute expectation value
        # pz = real(scalar(psi[ind]*proj_z*dag(prime(psi[ind], "n=$ind"))))
        # po = real(scalar(psi[ind]*proj_o*dag(prime(psi[ind], "n=$ind"))))
        pz = real(scalar(psi[ind]*proj_z*dag(prime(psi[ind], "Site"))))
        po = real(scalar(psi[ind]*proj_o*dag(prime(psi[ind], "Site"))))
        # println(" 0  : ", pz)
        # println(" 1  : ", po)
        # println(" sum : ", pz + po)
        if r < pz
            # println("bs : 0")
            bs_arr[ind] = 0
            # apply the projector
            pz_psi = proj_z*psi[ind]
            noprime!(pz_psi)
            psi[ind] = pz_psi
            orthogonalize!(psi, ind)
            psi[ind] ./= norm(psi[ind])
        # elseif r < pz + po
        else
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

function convert_bs_err_corr(bs)
    bs_reconst = []
    for i=1:length(bs)-1
        if bs[i] + bs[i+1] == 1
            push!(bs_reconst, 1)
        elseif bs[i] + bs[i+1] == 2
            push!(bs_reconst, 0)
        elseif bs[i] + bs[i+1] == 0
            push!(bs_reconst, 0)
        end
    end
    oi = findall(isequal(1), bs_reconst)
    return oi
end

function mc_time_stats_single_traj(psi, tc_non_opt, tc_mat, traj)
    psi_dc = deepcopy(psi)
    bs_straj = mc_syn_gen_single_traj(psi_dc)
    oi_N = convert_bs_err_corr(bs_straj)
    str_time = max_time_to_gs_1d_fr(oi_N .+ 1, tc_non_opt, tc_mat)
    return str_time
end


function mc_syn_gen_ntraj(psi, ntraj)
    bs_arr_traj = []
    for traj=1:ntraj
        psi_dc = deepcopy(psi)
        push!(bs_arr_traj, mc_syn_gen_single_traj(psi_dc))
    end
    return bs_arr_traj
end
