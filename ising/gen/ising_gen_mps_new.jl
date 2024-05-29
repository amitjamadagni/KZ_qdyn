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


function proj_op(psi, i, j)


    sind_i = siteind(psi, i)
    sind_j = siteind(psi, j)
    # construct the operator
    # two_qubit_gate = ITensor(sind_j', sind_i', sind_j, sind_j)

    return ITensor([1. 0. 0. 0.; 0. 0. 0. 0.; 0 0. 0. 0.; 0. 0. 0. 1.], sind_j', sind_i', sind_j, sind_i)

end

function expectation_proj(psi, i, j, op_proj)
    orthogonalize!(psi, i)
    psi_ij = psi[i]*psi[j]
    psi_ij_p = prime(psi_ij, "Site")
    # println(scalar(psi_ij*op_proj*dag(psi_ij_p)))
    return real(scalar(psi_ij*op_proj*dag(psi_ij_p)))
    # return scalar(psi_ij*op_proj*dag(psi_ij_p))
end

function expectation_proj_new(psi, i, j, op_proj)
    orthogonalize!(psi, i)
    psi_n = apply(ITensor[op_proj], psi)
    # return inner(psi, dag(psi_n))
    # println(inner(psi, dag(psi_n)))
    return real(inner(psi, psi_n))
end

function expect_gs(psi)
    N = length(psi)

    ind_arr = [k for k=collect(1:1:N-1)]

    expect_m = []

    for (bs_ind, j) in enumerate(ind_arr)
        proj_op_minus = proj_op(psi, j, j+1)
        push!(expect_m, 1-expectation_proj(psi, j, j+1, proj_op_minus))
    end

    # push!(expect_m, 1-expectation_proj_new(psi, 1, N, proj_op(psi, 1, N)))

    return sum(expect_m)/float(N)
end

function site_expect_gs(psi)
    N = length(psi)

    ind_arr = [k for k=collect(1:1:N-1)]

    expect_m = []

    for (bs_ind, j) in enumerate(ind_arr)
        proj_op_minus = proj_op(psi, j, j+1)
        push!(expect_m, 1-expectation_proj(psi, j, j+1, proj_op_minus))
    end

    # push!(expect_m, 1-expectation_proj_new(psi, 1, N, proj_op(psi, 1, N)))

    # return sum(expect_m)/float(N)
    return expect_m
end
