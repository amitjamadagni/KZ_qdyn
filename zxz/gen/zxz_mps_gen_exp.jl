using Printf, HDF5, ITensors
using JLD: save
# include("../exact/graph_primer.jl")
# include("../exact/find_rescue_1d_periodic.jl")
using SparseArrays

function run_zxz_pert(sites, N, J, h1)
    let
        ampo = AutoMPO()

        for j=collect(1:1:N-2)
            ampo += (-8*J, "Sz", j, "Sx", j+1, "Sz", j+2)
        end

        for k=collect(1:1:N)
            ampo += (-2*h1, "Sx", k)
        end

        H = MPO(ampo, sites)

        println("check pass ")

        state = [isodd(n) ? "Up" : "Dn" for n in 1:N]
        psi0 = randomMPS(sites, state, N)

        sweeps = Sweeps(10)
        # Set maximum MPS bond dimensions for each sweep
        sw_arr = []
        maxdim!(sweeps, 10,20,100,200)
        # Set maximum truncation error allowed when adapting bond dimensions
        cutoff!(sweeps, 1E-15)
        @show sweeps

        # Run the DMRG algorithm, returning energy and optimized MPS
        energy, psi = dmrg(H, psi0, sweeps)
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

###########################

function proj_op(psi, i)
    sind_i = siteind(psi, i)
    sind_j = siteind(psi, i+1)
    sind_k = siteind(psi, i+2)

    # construct the operator
    # two_qubit_gate = ITensor(sind_j', sind_i', sind_j, sind_j)

    mat = [ 0.5  0.   0.5  0.   0.   0.   0.   0.;
            0.   0.5  0.  -0.5  0.   0.   0.   0.;
            0.5  0.   0.5  0.   0.   0.   0.   0.;
            0.  -0.5  0.   0.5  0.   0.   0.   0.;
            0.   0.   0.   0.   0.5  0.  -0.5  0.;
            0.   0.   0.   0.   0.   0.5  0.   0.5;
            0.   0.   0.   0.  -0.5  0.   0.5  0.;
            0.   0.   0.   0.   0.   0.5  0.   0.5]

    return ITensor(mat, sind_k', sind_j', sind_i', sind_k, sind_j, sind_i)
end


function expectation_proj(psi, i, op_proj)

    orthogonalize!(psi, i)

    psi_n = apply(ITensor[op_proj], psi)
    expect_arr = real(inner(psi, psi_n))

    return expect_arr
end

function expect_gs(psi)
    N = length(psi)

    ind_arr = [k for k=collect(1:1:N-2)]

    expect_m = []

    for (bs_ind, j) in enumerate(ind_arr)
        proj_op_minus = proj_op(psi, j)
        push!(expect_m, 1-expectation_proj(psi, j, proj_op_minus))
    end

    # push!(expect_m, 1-expectation_proj_new(psi, 1, N, proj_op(psi, 1, N)))

    return sum(expect_m)/float(N)
end

function site_expect_gs(psi)
    N = length(psi)

    ind_arr = [k for k=collect(1:1:N-2)]

    expect_m = []

    for (bs_ind, j) in enumerate(ind_arr)
        proj_op_minus = proj_op(psi, j)
        push!(expect_m, 1-expectation_proj(psi, j, proj_op_minus))
    end

    # push!(expect_m, 1-expectation_proj_new(psi, 1, N, proj_op(psi, 1, N)))

    # return sum(expect_m)/float(N)

    return expect_m
end

#####################################

function evolution_tebd_gates_time_dependent_it(psi_init, dt, t, t_q)

    N = length(psi_init)
    s = siteinds(psi_init)

    gates = ITensor[]

    hj1 = -8*op("Sz", s[1])*op("Sx", s[2])*op("Sz", s[3]) - 2*(-t/t_q)*op("Sx", s[1])*op("Id", s[2])*op("Id", s[3]) - (1/2)*2*(-t/t_q)*op("Id", s[1])*op("Sx", s[2])*op("Id", s[3]) - (1/3)*2*(-t/t_q)*op("Id", s[1])*op("Id", s[2])*op("Sx", s[3])
    Gj1 = exp(-0.5im*dt*hj1)
    push!(gates, Gj1)

    hj2 = -8*op("Sz", s[2])*op("Sx", s[3])*op("Sz", s[4]) - (1/2)*2*(-t/t_q)*op("Sx", s[2])*op("Id", s[3])*op("Id", s[4]) - (1/3)*2*(-t/t_q)*op("Id", s[2])*op("Sx", s[3])*op("Id", s[4]) - (1/3)*2*(-t/t_q)*op("Id", s[2])*op("Id", s[3])*op("Sx", s[4])
    Gj2 = exp(-0.5im*dt*hj2)
    push!(gates, Gj2)

    for k=collect(3:1:N-4)
        s1 = s[k]
        s2 = s[k+1]
        s3 = s[k+2]

        hj = -8*op("Sz", s1)*op("Sx", s2)*op("Sz", s3) - (1/3)*2*(-t/t_q)*op("Sx", s1)*op("Id", s2)*op("Id", s3) - (1/3)*2*(-t/t_q)*op("Id", s1)*op("Sx", s2)*op("Id", s3) - (1/3)*2*(-t/t_q)*op("Id", s1)*op("Id", s2)*op("Sx", s3)
        # Gj = exp(-0.5im*dt*hj)
        Gj = exp(-im*dt*hj)
        push!(gates, Gj)
    end

    hj3 = -8*op("Sz", s[end-3])*op("Sx", s[end-2])*op("Sz", s[end-1]) - (1/3)*2*(-t/t_q)*op("Sx", s[end-3])*op("Id", s[end-2])*op("Id", s[end-1]) - (1/3)*2*(-t/t_q)*op("Id", s[end-3])*op("Sx", s[end-2])*op("Id", s[end-1]) - (1/2)*2*(-t/t_q)*op("Id", s[end-3])*op("Id", s[end-2])*op("Sx", s[end-1])
    Gj3 = exp(-0.5im*dt*hj3)
    push!(gates, Gj3)

    hj4 = -8*op("Sz", s[end-2])*op("Sx", s[end-1])*op("Sz", s[end]) - (1/3)*2*(-t/t_q)*op("Sx", s[end-2])*op("Id", s[end-1])*op("Id", s[end]) - (1/2)*2*(-t/t_q)*op("Id", s[end-2])*op("Sx", s[end-1])*op("Id", s[end]) - 2*(-t/t_q)*op("Id", s[end-2])*op("Id", s[end-1])*op("Sx", s[end])
    Gj4 = exp(-0.5im*dt*hj4)
    push!(gates, Gj4)

    append!(gates, reverse(gates))
    return gates
end
