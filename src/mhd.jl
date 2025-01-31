########################################################
####################  Functions  #######################
########################################################

# write the mesh to file for plotting
function write_mesh(mesh::Array{Float64, 2}, Bx, gamma, tl)
    filename = "data/mhd_TL$tl.out"

    open(filename, "w") do file
        for j in eachindex(mesh[1, :])
            write(file, cell_to_text(mesh[:, j], Bx, gamma))
            write(file, '\n')
        end
    end
end


# write a cell to text
function cell_to_text(cell::Array{Float64, 1}, Bx, gamma)
    rho, rho_vx, rho_vy, rho_vz, By, Bz, e = cell
    pressure = (gamma - 1)*(e - 0.5*sum([v^2 / rho for v in [rho_vx, rho_vy, rho_vz]]) -
                0.5*sum([b^2 for b in [Bx, By, Bz]]))
    return "$(rho),$(rho_vx / rho),$(rho_vy / rho),$(rho_vz / rho),$(By),$(Bz),$(e),$(pressure)"
end


# minmod function for approximating derivatives
function minmod(alpha, wjm, wj, wjp)
    if sign(wjm) == sign(wj) && sign(wj) == sign(wjp)
        delminus = alpha * (wj - wjm)
        delplus = alpha * (wjp - wj)
        delave = 0.5 * (delminus + delplus)

        return sign(wj) * min(abs(delminus), abs(delplus), abs(delave))
    else
        return 0
    end
end


# collect fluxes for all cells
function gather_fluxes!(fluxes::Array{Float64, 2}, mesh::Array{Float64, 2}, Bx, gamma)
    max_eig = 0
    for j in eachindex(mesh[1, :])
        eigj = gather_fluxes!(fluxes[:, j], mesh[:, j], Bx, gamma)
        if eigj > max_eig
            max_eig = eigj
        end
    end
    return max_eig
end


# calculate fluxes for one secc
function gather_fluxes!(flxs::Vector{Float64}, cell::Vector{Float64}, Bx, gamma)
    # variables
    rho, rho_vx, rho_vy, rho_vz, By, Bz, e = cell
    vx = rho_vx / rho
    vy = rho_vy / rho
    vz = rho_vz / rho

    Bmag2 = Bx^2 + By^2 + Bz^2
    vmag2 = vx^2 + vy^2 + vz^2

    # equation of state
    p = (gamma - 1)*(e - 0.5*rho*vmag2 - 0.5*Bmag2)

    pstar = p + 0.5*Bmag2

    # fluxes                                            # quantity
    flxs .= [
        rho_vx,                                         # rho
        rho_vx^2 / rho + pstar - Bx^2,                  # rho_vx
        rho_vx * vy - Bx * By,                          # rho_vy
        rho_vx * vz - Bx * Bz,                          # rho_vz
        By * vx - Bx * vy,                              # By
        Bz * vx - Bx * vz,                              # Bz
        (e + pstar)*vx - Bx*(Bx*vx + By*vy + Bz*vz)     # energy
       ]

    # We can get some of the CFL calculations for free doing it here
    return abs(vx) + eig(gamma * p / rho, Bmag2 / rho, Bz, rho)

end


eig(A, B, Bz, rho) = sqrt(0.5*(A + B + sqrt((A+B)^2 - 4*A*(Bz^2) / rho)))


########################################################
########################################################


function run(debug::Bool=false)
    # parameters
    nx = 800
    alpha = 1.4
    Bx = 0.75
    gamma = 2

    dx = 1/nx
    cfl = 0.4 * dx

    T_MAX = 0.2
    NUM_IMS = 30
    TOL = 0.0001
    numeq = 7

    # set initial and boundary conditions
    # first dimension is un/staggered
    mesh = Array{Float64}(undef, 2, numeq, nx)
    stag = 0
    # Brio-Wu shock tube
    e_init = 1 / (gamma - 1) - 0.5
    mesh[1, :, 1:(div(nx,2))] .= [1., 0., 0., 0., 1., 0., e_init]
    mesh[1, :, (div(nx,2)+1):nx] .= [0.125, 0., 0., 0., -1., 0., e_init / 10]

    mesh[2, :, :] .= copy(mesh[1, :, :])

    flxs = similar(mesh[1, :, :])

    # misc
    t = 0
    tl = 0
    update_partial = Vector{Float64}(undef, numeq)
    wjnhalf_arr = Vector{Float64}(undef, numeq)
    wjplusnhalf_arr = similar(wjnhalf_arr)
    fjnhalf = similar(wjnhalf_arr)
    fjplusnhalf = similar(wjnhalf_arr)

    while t <= T_MAX

        max_eig = gather_fluxes!(flxs, mesh[stag+1, :, :], Bx, gamma)

        dt = cfl / max_eig
        r = dt / dx
        rhalf = r / 2
        
        for j in 3:nx-2

            for eq = 1:numeq
                @views w = mesh[stag+1, eq, :]
                @views f = flxs[eq, :]
                wprime     = minmod(alpha, w[j-1], w[j],   w[j+1])
                wprimeplus = minmod(alpha, w[j-(2*stag)],   w[j-(2*stag)+1], w[j-(2*stag)+2])

                wjnhalf_arr[eq]     = w[j]   - rhalf * minmod(alpha, f[j-1], f[j],   f[j+1])
                wjplusnhalf_arr[eq] = w[j-(2*stag)+1] - rhalf * minmod(alpha, f[j-(2*stag)], f[j-(2*stag)+1], f[j-(2*stag)+2])

                update_partial[eq] = 0.5*(w[j] + w[j-(2*stag)+1]) + (-1*stag)*(wprime - wprimeplus) / 8
            end

            gather_fluxes!(fjnhalf,     wjnhalf_arr,     Bx, gamma)
            gather_fluxes!(fjplusnhalf, wjplusnhalf_arr, Bx, gamma)

            if stag == 0
                @. mesh[2, 1:numeq, j] = update_partial + r*(fjplusnhalf - fjnhalf)
            else
                @. mesh[1, 1:numeq, j] = update_partial - r*(fjplusnhalf - fjnhalf)
            end

        end

        # Track progress occasionally
        if tl % 1000 == 0
            print("t = $(round(t; digits=4))/$T_MAX\n")
        end

        # Check for nans (debug mode)
        if debug && sum(isnan.(mesh)) > 0
            print("NaN at time $t\n")
            return
        end

        # Output plot for animation
        if isapprox(t % (T_MAX / NUM_IMS), 0, atol=TOL)
            write_mesh(mesh[stag+1, :, :], Bx, gamma, tl)
        end

        # next time step we're back on the other grid
        if stag == 1
            stag = 0
        else
            stag = 1
        end
        t += dt
        tl += 1
    end
end

run(length(ARGS) > 0 && ARGS[1]=="debug")
