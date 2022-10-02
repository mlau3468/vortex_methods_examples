using WriteVTK
using Revise
using PAGE
include("../include_all.jl")

airfoil_data = ["naca0012.csv";]

# initialize wing
nchord = 6
wing = WingNLVLM(nsec=2, nchord=nchord, refloc=0.0, airfoil_data=airfoil_data)

camber = zeros(2, nchord+1)
camber[1,:] .= LinRange(0, 1, nchord+1)
dummy = zeros(nchord+1)
@views PAGE.meancamber2!(camber[2,:], dummy, 0.0, 0.25, camber[1,:])

span = 2.708707
chord = 0.69776
define_wing_section!(wing, 1, camber, [0.0;-span/2;0.0], chord, 5.0, 0.0)
define_wing_section!(wing, 2, camber, [0.0;span/2;0.0], chord, 5.0, 0.0)
set_wing_region!(wing, 1, 1, 30; spacetype="equal")

# generate mesh
pansys1 = generate_mesh(wing, true)

uinf = [30;0;0]
wakedist = chord*30

pansys = create_pansys(pansys1)

wakesys = create_te_vortlat!(pansys, uinf, wakedist, true)

vis_mesh(pansys.pan_con, pansys.pan_vert, "pansys")
vis_mesh(wakesys.te_con, wakesys.te_vert, "wake")

function solver(pansys, wakesys, uinf)
    r_rankine = 0.01
    r_cutoff = 0.001
    npan = pansys.npan
    pan_cpt = pansys.pan_cpt
    pan_norm = pansys.pan_norm
    pan_vert = pansys.pan_vert
    pan_con = pansys.pan_con
    pan_edgelen = pansys.pan_edgelen
    pan_edgeuvec = pansys.pan_edgeuvec
    pan_edgevec = pansys.pan_edgevec
    pan_numpt = pansys.pan_numpt
    pan_neigh_idx = pansys.pan_neigh_idx
    pan_area = pansys.pan_area
    airfoil_idx = pansys.airfoil_idx
    airfoil_cl = pansys.airfoil_cl
    airfoil_cd = pansys.airfoil_cd
    nte = wakesys.nte
    te_pan_idx = wakesys.te_pan_idx
    te_con = wakesys.te_con
    te_vert = wakesys.te_vert
    te_numpt = wakesys.te_numpt
    te_edgelen = wakesys.te_edgelen
    te_edgeuvec = wakesys.te_edgeuvec
    airfoil_strips = pansys.airfoil_strips
    nchord = size(airfoil_strips, 1)
    nstrips = size(airfoil_strips, 2)
    rhoinf = 1.225

    max_iter = 50

    # Compute strip properties
    alpha_strip_cpt = zeros(3, nstrips) # point to compute local strip alpha. 1/4 chord
    chord_strip = zeros(nstrips)
    norm_strip = zeros(3, nstrips)
    tan_strip = zeros(3, nstrips)
    area_strip = zeros(nstrips)
    for i = 1:nstrips
        p1 = pan_vert[:,pan_con[1,airfoil_strips[1,i]]]
        p2 = pan_vert[:,pan_con[2,airfoil_strips[nchord,i]]]
        p3 = pan_vert[:,pan_con[3,airfoil_strips[nchord,i]]]
        p4 = pan_vert[:,pan_con[4,airfoil_strips[1,i]]]
        
        #new_pt = 0.75.*(p1.+p4)./2 .+ 0.25.*(p2.+p3)./2
        new_pt = (p1 .+ p2 .+ p3 .+ p4)./4
        alpha_strip_cpt[:,i] .= new_pt
        chord_strip[i] = (norm(p2.-p1) + norm(p4.-p3))/2

        v3 = cross(p3.-p1, p4.-p2)
        norm_strip[:,i] .= v3./norm(v3)

        v1 = p2.-p1
        v2 = p3.-p4
        v1 = v1./norm(v1)
        v2 = v2./norm(v2)
        tan_strip[:,i] .= (v1 .+ v2)./2
        for j = 1:nchord
            idx = airfoil_strips[j,i]
            area_strip[i] += pan_area[idx]
        end
    end

    # Linear VLM solver
    A = zeros(npan, npan)
    RHS = zeros(npan)
    for i = 1:npan
        for j = 1:npan
            vel = zeros(3)
            #vel_dubpan!(vel, pan_cpt[:,i], 1.0, pan_vert, pan_con[:,j], pan_numpt[j], pan_edgelen[:,j], pan_edgeuvec[:,:,j], r_rankine, r_cutoff)
            vel_vring!(vel, pan_cpt[:,i], 1.0, pan_vert, pan_con[:,j], pan_numpt[j], pan_edgelen[:,j], pan_edgeuvec[:,:,j], r_rankine, r_cutoff)
            A[i,j] += dot(vel, pan_norm[:,i])
        end

        for j = 1:nte
            idx = te_pan_idx[j]
            vel = zeros(3)
            vel_vring!(vel, pan_cpt[:,i], 1.0, te_vert, te_con[:,j], te_numpt[j], te_edgelen[:,j], te_edgeuvec[:,:,j], r_rankine, r_cutoff)
            A[i,idx] += dot(vel, pan_norm[:,i])
        end
    end

    for i = 1:npan
        RHS[i] = -dot(uinf, pan_norm[:,i])
    end

    iter = 0
    converged = false
    tol = 1e-4
    rlx = 0.7
    while iter < max_iter && !converged
        iter += 1
        gam_sol = A\RHS

        # velocities at vlm bound vortex segments
        vels = zeros(3, npan)
        for i = 1:npan
            coord = (pan_vert[:,pan_con[1,i]] .+ pan_vert[:,pan_con[4,i]])./2
            for j = 1:npan
                if i == j
                    @views vel_vring2!(vels[:,i], coord, gam_sol[j], pan_vert, pan_con[:,j], pan_numpt[j], pan_edgelen[:,j], pan_edgeuvec[:,:,j], r_rankine, r_cutoff)
                elseif j == pan_neigh_idx[4,i] # influencing panel is directly upstream
                    @views vel_vring3!(vels[:,i], coord, gam_sol[j], pan_vert, pan_con[:,j], pan_numpt[j], pan_edgelen[:,j], pan_edgeuvec[:,:,j], r_rankine, r_cutoff)
                else
                    @views vel_vring!(vels[:,i], coord, gam_sol[j], pan_vert, pan_con[:,j], pan_numpt[j], pan_edgelen[:,j], pan_edgeuvec[:,:,j], r_rankine, r_cutoff)
                end
            end
            for j = 1:nte
                @views vel_vring!(vels[:,i], coord, gam_sol[te_pan_idx[j]], te_vert, te_con[:,j], te_numpt[j], te_edgelen[:,j], te_edgeuvec[:,:,j], r_rankine, r_cutoff)
            end

            vels[:,i] .+= uinf
        end

        # vlm forces
        forces = zeros(3, npan)
        vmags = zeros(npan)
        for i = 1:npan
            vmags[i] = norm(vels[:,i])
            if pan_neigh_idx[4,i] > 0
                gam = gam_sol[i] - gam_sol[pan_neigh_idx[4,i]]
            else
                gam = gam_sol[i]
            end
            forces[:,i] .= rhoinf.*gam.*cross(vels[:,i],pan_edgevec[:,4,i])
        end

        # strip invscid cl
        cl_invisc_strip = zeros(nstrips)
        mag_inv_strip = zeros(nstrips)
        vmag_strip = zeros(nstrips)
        for i = 1:nstrips
            lift = 0.0
            for j = 1:nchord
                idx = airfoil_strips[j,i]
                lift += dot(forces[:,idx], [0;0;1])
                vmag_strip[i] += vmags[idx]
            end
            vmag_strip[i] = vmag_strip[i]/nchord
            cl_invisc_strip[i] = lift/(0.5*rhoinf*vmag_strip[i]^2*area_strip[i])

            for j = 1:nchord
                if j == 1
                    mag_inv_strip[i] += gam_sol[airfoil_strips[j,i]]
                else
                    mag_inv_strip[i] += gam_sol[airfoil_strips[j,i]] - gam_sol[airfoil_strips[j-1,i]]
                end
            end
            cl_invisc_strip[i] = -2*mag_inv_strip[i]/vmag_strip[i]/chord_strip[i]
        end

        # Compute strip aoa
        vel_strip = zeros(3, nstrips)
        aoa_strip = zeros(nstrips)
        cl_visc_strip = zeros(nstrips)
        cl_diff = zeros(nstrips)
        residual = zeros(npan)
        for i = 1:nstrips
            @views compute_vel4!(vel_strip[:,i], alpha_strip_cpt[:,i], pansys, wakesys, gam_sol, uinf, r_rankine, r_cutoff)
            uperp = dot(vel_strip[:,i], norm_strip[:,i])
            utan = dot(vel_strip[:,i], tan_strip[:,i])
            aoa_strip[i] = atan(uperp, utan)

            # 2d correction to induced angle
            alpha_2d = mag_inv_strip[i] / (pi*chord_strip[i]*vmag_strip[i])
            aoa_strip[i] -= alpha_2d

            cl_visc_strip[i] = airfoil_cl[airfoil_idx[i]](rad2deg(aoa_strip[i]))

            cl_diff[i] = cl_visc_strip[i] - cl_invisc_strip[i]

            for j = 1:nchord
                idx = airfoil_strips[j,i]
                residual[idx] = cl_diff[i]*RHS[idx]
            end
        end

        if maximum(abs.(cl_diff)) < tol
            converged = true
        end

        for i = 1:npan
            RHS[i] = RHS[i] + rlx*residual[i]
        end

        f_total = sum(forces, dims=2)
        println(f_total)
    end

    if !converged
        println("not converged")
    end

    #vis_mesh_temp(pansys.pan_con, pansys.pan_vert, gam_sol, vmags, "sol")
end

solver(pansys, wakesys, uinf)