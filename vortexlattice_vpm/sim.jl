include("aeroel.jl")
include("geo.jl")
include("vis.jl")


function simulate(components, panels, te_idx, tsteps, dt, uinf, rho, particles, wakelines, wakerings, prefix)

    # test geometry motion
    test_geo(components, panels, dt, prefix)


    # Initialize
    setPointVels!(components, panels)
    panels_neigh, panels_neighside, panels_neighdir = calcneighbors(panels)
    tewidth = length(te_idx)
    maxwakelen = 1
    wakelen = 0

    # trailing edge variables
    te_neigh = Int64[]
    te_neighdir = Int64[]
    te_neighside = Int64[]

    A = zeros(Float64,length(panels), length(panels))
    RHS = zeros(Float64,length(panels))

    panels2vtk(panels, prefix * "_panels_0.vtu")
    particles2vtk(particles, prefix * "_particles_0.vtu")
    wakepanels2vtk(wakerings, prefix * "_wakerings_0.vtu")

    # timestep
    for t = 1:tsteps
        # calculate influence coefficients
        for i = 1:length(panels)
            for j = 1:length(panels)
                # influence of jth panel on ith collocation point
                vel = vrtxring(panels[j].pts, panels[i].cpt, 1)
                A[i,j] = dot(vel, panels[i].normal)
            end
        end

        # build rhs vector
        RHS[:] .= 0.0
        for i = 1:length(panels)
            RHS[i] = -dot(uinf, panels[i].normal)
            RHS[i] = RHS[i] - dot(panels[i].wake_vel, panels[i].normal)
            RHS[i] = RHS[i] + dot(panels[i].vcpt, panels[i].normal)
        end

        # solve matrix for panel gamma
        sol = A\RHS
        for i = 1:length(panels)
            newPanGam(panels[i], sol[i], dt)
        end

        # calculate new wake elements
        new_wakerings = []
        for i = 1:length(te_idx)
            idx = te_idx[i]
            p1 = panels[idx].pts[:,4]
            p2 = panels[idx].pts[:,3]

            # calculate where the points will be next
            omega = components[panels[idx].compidx][1].omega
            vbody = components[panels[idx].compidx][1].vel
            rotm = stepRotMat(omega, dt)
            p1_new = rotm*p1 .+ vbody.*dt
            p2_new = rotm*p2 .+ vbody.*dt

            vel1 = elemVel(panels, particles, wakelines, wakerings, p1) .+ uinf 
            vel2 = elemVel(panels, particles, wakelines, wakerings,p2) .+ uinf 
            gam = panels[idx].last_gam[1]
            p3 = p2 .+ vel2.*dt
            p4 = p1 .+ vel1.*dt
            new_wakering = initWakeRing([p1_new p2_new p3 p4])
            new_wakering.gam[1] = gam
            push!(new_wakerings, new_wakering)
        end

        # move wake
        stepWake!(panels, particles, wakelines, wakerings, uinf, dt)
        
        # move geometry
        stepGeometry!(components, panels, dt)
        
        # add new wakerings
        wakerings = cat(new_wakerings, wakerings, dims=1)

        wakelen = wakelen + 1

        # update panel wake_vel
        for i =1:length(panels)
            panels[i].wake_vel[:] = wakeElemVel(particles, wakelines, wakerings, panels[i].cpt)
        end
        
        # convert end wake rings into particles
        if wakelen > maxwakelen
            wakeend = wakerings[tewidth*(wakelen-1)+1:end]
            if length(te_neigh) == 0
                te_neigh, te_neighside, te_neighdir = calcneighbors(wakeend)
            end
            new_particles, new_wakelines = shedParticles(wakeend, te_neigh, te_neighdir, wakelines)

            wakelen = wakelen - 1
            wakerings = wakerings[1:tewidth*wakelen]
            
            # append wakelines to system 
            if length(wakelines) == tewidth
                for k = 1:length(wakelines)
                    wakelines[k] = new_wakelines[k]
                end
            else
                wakelines = new_wakelines
            end

            # append particles to system
            particles = cat(particles, new_particles, dims=1)
        end

        # pressure calculation
        for i = 1:length(panels)
            val = 0.0
            # i side, towards te
            if panels_neigh[1,i] > 0
                gam2 = panels[panels_neigh[1,i]].gam[1]
            else
                gam2 = 0.0
            end
            vcpt = panels[i].vcpt[:]
            val = val .+ dot(uinf.-vcpt+panels[i].wake_vel, panels[i].tani_uvec).* (panels[i].gam[1]-gam2)./panels[i].tani_len

            # j side, perpendicular to te direction
            if panels_neigh[4,i] > 0
                gam2 = panels[panels_neigh[4,i]].gam[1]
            else
                gam2 = 0.0
            end
            val = val .+ dot(uinf.-vcpt+panels[i].wake_vel, panels[i].tanj_uvec).* (panels[i].gam[1]-gam2)./panels[i].tanj_len

            val = val .+ panels[i].dgdt[1]
            panels[i].dp[1] = -rho*val[1]
            panels[i].df[:] = -panels[i].dp*panels[i].area[1].*panels[i].normal
        end

        # total forces
        total_force = [0;0;0]
        for i=1:length(panels)
            total_force = total_force .+ panels[i].df
        end

        Fz = total_force[3]
        println("Step: $t, Fz=$Fz")
        
        panels2vtk(panels, prefix * "_panels_$t.vtu")
        particles2vtk(particles, prefix * "_particles_$t.vtu")
        wakepanels2vtk(wakerings, prefix * "_wakerings_$t.vtu")

    end

end