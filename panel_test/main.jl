include("inf.jl")
include("aeroel.jl")
include("vis.jl")
include("linsys.jl")
import Base.Threads.@spawn

# solver input
uinf = [30; 0; 0]
Pinf = 0
rhoinf = 1.225
dt = 0.005
t_end = 100*dt
debug_dir = "./plane/debug/"
vis_dir = "./plane/vis/"

# Build reference frames
refs = []
ref_keys = []
newRef = refFrame("test", "0", [], false, false, [0;0;0], [0.9961947 0.0 -0.0871557; 0.0 1.0 0.0; 0.0871557 0.0 0.9961947]', [0;0;0], zeros(3,3))
# orient is defined as v_global = orient * v_local
refs = append!(refs, [newRef])
ref_keys = ["test"]

# Read geometry
@time begin
    panels, rr_all, wake_panels, rr_wake, wake_particles, end_vorts = get_geo("./dust_output/plane/geo_input.h5", uinf, refs, ref_keys)
end

function run(panels, rr_all, wake_panels, rr_wake, wake_particles, end_vorts)



t = 0.0
step_num = Int(1)

nWakePan = size(wake_panels, 1)
nPan = size(panels, 1)

lastpanidou = zeros(nWakePan) #Last vortex intensity from removed panels
endpanidou = zeros(nWakePan) #Last vortex intensity from removed panels

# visualize
panels2vtk(panels, rr_all, "mesh_$step_num.vtu", vis_dir)
panels2vtk(wake_panels, rr_wake, "wake_pan_$step_num.vtu", vis_dir)
particles2vtk(wake_particles, "wake_particles_$step_num.vtu", vis_dir)
step_num = step_num + 1


# Timestepping
while t < t_end
    # update uvort, influence of vortex particles and end vortex line
    panels = update_panel_uvorts(panels, wake_particles, end_vorts)
#=
    for i = 1:size(wake_panels,1)
        vel = uvortParticles(wake_particles, panels[i].center)
        vel2 = uvortVortLines(end_vorts, panels[i].center)
        vel = vel .+ vel2
        wake_panels[i].velVort[:] = vel[:]
    end
=#  


# panel influence matrix
@time begin
A, B, RHS = panel_influence(panels, rr_all)
end

# influence of trailing edge
A = te_influence(wake_panels, rr_wake, panels, A)

# assemble RHS (add uinf and uvort contribution)
RHS = calc_RHS(panels, B, RHS)

# debug: output
writedlm("B_static.csv", B)

# solve linear system
solution = A\RHS
# debug: output
debugMatrices(A, RHS, solution, step_num-1, debug_dir)

# update strengths
for i = 1:size(panels, 1)
    panels[i].mag[1] = solution[i]
end
for i = 1:size(wake_panels,1)
    wake_panels[i].mag[1] = panels[wake_panels[i].panIdx[1]].mag[1] - panels[wake_panels[i].panIdx[2]].mag[1]
    #println(wake_panels[i].mag[1])
end

# calculate velocities at particles and evolve by dt
@time begin
wake_particles = update_particles(panels, wake_panels, wake_particles, end_vorts, rr_all, rr_wake, uinf, dt)
end

# shed new particles
pts1 = zeros(3,nWakePan) # temporary vector to hold end vortex points
pts2 = zeros(3,nWakePan)
@Threads.threads for i = 1:size(wake_panels,1)
    posp1 = rr_wake[:,wake_panels[i].ee[3]]
    v1 = elemVel(panels, wake_panels, wake_particles, end_vorts, rr_all, rr_wake, posp1) .+ uinf
    posp2 = rr_wake[:,wake_panels[i].ee[4]]
    v2 = elemVel(panels, wake_panels, wake_particles, end_vorts, rr_all, rr_wake, posp2) .+ uinf
    pts1[:,i] = posp1 .+ v1.*dt
    pts2[:,i] = posp2 .+ v2.*dt

    #println(posp1)
    #println(pts1[:,i])
    #println(posp2)
    #println(pts1[:,i])
    #println("------")

end

for i = 1:size(wake_panels,1)
    partvec = zeros(3)
    # left side
    #println(wake_panels[i].neigh_te)
    dir = rr_wake[:,wake_panels[i].ee[4]] .- pts2[:,i]
    if wake_panels[i].neigh_te[1] > 0
        ave = wake_panels[i].mag[1] - wake_panels[i].neigh_orient[1] * wake_panels[wake_panels[i].neigh_te[1]].mag[1]
        ave = ave/2
    else
        ave = wake_panels[i].mag[1]
    end
    partvec = partvec .+ dir.*ave

    # right side
    dir = -rr_wake[:,wake_panels[i].ee[3]] .+ pts1[:,i]
    if wake_panels[i].neigh_te[2] > 0
        ave = wake_panels[i].mag[1] - wake_panels[i].neigh_orient[2] * wake_panels[wake_panels[i].neigh_te[2]].mag[1]
        ave = ave/2
    else
        ave = wake_panels[i].mag[1]
    end
    partvec = partvec .+ dir.*ave

    # end side
    dir =  pts2[:,i] .- pts1[:,i]
    ave = wake_panels[i].mag[1]-lastpanidou[i]
    lastpanidou[i] = wake_panels[i].mag[1]
    partvec = partvec .+ dir.*ave
    
    # calculate the center
    posp = (pts1[:,i] .+ pts2[:,i] .+ rr_wake[:,wake_panels[i].ee[4]].+ rr_wake[:,wake_panels[i].ee[3]])./4

    # add wake particle
    cent = posp
    mag = norm(partvec)
    if wake_panels[i].mag[1] > 1e-13
        dir = partvec./mag
    else
        dir = partvec
    end
    vel = [0;0;0] # panel doesn't move
    new_part = wake_part(dir, [mag], cent, vel)
    #println(cent)
    push!(wake_particles, new_part)
    endpanidou[i] = wake_panels[i].mag[1]

    # attach end vortex
    ee = [wake_panels[i].ee[4]; wake_panels[i].ee[3]]
    p1 = rr_wake[:,wake_panels[i].ee[4]]
    p2 = rr_wake[:,wake_panels[i].ee[3]]
    mag = lastpanidou[i]
    mag = [mag]
    ver_vel = [0;0;0]
    center, n_ver, n_sides, edge_vec, edge_len, edge_uni = getLineProp([p1 p2])
    newvortline = vortline(ee,  [p1 p2], center, edge_vec, edge_uni, edge_len, ver_vel, mag)

    if size(end_vorts,1) < size(wake_panels,1)
        push!(end_vorts, newvortline)
    else
        end_vorts[i] = newvortline
    end

end

# visualize
panels2vtk(panels, rr_all, "mesh_$step_num.vtu", vis_dir)
panels2vtk(wake_panels, rr_wake, "wake_pan_$step_num.vtu", vis_dir)
particles2vtk(wake_particles, "wake_particles_$step_num.vtu", vis_dir)

# attach end vortex

t = t + dt
step_num = step_num + Int(1)

println("Time: $t")
end # timestepping
end

run(panels, rr_all, wake_panels, rr_wake, wake_particles, end_vorts)