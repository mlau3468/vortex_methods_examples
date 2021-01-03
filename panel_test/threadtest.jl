include("vis.jl")
include("inf.jl")
include("aeroel.jl")
# Read input file
fname = "./dust_output/plane/geo_input.h5"
fid = h5open(fname, "r")

display(fid["Components"])
display(read(fid["Components"]["Comp003"]["Geometry"]["neigh"]))

# Build reference frames
refs = []
ref_keys = []
newRef = refFrame("test", "0", [], false, false, [0;0;0], [0.9961947 0.0 -0.0871557; 0.0 1.0 0.0; 0.0871557 0.0 0.9961947]', [0;0;0], zeros(3,3))
# orient is defined as v_global = orient * v_local
refs = append!(refs, [newRef])
ref_keys = ["test"]

# solver input
uinf = [30; 0; 0]
Pinf = 0
rhoinf = 1.225
dt = 0.005
t_end = 100*dt
debug_dir = "./debug2/"
vis_dir = "./vis2/"
t = 0.0
step_num = Int(1)

#=
# Read geometry
@time begin
    panels, rr_all, wake_panels, rr_wake, wake_particles, end_vorts = get_geo("./dust_output/plane/geo_input.h5", uinf, refs, ref_keys)
end



# visualize
panels2vtk(panels, rr_all, "mesh_$step_num.vtu", vis_dir)
panels2vtk(wake_panels, rr_wake, "wake_pan_$step_num.vtu", vis_dir)
particles2vtk(wake_particles, "wake_particles_$step_num.vtu", vis_dir)
=#