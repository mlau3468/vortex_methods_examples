using Interpolations
using DelimitedFiles
using Plots
using Printf

# strongly coupled algorithm from Gallay, Nonlinear Lifting Line for PrePost-stall Flows.
# Applied to a discrete wing modelled using vortex lattices
# Borrows geometry and visualization functions from VPM code

include("geo.jl")
include("vis.jl")
include("aeroel.jl")

# operating conditions
alpha = 5
U = 50
uinf = [U*cos(deg2rad(alpha)); 0; U*sin(deg2rad(alpha))]
rho = 1.225

# read in airfoil data
c81 = readdlm("naca0012.csv", ',', Float64)
cl_interp =  LinearInterpolation(c81[:,1], c81[:,2])
cd_interp =  LinearInterpolation(c81[:,1], c81[:,3])

# read in panels
panVert, panCon = createRect(10, 1, 20, 1,[0;0;0],0)
panCpt, panNorm, panNPt, panEdgeVec, panEdgeLen, panEdgeUVec, panArea, panTang, panSinTi, panCosTi = calcPanProps(panVert, panCon)
panNeighIdx, panNeighSide, panNeighDir, panNNeigh = calcneighbors(panCon, panNPt)


npts = size(panVert,2)
npan = size(panCon,2)

panPres = zeros(Float32, npan)
panVelSurf = zeros(Float32, 3, npan)
panGam = zeros(Float64, npan)

A = zeros(Float64, npan, npan)
RHS = zeros(Float64, npan)


for i = 1:npan
    for j = 1:npan
        # influence of jth panel on ith collocation point
        vel = vrtxring(panVert[:,panCon[:,j]], panCpt[:,i], 1.0)
        A[i,j] = dot(vel, panNorm[:,i])
    end
end

# build rhs vector
for i = 1:npan
    RHS[i] = -dot(uinf, panNorm[:,i])
    #RHS[i] = RHS[i] - dot(panels[i].wake_vel, panNorm[:,i])
    #RHS[i] = RHS[i] + dot(panels[i].vcpt, panNorm[:,i])
end

# solve matrix for panel gamma
sol = A\RHS
for i = 1:npan
    panGam[i] = sol[i]
end


pan2Vtu(panCon, panVert, panGam, panPres, panVelSurf, "test.vtu")