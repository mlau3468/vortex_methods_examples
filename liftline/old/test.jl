include("aeroel.jl")
x1 = 0
y1 = 1
z1 = 5
x2 = 2
y2 = -5
z2 = 3

x = 1
y = 7
z = 8

vel2 = zeros(3)
vortxl!(x,y,z,x1,y1,z1,x2,y2,z2, 1, vel2)
println(vel2)

vel1 = vrtxline([x1;y1;z1], [x2;y2;z2], [x;y;z], 1)
println(vel1)