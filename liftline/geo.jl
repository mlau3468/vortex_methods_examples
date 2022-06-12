using Interpolations

"""
    cosine_spacing(startval, endval, num, weight, spacetype)
Cosine space num numbers between startval and endval.
"""
function cosine_spacing(startval, endval, num, weight, spacetype)
    if spacetype == "cosOB"
        w = weight
        k = sin.(LinRange(0,pi/2,num))
        k2 = LinRange(0,1,num)
        b = k.*w .+ k2.*(1-w)
        res = startval .+ (endval-startval).*b
    elseif spacetype == "cosIB"
        w = weight
        k = 1 .- sin.(LinRange(pi/2,pi,num))
        k2 = LinRange(0,1,num)
        b = k.*w .+ k2.*(1-w)
        res = startval .+ (endval-startval).*b
    elseif spacetype == "equal"
        res = LinRange(startval, endval, num)
    elseif spacetype == "cos"
        w = weight
        beta = LinRange(0, pi, num)
        k = (1 .-cos.(beta))/2
        k2 = LinRange(0,1,num)
        b = k.*w .+ k2.*(1-w)
        res =startval .+ (endval-startval).*b
    end
    return res
end

function refine_wing(xloc, yloc, zloc, chord, pitch, nspan::Int, spanspace::String)
    nsec = length(xloc)

    dist = zeros(nsec)
    for i in 2:nsec
        dist[i] = dist[i-1] + sqrt((xloc[i]-xloc[i-1])^2 + (yloc[i]-yloc[i-1])^2 + (zloc[i]-zloc[i-1])^2 )
    end

    # interpolations for properties
    xref_int = LinearInterpolation(dist, xloc)
    yref_int = LinearInterpolation(dist, yloc)
    zref_int = LinearInterpolation(dist, zloc)
    chord_int = LinearInterpolation(dist, chord)
    pitch_int = LinearInterpolation(dist, pitch)

    # build where to query interpolations
    distlocs = cosine_spacing(dist[1], dist[end], nspan+1, 0.8, spanspace)

    # query interpolations
    new_xloc = xref_int(distlocs)
    new_yloc = yref_int(distlocs)
    new_zloc = zref_int(distlocs)
    new_chord = chord_int(distlocs)
    new_pitch = pitch_int(distlocs)
    return new_xloc, new_yloc, new_zloc, new_chord, new_pitch
end

function refine_wing(xloc, yloc, zloc, chord, pitch, nspan::Array{Int})
    nsec = length(xloc)

    dist = zeros(nsec)
    for i in 2:nsec
        dist[i] = dist[i-1] + sqrt((xloc[i]-xloc[i-1])^2 + (yloc[i]-yloc[i-1])^2 + (zloc[i]-zloc[i-1])^2 )
    end

    # interpolations for properties
    xref_int = LinearInterpolation(dist, xloc)
    yref_int = LinearInterpolation(dist, yloc)
    zref_int = LinearInterpolation(dist, zloc)
    chord_int = LinearInterpolation(dist, chord)
    pitch_int = LinearInterpolation(dist, pitch)

    # build where to query interpolations
    distlocs = zeros(sum(nspan) + 1)
    idx = 1
    for i = 2:nsec
        new_distloc = collect(LinRange(dist[i-1],dist[i],nspan[i-1]+1))
        if i == nsec
            distlocs[idx:idx+nspan[i-1]] .= new_distloc
        else
            distlocs[idx:idx+nspan[i-1]-1] .= new_distloc[1:end-1]
            idx += nspan[i-1]
        end
    end

    # query interpolations
    new_xloc = xref_int(distlocs)
    new_yloc = yref_int(distlocs)
    new_zloc = zref_int(distlocs)
    new_chord = chord_int(distlocs)
    new_pitch = pitch_int(distlocs)
    return new_xloc, new_yloc, new_zloc, new_chord, new_pitch
end

function wing2vortlat(xloc, yloc, zloc, chord, pitch, nchord, chordspace, xref)
    floatType = Float32
    intType = Int32

    nspan = length(xloc) - 1 # number of spanwise elements
    npt = (nspan+1)*(nchord+1) # number of points
    npan = nspan*nchord # number of panels

    # create nondimensional chordwise points from 0 to 1 
    x_rel = cosine_spacing(0, 1, nchord+1, 0.75, chordspace)

    # Build points and connectivity
    pan_vert = zeros(floatType, 3, npt)
    pan_con = zeros(intType, 4, npan)
    count = 1
    for i = 1:(nchord+1)
        for j = 1:(nspan+1)
            new_x = xloc[j] + (x_rel[i]-xref)*chord[j]*cosd(pitch[j])
            new_z = zloc[j] - (x_rel[i]-xref)*chord[j]*sind(pitch[j])
            pan_vert[1,count] = new_x
            pan_vert[2,count] = yloc[j]
            pan_vert[3,count] = new_z
            count = count + 1
        end
    end

    # Build panel connectivity
    # Follow counter clockwise paneling
    for i = 1:nchord
        for j = 1:nspan
            pan_con[1,(i-1)*(nspan)+j] = (i-1)*(nspan+1)+j
            pan_con[2,(i-1)*(nspan)+j] = (i)*(nspan+1)+j
            pan_con[3,(i-1)*(nspan)+j] = (i)*(nspan+1)+j+1
            pan_con[4,(i-1)*(nspan)+j] = (i-1)*(nspan+1)+j+1
        end
    end

    println("Number of Panels: $npan")
    println("Number of Vertices: $npt")

    return (pan_vert=pan_vert, pan_con=pan_con)
end

function wing2liftline(xle, yle, zle, chord, pitch, tevec, telen)
    floatType = Float32
    intType = Int32

    nspan = length(xle) - 1 # number of spanwise elements
    npt = (nspan+1)*3 # number of points
    npan = nspan # number of panels

    xref = 0.25

    # Build points and connectivity
    pan_vert = zeros(floatType, 3, npt)
    count = 1
    pan_con = zeros(intType, 6, npan)

    for j = 1:nspan+1
        # leading edge point
        ple_x = xle[j] + (xref)*chord[j]*cosd(pitch[j])
        ple_z = zle[j] - (xref)*chord[j]*sind(pitch[j])
        ple_y = yle[j]

        # trailing edge point
        pte_x = xle[j] + (1+xref)*chord[j]*cosd(pitch[j])
        pte_z = zle[j] - (1+xref)*chord[j]*sind(pitch[j])
        pte_y = yle[j]

        # wake point
        pwake_x = pte_x + tevec[1]*telen
        pwake_z = pte_z + tevec[2]*telen
        pwake_y = pte_y + tevec[3]*telen

        # store points
        pan_vert[1,count] = ple_x
        pan_vert[2,count] = ple_y
        pan_vert[3,count] = ple_z
        pan_vert[1,count+1] = pte_x
        pan_vert[2,count+1] = pte_y
        pan_vert[3,count+1] = pte_z
        pan_vert[1,count+2] = pwake_x
        pan_vert[2,count+2] = pwake_y
        pan_vert[3,count+2] = pwake_z
        count += 3
    end

    # Build panel connectivity
    # Follow counter clockwise paneling
    for i = 1:nspan
        pan_con[1,i] = 3*i + 3 # wake
        pan_con[2,i] = 3*i + 2 # te 
        pan_con[3,i] = 3*i + 1 # le
        pan_con[4,i] = 3*(i-1) + 1 # le
        pan_con[5,i] = 3*(i-1) + 2 # te
        pan_con[6,i] = 3*(i-1) + 3 # wake
    end

    return (pan_vert=pan_vert, pan_con=pan_con)
end

function calc_liftline_props(pan_vert, pan_con)
    floatType = Float32

    npan = size(pan_con, 2)
    pan_cpt = zeros(floatType, 3, npan) # collocation points
    pan_norm = zeros(floatType, 3, npan)
    pan_area = zeros(floatType, npan)
    pan_edgeuni = zeros(floatType, 3, 5, npan)
    pan_edgelen = zeros(floatType, 5, npan)
    pan_edgevec = zeros(floatType, 3, 5, npan)

    for i = 1:npan
        pan_cpt[:,i] = (pan_vert[:,pan_con[2,i]] .+ pan_vert[:,pan_con[3,i]] .+ pan_vert[:,pan_con[4,i]] .+ pan_vert[:,pan_con[5,i]])./4
        v1 = pan_vert[:,pan_con[4,i]] .- pan_vert[:,pan_con[2,i]]
        v2 = pan_vert[:,pan_con[5,i]] .- pan_vert[:,pan_con[3,i]]
        v3 = cross(v1,v2)
        pan_norm[:,i] = v3./norm(v3)
        pan_area[i] = 0.5*norm(v3)
        for j = 1:5
            pan_edgevec[:,j,i] .= pan_vert[:,pan_con[j+1,i]] .- pan_vert[:,pan_con[j,i]]
            pan_edgelen[j,i] = norm(pan_edgevec[:,j,i])
            pan_edgeuni[:,j,i] .= pan_edgevec[:,j,i]./pan_edgelen[j,i]
        end
    end

    return pan_cpt, pan_norm, pan_area, pan_edgeuni, pan_edgelen
end