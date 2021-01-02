function panel_influence(panels, rr_all)
    nPan = size(panels,1)
    A = zeros(nPan, nPan)
    B = zeros(nPan, nPan) #Bstatic
    RHS = zeros(nPan)
    @inbounds @Threads.threads for i = 1:nPan
    @inbounds for j = 1:nPan
            if i == j
                dou = -2*pi
            else
                dou = dub(panels[j], panels[i].center, rr_all)
            end
            A[i,j] = -dou
    
            sou = sourc(panels[j], panels[i].center, dou, rr_all)
            B[i,j] = sou
        end
        RHS[i] = 0.0
    end
    return A, B, RHS
end

function te_influence(wake_panels, rr_wake, panels, A)
    @Threads.threads for i = 1:size(panels,1)
        # go through each wake panel
        for j = 1:size(wake_panels,1)
            # effect of wake panel on panel i
            a = -dub(wake_panels[j], panels[i].center, rr_wake)
            #=Modify influence coefficients of trailing edge panels 
            associated with this trailing edge wake on panel i =#
            A[i, wake_panels[j].panIdx[1]] = A[i, wake_panels[j].panIdx[1]] + a
            A[i, wake_panels[j].panIdx[2]] = A[i, wake_panels[j].panIdx[2]] - a
        end
    end
    return A
end

function calc_RHS(panels, B, RHS)
    @inbounds @Threads.threads for i = 1:size(panels,1)
    @inbounds for j = 1:size(panels,1)
            RHS[i] = RHS[i] + B[i,j] .* sum(panels[j].norm.*(-uinf.-panels[j].velVort))
        end
    end
    return RHS
end