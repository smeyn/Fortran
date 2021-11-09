"""
Evaluate a spline for a grid.
Scratch pad buffers are preaocated to reduce garbage
collection.     

# Arguments:
- `spline::Spline`: the bspline type to use 
- `x:vector{Float64}` :  x coordinates of grid
- `y:vector{Float64}` :  y coordinates of grid
- `z:Matrix{Float64}` :  output! grid results
- `wx:Matrix{Float64}` :  scratch pad
- `wy:Matrix{Float64}` :  scratchpad
- `lx:Vector{Int}` : scratch pad
- `ly:Vector{Int}` : scratch pad
"""
function fpbisp!(spline::Spline, x::Vector{Float64}, y::Vector{Float64}, z::Matrix{Float64}, wx::Matrix{Float64}, wy::Matrix{Float64}, lx::Vector{Int}, ly::Vector{Int})
   
    kx1 = 0;ky1 = 0;l = 0;l1 = 0;l2 = 0;m = 0;nkx1 = 0;nky1 = 0
    arg = 0.0;sp = 0.0;tb = 0.0;te = 0.0
    h = zeros(Float64, 6)
 
    mx = length(x)
    my = length(y)
    kx1 = spline.kx + 1
    nx = spline.nx
    nkx1 = nx - kx1
    tb = spline.tx[kx1]
    te = spline.tx[nkx1 + 1]
    l = kx1
    l1 = l + 1
    for i in 1:mx
        arg = x[i]
        if (arg < tb) arg = tb end
        if (arg > te) arg = te end
        while !(arg < spline.tx[l1] || l == nkx1)
            l = l1
            l1 = l + 1
        end        
        fpbspl!(spline.tx, nx, spline.kx, arg, l, h)
        lx[i] = l - kx1
        for j in 1:kx1 
            wx[i,j] = h[j]
        end        
    end
    ky1 = spline.ky + 1
    ny = spline.ny
    nky1 = ny - ky1
    tb = spline.ty[ky1]
    te = spline.ty[nky1 + 1]
    l = ky1
    l1 = l + 1
    for i in 1:my 
        arg = y[i]
        if (arg < tb) arg = tb end
        if (arg > te) arg = te end       
        while !(arg < spline.ty[l1] || l == nky1)
            l = l1
            l1 = l + 1
        end          
        fpbspl!(spline.ty, ny, spline.ky, arg, l, h)
        ly[i] = l - ky1
        for j in 1:ky1
            wy[i,j] = h[j]            
        end
     end
    for i in 1:mx 
        l = lx[i] * nky1
        for i1 in 1:kx1 
            h[i1] = wx[i,i1]
        end
        for j in 1:my
            l1 = l + ly[j]
            sp = 0.
            for i1 in 1:kx1 
                l2 = l1
                for j1 in 1:ky1
                    l2 = l2 + 1
                    sp = sp + spline.c[l2] * h[i1] * wy[j,j1]
                end
                l1 = l1 + nky1
            end
            z[i,j] = sp
        end
    end
    return
end

    