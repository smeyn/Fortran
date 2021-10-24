"""evaluate a spline for a grid  """
function fpbisp(spline::Spline, x::Vector, y::Vector, z::Matrix, wx::Matrix, wy::Matrix, lx::Vector, ly::Vector)
   #= 
   c  ..scalar arguments..
         integer nx,ny,kx,ky,mx,my
   c  ..array arguments..
      integer lx(mx),ly(my)
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),
     * wx(mx,kx+1),wy(my,ky+1) =#
    # ..local scalars..
    #    integer kx1,ky1,l,l1,l2,m,nkx1,nky1
    kx1 = 0;ky1 = 0;l = 0;l1 = 0;l2 = 0;m = 0;nkx1 = 0;nky1 = 0
    #    real arg,sp,tb,te
    arg = 0.0;sp = 0.0;tb = 0.0;te = 0.0
    # c  ..local arrays..
    #      real h(6)
    h = zeros(Float64, 6)
    # c  ..subroutine references..
    # c    fpbspl
    # c  ..
    mx = length(x)
    my = length(y)
    kx1 = spline.kx + 1
    nx = spline.nx
    nkx1 = nx - kx1
    tb = spline.tx[kx1]
    te = spline.tx[nkx1 + 1]
    l = kx1
    l1 = l + 1
    for i in 1:mx # do 40 i=1,mx
        arg = x[i]
        if (arg < tb) arg = tb end
        if (arg > te) arg = te end
        @label L10
        if (arg < spline.tx[l1] || l == nkx1) @goto L20 end
        l = l1
        l1 = l + 1
        @goto L10
        @label L20 
        fpbspl(spline.tx, nx, spline.kx, arg, l, h)
        lx[i] = l - kx1
        for j in 1:kx1 # do 30 j=1,kx1
            wx[i,j] = h[j]

            @label L30#    continue
        end
        @label L40 # continue
    end
    ky1 = spline.ky + 1
    ny = spline.ny
    nky1 = ny - ky1
    tb = spline.ty[ky1]
    te = spline.ty[nky1 + 1]
    l = ky1
    l1 = l + 1
    for i in 1:my # do 80 i=1,my
        arg = y[i]
        if (arg < tb) arg = tb end
        if (arg > te) arg = te end
        @label L50
        if (arg < spline.ty[l1] || l == nky1) @goto L60 end
        l = l1
        l1 = l + 1
        @goto L50
        @label L60 
        fpbspl(spline.ty, ny, spline.ky, arg, l, h)
        ly[i] = l - ky1
        for j in 1:ky1# do 70 j=1,ky1
            wy[i,j] = h[j]
            @label L70 #   continue
        end
        @label L80 # continue
    end
    # m = 0
    for i in 1:mx # do 130 i=1,mx
        l = lx[i] * nky1
        for i1 in 1:kx1 # do 90 i1=1,kx1
            h[i1] = wx[i,i1]
            @label L90 #   continue
        end
        for j in 1:my # do 120 j=1,my
            l1 = l + ly[j]
            sp = 0.
            for i1 in 1:kx1 # do 110 i1=1,kx1
                l2 = l1
                for j1 in 1:ky1 # do 100 j1=1,ky1
                    l2 = l2 + 1
                    sp = sp + spline.c[l2] * h[i1] * wy[j,j1]
                    @label L100 #       continue
                end
                l1 = l1 + nky1
                @label L110 #     continue
            end
            # m = m + 1
            z[j,i] = sp
            @label L120 #   continue
        end
        @label L130 # continue
    end
    return
end

    