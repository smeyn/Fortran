
using Printf

  function prepareSpline(x, y, z, w, kx, ky,max_x_knots, max_y_knots;
      xb=nothing, xe=nothing, yb=nothing, ye=nothing, io=nothing)
      if isnothing(xb)  xb=min(x...) end
      if isnothing(xe)  xe=max(x...) end
      if isnothing(yb)  yb=min(y...) end
      if isnothing(ye)  ye=max(y...) end
      tx = zeros(Float64, max_x_knots)
      ty = zeros(Float64, max_y_knots)
      nx=0
      ny=0
      fp=0.0
      c = zeros(Float64, (max_x_knots-kx-1)*(max_y_knots-ky-1))
      return Spline(x, y, z, w,
                    xb, xe, yb, ye, 
                    kx, ky, nx, tx, ny,ty, c, fp, "")  
  end
  
function write_c(io, c::Vector, n)
    if isnothing(io) return end
    #io = IOBuffer()
    write(io, "debug: b-spline coefficients \n")
    for i in 1:n
        @printf(io, "%8.4f ",c[i])
        if mod(i,8) == 0
            write(io, "\n")
        end
    end
    write(io, "\n")       
end

function wrt_fpint(io, fpint, l1, l2, fpmax)
    if isnothing(io) return end
    @printf(io, "debug fpint=======\n")
    @printf(io, "debug l1, l2, fpmax %4i %4i %.4f\n", l1, l2, fpmax)
    for i in 1:5
        @printf(io, "%.4f ", fpint[i])
    end
    write(io, '\n')
end

function wrt_index(io,index, nreg)
    if isnothing(io) return end
    @printf(io, "debug index =======\n")
    @printf(io, "debug %4i\n", nreg)
    for i in 1:nreg
        @printf(io, "%4i ", index[i])
    end
    write(io, '\n')
end   



function write_sp(io,sp, in, k1)
    if isnothing(io) return end
    @printf(io, "debug sp =======\n")
    @printf(io, "debug in=%4i  k1=%4i\n", in, k1)
    for i in 1:k1
        @printf(io, "%8.4f ", sp[in, i])
    end
    write(io, '\n')
end   

"""
calculate surface splines.

"""
function fpsurf!(iopt,m, spline::Spline,s,nxest,nyest,
        eta,tol,maxit, 
        wrk1::lwork1 ,wrk2::lwork2;io=nothing)

    #@info "fpsurf start"
    #@debug "wrk1.h= $(wrk1.h)"
    #@debug "nxest: $nxest, nyest:$nyest, nc:$nc. Size c:$(length(spline.c))"
    # set the scope of these variables
    iband = 0; iband1 = 0; iband3 = 0; iband4 = 0; 
    nreg = 0
    nk1x = 0
    nk1y = 0
    ncof = 0
    nminx = 0; nminy = 0; nreg = 0; 
    nyy = 0; nxx = 0
    rank = 0
    sq = 0.0
    fpms = 0

    hx = zeros(Float64, 6)
    hy = zeros(Float64, 6)
    one = 0.1e+01
    con1 = 0.1e0
    con9 = 0.9e0
    con4 = 0.4e-01
    half = 0.5e0
    ten = 0.1e+02
    v_sin = 0.0
    v_cos = 0.0
        #= 
        cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        c part 1: determination of the number of knots and their position.     c
        c ****************************************************************     c
        c given a set of knots we compute the least-squares spline sinf(x,y),  c
        c and the corresponding weighted sum of squared residuals fp=f(p=inf). c
        c if iopt=-1  sinf(x,y) is the requested approximation.                c
        c if iopt=0 or iopt=1 we check whether we can accept the knots:        c
        # if fp <=s we will continue with the current set of knots.          c
        # if fp > s we will increase the number of knots and compute the     c
        #    corresponding least-squares spline until finally  fp<=s.        c
        c the initial choice of knots depends on the value of s and iopt.      c
        # if iopt=0 we first compute the least-squares polynomial of degree  c
        #   kx in x and ky in y; nx=nminx=2*kx+2 and ny=nminy=2*ky+2.        c
        #   fp0=f(0) denotes the corresponding weighted sum of squared       c
        #   residuals                                                        c
        # if iopt=1 we start with the knots found at the last call of the    c
        #   routine, except for the case that s>=fp0; then we can compute    c
        #   the least-squares polynomial directly.                           c
        c eventually the independent variables x and y (and the corresponding  c
        c parameters) will be switched if this can reduce the bandwidth of the c
        c system to be solved.                                                 c
        cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        #ichang denotes whether(1) or not(-1) the directions have been inter-
        #changed. =#
    ichang = -1
    ier = 0
    x0 = spline.xb
    x1 = spline.xe
    y0 = spline.yb
    y1 = spline.ye
    kx = spline.kx
    ky = spline.ky#kyy
    kx1 = spline.kx + 1
    ky1 = spline.ky + 1
    nxe = nxest
    nye = nyest
    eps = sqrt(eta)
    
    if (iopt < 0) @goto  L20 end
    #calculation of acc, the absolute tolerance for the root of f(p)=s.
    acc =  tol*s
    if (iopt == 0) @goto  L10 end
    if (wrk1.fp0 > s) @goto  L20 end
    # initialization for the least-squares polynomial.
    @label L10 
    @info "IOPT >= 0"
    nminx = 2 * kx1
    nminy = 2 * ky1
    nx = nminx
    ny = nminy
    ier = -2
    @goto  L30  
    @label L20  
    @info "fp0 > s or iopt < 0"
    nx = spline.nx #nx0
    ny = spline.ny # ny0
   # main loop for the different sets of knots. m is a save upper bound
   # for the number of trials.

    @label L30  
    for iter in 1:m # do 420 iter=1,m
        @info "iter: $iter"
      # find the position of the additional knots which are needed for the
       # b-spline representation of s(x,y).
        l = nx
        for i in 1:kx1  # do 40 i=1,kx1
            spline.tx[i] = x0
            spline.tx[l] = x1
            l = l - 1
          
        end # @label L40    continue
        l = ny
        for i in 1:ky1 # do 50 i=1,ky1
            spline.ty[i] = y0
            spline.ty[l] = y1
            l = l - 1
        end#  @label L50    continue
        # find nrint, the total number of knot intervals and nreg, the number
        # of panels in which the approximation domain is subdivided by the
        # intersection of knots.
        nxx = nx - 2 * kx1 + 1
        nyy = ny - 2 * ky1 + 1
        nrint = nxx + nyy
        nreg = nxx * nyy
        @debug "nxx ($nxx) = nx($nx) - 2 * kx1($kx1) + 1"
        @debug "nyy ($nyy) = ny($ny) - 2 * ky1($ky1) + 1"
        @debug "nrint ($nrint) = $nxx + $nyy"
        @debug "nreg ($nreg) = $nxx * $nyy"
        # find the bandwidth of the observation matrix a.
        # if necessary, interchange the variables x and y, in order to obtain
        # a minimal bandwidth.
        iband1 = kx * (ny - ky1) + ky
        l = ky * (nx - kx1) + kx
        if (iband1 <= l) @goto  L130 end
        iband1 = l
        ichang = -ichang
        for i in 1:m # do 60 i=1,m
            store = spline.x[i]
            spline.x[i] = spline.y[i]
            spline.y[i] = store
        end #  @label L60    continue
        store = x0
        x0 = y0
        y0 = store
        store = x1
        x1 = y1
        y1 = store
        n = min(nx, ny)
        for i in 1:n # do 70 i=1,n
            store = spline.tx[i]
            spline.tx[i] = spline.ty[i]
            spline.ty[i] = store
        end#  @label L70    continue
        n1 = n + 1
        #  if(nx-ny) 80,120,100
        if (nx - ny) < 0 @goto L80
        elseif (nx - ny) == 0 @goto L120
        else @goto L100 end
        @label L80 
        for i in n1:ny#   do 90 i=n1,ny
            spline.tx[i] = spline.ty[i]
        end#  @label L90    continue
        @goto  L120
        @label L100  #  do 110 i=n1,nx
        for i in n1:nx
            spline.ty[i] = spline.tx[i]
        end#  @label L110    continue
        @label L120 
        l = nx
        nx = ny
        ny = l
        l = nxe
        nxe = nye
        nye = l
        l = nxx
        nxx = nyy
        nyy = l
        l = kx
        kx = ky
        ky = l
        kx1 = kx + 1
        ky1 = ky + 1
        @label L130
        iband = iband1 + 1
        # arrange the data points according to the panel they belong to.
        fporde!(spline.x, spline.y, m, kx, ky, spline.tx, nx, spline.ty, ny, wrk1.nummer, wrk1.index, nreg)
        # find ncof, the number of b-spline coefficients.
        wrt_index(io, wrk1.index, nreg)
        nk1x = nx - kx1
        nk1y = ny - ky1
        ncof = nk1x * nk1y
        #@debug "nk1x:$nk1x  nk1y:$nk1y  ncof:$ncof"
        # initialize the observation matrix a.
        for i in 1:ncof # do 140 i=1,ncof
            wrk1.f[i] = 0.
            for j in 1:iband # do 140 j=1,iband
                wrk1.a[i,j] = 0.
            end
        end #  @label L140    continue
        # initialize the sum of squared residuals.
        fp = 0.
        # fetch the data points in the new order. main loop for the
        # different panels.
        for num in 1:nreg # do 250 num=1,nreg
          # fix certain constants for the current panel; jrot records the column
          # number of the first non-zero element in a row of the observation
          # matrix according to a data point of the panel.
            num1 = num - 1
            lx = Int(trunc(num1 / nyy))
            #@debug "lx ($lx) = Int(trunc($num1 / $nyy))"
            l1 = lx + kx1
            ly = num1 - lx * nyy
            l2 = ly + ky1
            jrot = lx * nk1y + ly
                # test whether there are still data points in the panel.
            in = wrk1.index[num]
            @label L150  
            if (in == 0) @goto  L250 end 
          # fetch a new data point.
            wi = spline.w[in]
            zi = spline.z[in] * wi
            @debug "spline.z[in] * wi = $(spline.z[in]) * $wi"
          # evaluate for the x-direction, the (kx+1) non-zero b-splines at x(in).
            fpbspl!(spline.tx, nx, kx, spline.x[in], l1, hx)
            #@debug hx
          # evaluate for the y-direction, the (ky+1) non-zero b-splines at y(in).
            fpbspl!(spline.ty, ny, ky, spline.y[in], l2, hy)
            #@debug hy
          # store the value of these b-splines in spx and spy respectively.
            for i in 1:kx1 # do 160 i=1,kx1
                #@debug "wrk1.spx[$in,$i] = $(hx[i])"
                wrk1.spx[in,i] = hx[i]
            end #  @label L160      continue
            for i in 1:ky1# do 170 i=1,ky1
                #@debug "wrk1.spy[$in,$i] = $(hy[i])"
                wrk1.spy[in, i] = hy[i]
            end#  @label L170      continue
          # initialize the new row of observation matrix.
            for i in 1:iband # do 180 i=1,iband
                wrk1.h[i] = 0.
            end#  @label L180      continue
          # calculate the non-zero elements of the new row by making the cross
          # products of the non-zero b-splines in x- and y-direction.
            i1 = 0
            for i in 1:kx1 # do 200 i=1,kx1
                hxi = hx[i]
                j1 = i1
                for j in 1:ky1 # do 190 j=1,ky1
                    j1 = j1 + 1
                    wrk1.h[j1] = hxi * hy[j] * wi
                    #@debug "wrk1.h= $(wrk1.h)"
                end#  @label L190        continue
                i1 = i1 + nk1y
                @label L200  #    continue
            end

            #@debug "rotate the row into triangle by givens transformations ."
            # rotate the row into triangle by givens transformations .
            irot = jrot
            for i in 1:iband # do 220 i=1,iband
                irot = irot + 1
                piv = wrk1.h[i]
                if (piv == 0.0) @goto  L220 end
                # calculate the parameters of the givens transformation.
                v_cos, v_sin, wrk1.a[irot, 1] = fpgivs(piv, wrk1.a[irot, 1])
                # apply that transformation to the right hand side.
                zi, wrk1.f[irot] = fprota(v_cos, v_sin, zi, wrk1.f[irot])
                if (i == iband) @goto  L230 end
                # apply that transformation to the left hand side.
                i2 = 1
                i3 = i + 1
                for j in i3:iband # do 210 j=i3,iband
                    i2 = i2 + 1
                    wrk1.h[j], wrk1.a[irot, i2] = fprota(v_cos, v_sin, wrk1.h[j], wrk1.a[irot, i2])
                end#  @label L210        continue
                @label L220 
            end #        @label L220      continue
            #@debug "add the contribution of the row to the sum of squares of residual right hand sides."
            # add the contribution of the row to the sum of squares of residual
            # right hand sides.
            @label L230  
            @debug "zi[in=$in]=$zi"
            fp = fp + zi * zi
           # @debug "find the number of the next data point in the panel."
            # find the number of the next data point in the panel.
            @label L240     
            in = wrk1.nummer[in]
            @goto  L150
            @label L250  #  continue
        end
        @debug "fp=$fp"
        #@debug "find dmax, the maximum value for the diagonal elements in the reduced traingle"
        # find dmax, the maximum value for the diagonal elements in the reduced
        # triangle.
        dmax = 0.
        #@debug "ncof: $ncof"
        for i in 1:ncof # do 260 i=1,ncof
            if (wrk1.a[i, 1] <= dmax) @goto  L260 end
            dmax = wrk1.a[i, 1]
            @label L260   # continue
        end
       # @debug "check whether the observation matrix is rank deficient."
        # check whether the observation matrix is rank deficient.
        sigma = eps * dmax
        for i in 1:ncof # do 270 i=1,ncof
            if (wrk1.a[i, 1] <= sigma) @goto  L280 end
        end #  @label L270    continue
        # backward substitution in case of full rank.
        #@debug "backward substitution in case of full rank."
        fpback!(wrk1.a, wrk1.f, ncof, iband, spline.c)
        #@debug "after fpback c:" spline.c[1:ncof]
        write_c(io, spline.c, ncof)
        rank = ncof
        for i in 1:ncof #        do 275 i=1,ncof
            wrk1.q[i,1] = wrk1.a[i,1] / dmax
            @label L275   # continue
        end
        @goto  L300
        # in case of rank deficiency, find the minimum norm solution.
        # check whether there is sufficient working space
        @label L280   
        #@debug "in case of rank deficiency, find the minimum norm solution."
        lwest = ncof * iband + ncof + iband
        #TODO Fix this size check
        #if (lwrk < lwest) @goto  L780 end
        for i in 1:ncof# do 290 i=1,ncof
            wrk1.ff[i] = wrk1.f[i]
            for j in 1:iband # do 290 j=1,iband
                wrk1.q[i, j] = wrk1.a[i, j]
                @label L290 #   continue
            end
        end
        lf = 1
        lh = lf + ncof
        la = lh + iband
        #@debug "Observation matrix" wrk1.q[1:ncof, 1:iband]
        sq, rank = fprank!(wrk1.q, wrk1.ff, ncof, iband,  sigma, spline.c,  wrk2)#wrk(la), wrk(lf), wrk(lh))
        #@debug "after fprank c:" spline.c[1:ncof]
        write_c(io, spline.c, ncof)
        #@debug "Observation matrix psot rank" wrk1.q[1:ncof, 1:iband]
        
        #@debug "fprank determined rank=$rank with $sq as squared residual"
        for i in 1:ncof # do 295 i=1,ncof
            wrk1.q[i, 1] = wrk1.q[i, 1] / dmax
        end #  @label L295    continue
        # add to the sum of squared residuals, the contribution of reducing
        # the rank.
        fp = fp + sq
        @label L300  
       
        #@debug "L300"
        if (ier == (-2))
             wrk1.fp0 = fp end
        # test whether the least-squares spline is an acceptable solution.
        if (iopt < 0) @goto  L820 end
        @debug "fp=$fp"
        fpms = fp - s
        if (abs(fpms) <= acc) # if(fp) 815,815,820
            if fp < 0 @goto L815
            elseif fp == 0 @goto L815
                else @goto L820
            end        
        end
        # test whether we can accept the choice of knots.
        if (fpms < 0.) 
            @info "fpms($fpms) < 0.0"
            @goto  L430 
        end
        # test whether we cannot further increase the number of knots.
        if (ncof > m)
            #@debug "ncof ($ncof) > m ($m)"
            #@debug "nk1x: $nk1x"

            @goto  L790 
        end
        ier = 0
        # search where to add a new knot.
        # find for each interval the sum of squared residuals fpint for the
        # data points having the coordinate belonging to that knot interval.
        # calculate also coord which is the same sum, weighted by the position
        # of the data points considered.
        @label L310 
        for i in 1:nrint #   do 320 i=1,nrint
            wrk1.fpint[i] = 0.
            wrk1.coord[i] = 0.
                #  @label L320    continue
        end
        for num in 1:nreg # do 360 num=1,nreg
            num1 = num - 1
            lx = Int(trunc(num1 / nyy))
            l1 = lx + 1
            ly = num1 - lx * nyy
            l2 = ly + 1 + nxx
            jrot = lx * nk1y + ly
            in = wrk1.index[num]
            @label L330  
            if (in == 0) @goto  L360 end
           # @debug "wrk1.spx[$in,1:$kx1] $(wrk1.spx[in,1:kx1])"
            write_sp(io, wrk1.spx, in, kx1)
            #@debug "wrk1.spy[$in,1:$ky1] $(wrk1.spy[in,1:ky1])"
            write_sp(io, wrk1.spy, in, ky1)
            store = 0.
            i1 = jrot
            for i in 1:kx1 # do 350 i=1,kx1
                hxi = wrk1.spx[in,i]
                j1 = i1
                for j in 1:ky1 # do 340 j=1,ky1
                    j1 = j1 + 1
                    store = store + hxi * wrk1.spy[in, j] * spline.c[j1]
                    @label L340   #     continue
                end
                i1 = i1 + nk1y
                @label L350  #    continue
            end
            store = (spline.w[in] * (spline.z[in] - store))^2
            wrk1.fpint[l1] = wrk1.fpint[l1] + store
            wrk1.coord[l1] = wrk1.coord[l1] + store * spline.x[in]
            wrk1.fpint[l2] = wrk1.fpint[l2] + store
            wrk1.coord[l2] = wrk1.coord[l2] + store * spline.y[in]
            in = wrk1.nummer[in]
            @goto  L330
            @label L360  #  continue
        end
        # find the interval for which fpint is maximal on the condition that
        # there still can be added a knot.
        @label L370    
        l = 0
        fpmax = 0.
        l1 = 1
        l2 = nrint
        if (nx == nxe) l1 = nxx + 1 end
        if (ny == nye) l2 = nxx end
        @debug "prior knot increase nx:$nx, nxe:$nxe l1:$l1  ny:$ny nye:$nye  l2:$l2"
        if (l1 > l2) @goto  L810 end
        for i in l1:l2 # do 380 i=l1,l2
            if (fpmax >= wrk1.fpint[i]) @goto  L380 end
            l = i
            fpmax = wrk1.fpint[i]
            @label L380  #  continue
        end
        wrt_fpint(io,wrk1.fpint,  l1, l2,  fpmax)
        @debug "l1=$l1, l2=$l2, fpmax=$fpmax"
        @debug wrk1.fpint[l1:l2]
        # test whether we cannot further increase the number of knots.
        if (l == 0) @goto  L785 end
        # calculate the position of the new knot.
        arg = wrk1.coord[l] / wrk1.fpint[l]
        # test in what direction the new knot is going to be added.
        if (l > nxx) @goto  L400 end
        # addition in the x-direction.
        jxy = l + kx1
        wrk1.fpint[l] = 0.
        fac1 = spline.tx[jxy] - arg
        fac2 = arg - spline.tx[jxy - 1]
        if ((fac1 > (ten * fac2)) || (fac2 > (ten * fac1))) @goto  L370 end
        j = nx
        for i in jxy:nx # do 390 i=jxy,nx
            spline.tx[j + 1] = spline.tx[j]
            j = j - 1
            @label L390 #   continue
        end
        spline.tx[jxy] = arg
        @debug "New knot in x at $jxy = $arg"
        nx = nx + 1
        @goto  L420
        # addition in the y-direction.
        @label L400
        jxy = l + ky1 - nxx
        wrk1.fpint[l] = 0.
        fac1 = spline.ty[jxy] - arg
        fac2 = arg - spline.ty[jxy - 1]
        if ((fac1 > (ten * fac2)) || ((fac2 > (ten * fac1)))) @goto  L370 end
        j = ny
        for i in jxy:ny # do 410 i=jxy,ny
            spline.ty[j + 1] = spline.ty[j]
            j = j - 1
                                #   @label L410   # continue
        end
        spline.ty[jxy] = arg
        @debug "New knot in y at $jxy = $arg"
        ny = ny + 1
        # restart the computations with the new set of knots.
        @label L420  # continue
    end
   # test whether the least-squares polynomial is a solution of our
   # approximation problem.
    @label L430  
    @debug "@ label 430"
    if (ier == (-2)) 
        @debug "L430. as ier = -2 jumping to L830. nk1x = $nk1x"
        @goto  L830 end
        #= 
        cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        c part 2: determination of the smoothing spline sp(x,y)                c
        c *****************************************************                c
        c we have determined the number of knots and their position. we now    c
        c compute the b-spline coefficients of the smoothing spline sp(x,y).   c
        c the observation matrix a is extended by the rows of a matrix,        c
        c expressing that sp(x,y) must be a polynomial of degree kx in x and   c
        c ky in y. the corresponding weights of these additional rows are set  c
        c to 1./p.  iteratively we than have to determine the value of p       c
        c such that f(p)=sum((w[i]*(z[i]-sp(x[i],y[i])))**2) be = s.           c
        c we already know that the least-squares polynomial corresponds to     c
        c p=0  and that the least-squares spline corresponds to p=infinity.    c
        c the iteration process which is proposed here makes use of rational   c
        c interpolation. since f(p) is a convex and strictly decreasing        c
        c function of p, it can be approximated by a rational function r(p)=   c
        c (u*p+v)/(p+w). three values of p(p1,p2,p3) with corresponding values c
        c of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the c
        c new value of p such that r(p)=s. convergence is guaranteed by taking c
        c f1 > 0 and f3 < 0.                                                   c
        cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc =#
    @debug "Entering Part 2. smoothing spline until fp($fp) -> s($s)"        
    kx2 = kx1 + 1
        # test whether there are interior knots in the x-direction.
    if (nk1x == kx1) @goto  L440 end
        # evaluate the discotinuity jumps of the kx-th order derivative of
        # the b-splines at the knots spline.tx[l),l=kx+2,...,nx-kx-1.
    fpdisc!(spline.tx, nx, kx2, wrk1.bx)
    @label L440
    ky2 = ky1 + 1
        # test whether there are interior knots in the y-direction.
    if (nk1y == ky1) @goto  L450 end
        # evaluate the discontinuity jumps of the ky-th order derivative of
        # the b-splines at the knots spline.ty[l),l=ky+2,...,ny-ky-1.
    fpdisc!(spline.ty, ny, ky2, wrk1.by)
        # initial value for p.
    @label L450  
    p1 = 0.
    f1 = wrk1.fp0 - s
    p3 = -one
    f3 = fpms
    p = 0.
    @debug "smoothing initial values: p1:$p1 f1:$f1 p3:$p3 f3:$f3"
    for i in 1:ncof# do 460 i=1,ncof
        p = p + wrk1.a[i, 1]
        # @label L460  #continue
    end
    @debug "p=$p"
    rn = ncof
    @assert p != 0.0
    p = rn / p
    @assert !isnothing(p) 
    # find the bandwidth of the extended observation matrix.
    iband3 = kx1 * nk1y
    iband4 = iband3 + 1
    ich1 = 0
    ich3 = 0
    # iteration process to find the root of f(p)=s.
    for iter in 1:maxit # do 770 iter=1,maxit
        @debug "Start of iteration. p=$p"
        pinv = one / p
        # store the triangularized observation matrix into q.
        for i in 1:ncof# do 480 i=1,ncof
            wrk1.ff[i] = wrk1.f[i]
            for j in 1:iband# do 470 j=1,iband
                wrk1.q[i, j] = wrk1.a[i, j]
                @label L470   #   continue
            end
            ibb = iband + 1
            for j in ibb:iband4 # do 480 j=ibb,iband4
                wrk1.q[i, j] = 0.
                # @label L480    #continue
            end
        end
        if (nk1y == ky1) @goto  L560 end
                # extend the observation matrix with the rows of a matrix, expressing
                # that for x=cst. sp(x,y) must be a polynomial in y of degree ky.
        for i in ky2:nk1y # do 550 i=ky2,nk1y
            ii = i - ky1
            for j in 1:nk1x # do 550 j=1,nk1x
                # initialize the new row.
                for l in 1:iband # do 490 l=1,iband
                    wrk1.h[l] = 0.
                    @label L490  #      continue
                end
                # fill in the non-zero elements of the row. jrot records the column
                # number of the first non-zero element in the row.
                for l in 1:ky2 # do 500 l=1,ky2
                    wrk1.h[l] = wrk1.by[ii,l] * pinv
                end#   @label L500    #    continue

                zi = 0.
                jrot = (j - 1) * nk1y + ii
                # rotate the new row into triangle by givens transformations without
                # square roots.
                for irot in jrot:ncof# do 540 irot=jrot,ncof
                    piv = wrk1.h[1]
                    i2 = min(iband1, ncof - irot)
                    if (piv == 0.) 
                        if (i2 <= 0) @goto L550 
                        else @goto L520 
                        end # 550,550,520
                    end
                    # calculate the parameters of the givens transformation.
                    v_cos, v_sin, wrk1.q[irot, 1] =  fpgivs(piv, wrk1.q[irot, 1])
                    # apply that givens transformation to the right hand side.
                    zi, wrk1.ff[irot]= fprota(v_cos, v_sin, zi, wrk1.ff[irot])
                    if (i2 == 0) @goto  L550 end
                # apply that givens transformation to the left hand side.
                    for l in 1:i2 # do 510 l=1,i2
                        l1 = l + 1
                        wrk1.h[l1], wrk1.q[irot, l1]=fprota(v_cos, v_sin, wrk1.h[l1], wrk1.q[irot, l1])
                        @label L510      #    continue
                    end
                    @label L520  
                    for l in 1:i2 # do 530 l=1,i2
                        wrk1.h[l] = wrk1.h[l + 1]
                    end #  @label L530          continue
                    wrk1.h[i2 + 1] = 0.
                    @label L540    #    continue
                end 
                @label L550    # continue
            end
        end
        @label L560    
        if (nk1x == kx1) @goto  L640 end
        # extend the observation matrix with the rows of a matrix expressing
        # that for y=cst. sp(x,y) must be a polynomial in x of degree kx.
        for i in kx2:nk1x # do 630 i=kx2,nk1x
            ii = i - kx1
            for j in 1:nk1y # do 630 j=1,nk1y
                # initialize the new row
                for l in 1:iband4 # do 570 l=1,iband4
                    wrk1.h[l] = 0.
                end# @label L570        continue
                # fill in the non-zero elements of the row. jrot records the column
                # numb        er of the first non-zero element in the row.
                j1 = 1
                for l in 1:kx2 # do 580 l=1,kx2
                    wrk1.h[j1] = wrk1.bx[ii,l] * pinv
                    j1 = j1 + nk1y
                end#  @label L580        continue
                zi = 0.
                jrot = (i - kx2) * nk1y + j
                # rotate the new row into triangle by givens transformations .
                for irot in jrot:ncof# do 620 irot=jrot,ncof
                    piv = wrk1.h[1]
                    i2 = min(iband3, ncof - irot)
                    if (piv == 0.) 
                        if (i2 <= 0) @goto L630 else @goto L600 end
                    end
                    # calculate the parameters of the givens transformation.
                    v_cos, v_sin, wrk1.q[irot, 1]= fpgivs(piv, wrk1.q[irot, 1])
                    # apply that givens transformation to the right hand side.
                    zi, wrk1.ff[irot] =fprota(v_cos, v_sin, zi, wrk1.ff[irot])
                    if (i2 == 0) @goto  L630 end
                    # apply that givens transformation to the left hand side.
                    for l in 1:i2 # do 590 l=1,i2
                        l1 = l + 1
                        wrk1.h[l1], wrk1.q[irot, l1] = fprota(v_cos, v_sin, wrk1.h[l1], wrk1.q[irot, l1])
                    end #  @label L590          continue
                    @label L600 
                    for l in 1:i2 # do 610 l=1,i2
                        wrk1.h[l] = wrk1.h[l + 1]
                    end#  @label L610          continue
                    wrk1.h[i2 + 1] = 0.
                    @label L620  #       continue
                end
                @label L630  #  continue
            end  
        end
        # find dmax, the maximum value for the diagonal elements in the
        # reduced triangle.
        @label L640    
        dmax = 0.
        for i in 1:ncof # do 650 i=1,ncof
            if (wrk1.q[i, 1] <= dmax) @goto  L650 end
            dmax = wrk1.q[i, 1]
            @label L650 #   continue
        end# 
        # check whether the matrix is rank deficient.
        sigma = eps * dmax
        for i in 1:ncof # do 660 i=1,ncof
            if (wrk1.q[i, 1] <= sigma) @goto  L670 end
        end  # @label L660 #   continue
        # backward substitution in case of full rank.
        @debug "calling fpback" 
        fpback!(wrk1.q, wrk1.ff, ncof, iband4,  spline.c)
        @debug "after fpback c:" spline.c[1:ncof]
        rank = ncof
        @goto  L675
        # in case of rank deficiency, find the minimum norm solution.
        @label L670    
        lwest = ncof * iband4 + ncof + iband4
        #if (lwrk < lwest) @goto  L780 end
        lf = 1
        lh = lf + ncof
        la = lh + iband4
        #        sq, rank = fprank(wrk1.q, wrk1.ff, ncof, iband,  sigma, spline.c,  wrk2)#wrk(la), wrk(lf), wrk(lh))
        @debug "calling fprank" 
        sq, rank = fprank(wrk1.q, wrk1.ff, ncof, iband4,   sigma, spline.c,  wrk2)#wrk[la],   wrk[lf], wrk[lh])
        @debug "after fprank c:" spline.c[1:ncof]
        @label L675
        for i in 1:ncof #    do 680 i=1,ncof
            wrk1.q[i, 1] = wrk1.q[i, 1] / dmax
        end#  @label L680    continue
        @debug "compute f(p)."
        fp = 0.
        for num in 1:nreg # do 720 num = 1,nreg
            num1 = num - 1
            lx = Int(trunc(num1 / nyy))
            ly = num1 - lx * nyy
            jrot = lx * nk1y + ly
            in = wrk1.index[num]
            @label L690
            if (in == 0) @goto  L720 end
            store = 0.
            i1 = jrot
            for i in 1:kx1 # do 710 i=1,kx1
                hxi = wrk1.spx[in,i]
                j1 = i1
                for j in 1:ky1 # do 700 j=1,ky1
                    j1 = j1 + 1
                    @debug "store = $store + $hxi * $(wrk1.spy[in,j]) * $(spline.c[j1])"
                    store = store + hxi * wrk1.spy[in,j] * spline.c[j1]
                end#  @label L700        continue
                i1 = i1 + nk1y
            end#  @label L710      continue
            @debug "fp calc num=$num, fp= $fp + ($(spline.w[in]) *($(spline.z[in]) - $store))^2"
            fp = fp + (spline.w[in] * (spline.z[in] - store))^2
            
            in = wrk1.nummer[in]
            @goto  L690

            @label L720  #  continue
        end
        @debug "fp now $fp"
        # test whether the approximation sp(x,y) is an acceptable solution.
        fpms = fp - s
        @debug "testing is fp($fp) -> s($s) : delta $(abs(fpms)) <= $acc"
        if (abs(fpms) <= acc) @goto  L820 end
        # test whether the maximum allowable number of iterations has been
        # reached.
        if (iter == maxit) @goto  L795 end
        # carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        
        @info ("smothing p1:$p1  f1:$f1, p:$p, f2=$f2, p3:$p3, f3=$f3")
        if (ich3 != 0) @goto  L740 end
        if ((f2 - f3) > acc) @goto  L730 end
        @debug "our initial choice of p is too large."
        
        p3 = p2
        f3 = f2
        p = p * con4
       
        if (p <= p1) 
            @debug " p($p) < p1($p1)"
            p = p1 * con9 + p2 * con1  
            @debug "p now $p"
        end
        @goto  L770
        @label L730
        if (f2 < 0.) ich3 = 1 end
        @label L740    
        if (ich1 != 0) @goto  L760 end
        if ((f1 - f2) > acc) @goto  L750 end
        @debug "our initial choice of p is too small"
        p1 = p2
        f1 = f2
        p = p / con4
        @debug "p=$p"
        if (p3 < 0.) @goto  L770 end
        if (p >= p3)
            @debug " p($p) >= p3($p3)"
             p = p2 * con1 + p3 * con9
            @debug "p now $p"
        end
       
        @goto  L770
        @label L750    
        if (f2 > 0.) ich1 = 1 end
        # test whether the iteration process proceeds as theoretically
        # expected.
        @label L760
        if ((f2 >= f1) || (f2 <= f3)) 
            @error "impossible result ((f2($f2) >= f1($f1)) || (f2($f2) <= f3($f3))"
            @goto  L800 end
        # find the new value of p.

        @info ("calling fprati $p,$p1,$f1,$p2,$f2,$p3,$f3")
        p, p1, p2, p3, f1, f2, f3 = fprati(p1, f1, p2, f2, p3, f3)
        @info ("return fprati $p,$p1,$f1,$p2,$f2,$p3,$f3")
        @label L770 # continue
        @debug " New p=$p"
    end
        # error codes and messages.
    @label L780  
    @debug "L780"
    ier = lwest
    @goto  L830
    @label L785  
    @debug "L785"
    ier = 5
    @goto  L830
    @label L790  
    @debug "L790"
    ier = 4
    @goto  L830
    @label L795  
    @debug "L795"
    ier = 3
    @goto  L830
    @label L800  
    @debug "L800"
    ier = 2
    @goto  L830
    @label L810  
    @debug "L810"
    ier = 1
    @goto  L830
    @label L815  
    @debug "L830"
    ier = -1
    fp = 0.
    @label L820 
    @debug "L820"

    if (ncof != rank) ier = -rank end
        # test whether x and y are in the original order.
    @label L830  
    @debug "Label 830. ichang = $ichang"
    if (ichang < 0) @goto  L930 end
        # if not, interchange x and y once more.
    l1 = 1
    for i in 1:nk1x# do 840 i=1,nk1x
        l2 = i
        for i in 1:nk1y # do 840 j=1,nk1y
            wrk1.f[l2] = spline.c[l1]
            l1 = l1 + 1
            l2 = l2 + nk1x
        end #  @label L840  continue
    end
    for i in 1:ncof # do 850 i=1,ncof
        spline.c[i] = wrk1.f[i]
    end #      @label L850  continue
    for i in 1:m # do 860 i=1,m
        store = spline.x[i]
        spline.x[i] = spline.y[i]
        spline.y[i] = store
    end #  @label L860  continue
    n = min(nx, ny)
    for i in 1:n # do 870 i=1,n
        store = spline.tx[i]
        spline.tx[i] = spline.ty[i]
        spline.ty[i] = store
    end#  @label L870  continue
    n1 = n + 1
    if (nx - ny) < 0 @goto L880
    elseif nx == ny @goto  L920
        else @goto  L900 end
    @label L880  
    for i in n1:ny # do 890 i=n1,ny
        spline.tx[i] = spline.ty[i]
    end#  @label L890  continue
    @goto  L920
    @label L900 
    for i in n1:nx # do 910 i=n1,nx
        spline.ty[i] = spline.tx[i]
    end #  @label L910  continue
    @label L920  
    l = nx
    nx = ny
    ny = l
    @label L930  
    if (iopt < 0) @goto  L940 end
    spline.nx = nx
    spline.ny = ny
    @label L940  
    spline.fp = fp
    @debug "returning with ier=$ier and fp = $(spline.fp)"
    return ier
end

