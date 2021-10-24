# cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# cc                                                                    cc
# cc        mnsurf : surfit test program                                cc3
# cc                                                                    cc
# cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# real x(80),y(80),z(80),w(80),tx(15),ty(15),c(200),wrk1(12000),
# * wrk2(6000),xx(11),yy(11),zz(121)

using Printf

function readInt(io)
    txt = readline(io)
    return parse(Int, txt)
end

function read3Floats(io)
    txts = [s for s in split(readline(io), " ") if length(s) > 0]
    floats = [parse(Float64, t) for t in txts]
    return floats
end

function readFloat(io)
    txt = readline(io)
    return parse(Float64, txt)
end

function readData(filepath)
    open(filepath) do io
      # c  we fetch the number of data points
        m1 = readInt(io)
        m = m1
        println("m: $m")
      #  we fetch the co-ordinate and function values of each data point.
        x = zeros(m)
        y = zeros(m)
        z = zeros(m)
        for i in 1:m
            x[i], y[i], z[i] = read3Floats(io)
         # y[i]=read(io, Float64)
         # z[i]=read(io, Float64)
            println("$i: $(x[i]), $(y[i]), $(z[i])")
        end
        delta = readFloat(io)
        println("delta: $delta")
        return x, y, z, delta
    end
end

function dumpResults(path, spline, s, iopt, ier)
    open(path, "w") do io
        showIER(io, ier)
        println(io, spline.errMsg)
        if (iopt >= 0) @goto L210 end
        println(io, "least-squares spline of degrees $(spline.kx), $(spline.ky)") # write(6,935) kx,ky
        @goto L220
        @label L210
        println(io, "smoothing spline of degrees $(spline.kx), $(spline.ky)")# write(6,940) kx,ky
        println(io, "smoothing factor s=$s")# write(6,945) s
        @label L220
        println(io, "sum squared residuals =$(spline.fp), error flag=$ier")# write(6,950) fp,ier
        println(io, "total number of knots in the x-direction = $(spline.nx)") # write(6,955) nx
        println(io, "position of the knots ") # write(6,960)
        println(io, "  ") # write(6,965) (ty(i),i=1,ny)
        for i in 1:spline.nx 
            @printf(io, " % .3f", spline.tx[i]) 
        end
        println(io, " ") # write(6,965) (tx(i),i=1,nx)
        println(io, "total number of knots in the y-direction =$(spline.nx)") # write(6,970) ny
        println(io, "position of the knots ") # write(6,960)
        println(io, "  ") # write(6,965) (ty(i),i=1,ny)
        for i in 1:spline.ny 
            @printf(io, " % .3f", spline.ty[i]) 
        end
        println(io, " ") # write(6,965) (tx(i),i=1,nx)

        nc = (spline.nx - spline.kx - 1) * (spline.ny - spline.ky - 1)
        println(io, "b-spline coefficients") # write(6,975)
        for i in 1:nc 
            @printf(io, " % .3f", spline.c[i])
        end
        println(io, "") # write(6,980) (c(i),i=1,nc)
   # c  evaluation of the spline approximation.
        zz = zeros(length(xx), length(yy))
        ier = Bispev.bispev(spline, xx, yy, zz)
        println(io, "spline evaluation on a given grid") # write(6,1000)
        print(io, "          ")
        for i in 1:mx 
            @printf(io," % .3f", xx[i]) 
        end
        println(io, "") # write(6,985) (xx(i),i=1,mx)
 
        for j in 1:my # do 230 j=1,my
            @printf(io," % .3f", yy[j])
            print(io, ":  ")
            for i in j:11:121 
                @printf(io," % .3f",  zz[i]) 
            end
            println(io, "") # write(6,995) yy(j),(zz(i),i=j,121,11)
        end
    end

end   
x, y, z, delta = readData("data/dasurf.txt")
m = length(x)
w = zeros(length(x))


# integer iwrk(300)
# integer i,ier,iopt,is,j,kwrk,kx,ky,lwrk1,lwrk2,m,mx,my,nc,
# * nmax,nx,nxest,ny,nyest
# real ai,delta,eps,fp,s,ww,xb,xe,yb,ye
# c  we fetch the number of data points
 

println("m = $m") 
@assert m > 0
#  the weights are set equal to delta**(-1)
ww = 1.0 / delta
for i in 1:m
    w[i] = ww
end

#  we set up the boundaries of the approximation domain.
xb = -2.
xe = 2.
yb = -2.
ye = 2.
# we generate a rectangular grid for evaluating the splines.
mx = 11
my = 11
xx = zeros(mx)
yy = zeros(my)
for i in 1:11
    ai = i - 6
    xx[i] = ai * 0.4
    yy[i] = xx[i]
end
#  we set up the dimension information
nxest = 15
nyest = 15
nmax = 15
kwrk = 300
lwrk1 = 12000
lwrk2 = 6000
#  we choose a value for eps
eps = 0.1e-05
#  main loop for the different spline approximations.
println("Spline with $m elements")
spline = Surfit.prepareSpline(copy(x[1:m]), copy(y[1:m]), copy(z[1:m]), w, 3, 3, nxest, nyest;
   xb=xb, xe=xe, yb=yb, ye=ye)

iopt = 0   
for is in 1:6 # do 300 is=1,6
    if is == 1 @goto L110
    elseif is == 2 @goto L120
   elseif is == 3 @goto L130
   elseif is == 4 @goto L140
   elseif is == 5 @goto L150
   elseif is == 6 @goto L160
   else @goto L300 end
   # @goto L(110,120,130,140,150,160),is
   #  we start computing the least-squares bicubic polynomial (large s)
    @label L110 
    global iopt = 0
    spline.kx = 3
    spline.ky = 3
    s =  900000.
    global expected_x_knots = [-2.000, -2.000, -2.000, -2.000,  2.000,  2.000,  2.000,  2.000]
    global expected_y_knots = [-2.000, -2.000, -2.000, -2.000,  2.000,  2.000,  2.000,  2.000] 
    global expected_fp =0.823754E+04 
    
    @goto L200
   # c  iopt=1 from the second call on.
    @label L120
    global iopt = 1
    s = 200.
    global expected_x_knots = [-2.000, -2.000, -2.000, -2.000, -0.456, -0.010,  2.000,  2.000,  2.000,  2.000]
    global expected_y_knots = [-2.000, -2.000, -2.000, -2.000, -0.177,  0.581,  2.000,  2.000,  2.000,  2.000] 
    global expected_fp =0.200072E+03
    
    @goto L200
   #  a value for s within its confidence interval
    @label L130
    global iopt = 1
    s = m
    global expected_x_knots = [-2.000, -2.000, -2.000, -2.000,  -0.456, -0.010 , 0.893,   2.000,  2.000,  2.000,  2.000]
    global expected_y_knots = [-2.000, -2.000, -2.000, -2.000, -1.275, -0.177,  0.581,  2.000,  2.000,  2.000,  2.000] 
    global expected_fp = 0.799988E+02
    @goto L200
   #  overfitting (s too small)
    @label L140
    global iopt = 1
    s = 20.
    global expected_x_knots = [-2.000, -2.000, -2.000, -2.000,  -0.456, -0.010, 0.469 , 0.893,   2.000,  2.000,  2.000,  2.000]
    global expected_y_knots = [-2.000, -2.000, -2.000, -2.000,  -1.275, -0.177,  0.581, 1.059,   2.000,  2.000,  2.000,  2.000] 
    global expected_fp = 0.200086E+02

    @goto L200
   # c  we change the degrees of the spline
    @label L150 
    global iopt = 0
    spline.kx = 5
    spline.ky = 5
    s = m
    global expected_x_knots = [-2.000, -2.000, -2.000, -2.000,-2.000, -2.000,  0.052,  2.000,  2.000,  2.000,   2.000,  2.000,  2.000]
    global expected_y_knots =[-2.000, -2.000, -2.000, -2.000, -2.000, -2.000,   -0.166 ,  2.000,  2.000,  2.000,   2.000,  2.000,  2.000]
    global expected_fp =  0.800291E+02

    @goto L200
   # c  finally, we also calculate a least-squares spline approximation
   # c  with specified knots.
    @label L160
    s = m
    global iopt = -1
    spline.kx = 3
    spline.ky = 3
    spline.nx = 11
    spline.ny = 11
    spline.tx = zeros(spline.nx)
    spline.ty = zeros(spline.ny)

    j = spline.kx + 2
    for i in 1:3 # do 170 i=1,3
        ai = i - 2
        spline.tx[j] = ai
        spline.ty[j] = ai
        j = j + 1
    end # 170    continue
    global expected_x_knots = [-2.000, -2.000, -2.000, -2.000, -1.000,  0.000,  1.000,   2.000,  2.000,   2.000,  2.000]
    global expected_y_knots = [-2.000, -2.000, -2.000, -2.000, -1.000,  0.000,  1.000,   2.000,  2.000,   2.000,  2.000]
    global expected_fp =  0.241781E+02   
   #  determination of the spline approximation.
    @label L200
    ier = 0
    basicLogger = MinLevelLogger(FileLogger("mnsurf-$is.log"), Logging.Debug)
    with_logger(basicLogger) do
        open("mnsurfout-$is.txt", "w") do io
            ier = Surfit.surfit(iopt, spline,  s,   eps;io=io)
        end
        println("Test case $is")
        showIER(ier)
        @test ier <= 0
        # @test all(spline.tx[1:spline.nx] .≈ expected_x_knots ) atol=0.01
        # @test all(spline.ty[1:spline.ny] .≈ expected_y_knots) atol=0.01
        @test all(isapprox.(spline.tx[1:spline.nx],  expected_x_knots, atol=0.001 ))
        @test all(isapprox.(spline.ty[1:spline.ny],  expected_y_knots, atol=0.001))
        # fp within 1%
        fp_err = (expected_fp-spline.fp)/spline.fp
        @test isapprox(fp_err, 0.0 , atol=0.001) 
       # surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
      #        nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
   #  printing of the fitting results.
        dumpResults("mnsurf_$is.txt", spline, s, iopt, ier)
    end   
   #= 
   if(iopt >= 0) @goto L210 end
   println("least-squares spline of degrees $(spline.kx), $(spline.ky)") # write(6,935) kx,ky
   @goto L220
   @label L210
       println("smoothing spline of degrees $(spline.kx), $(spline.ky)")# write(6,940) kx,ky
   println(" smoothing factor s=$s")# write(6,945) s
   @label L220
   println("sum squared residuals =$(spline.fp), error flag=$ier")#write(6,950) fp,ier
   println("total number of knots in the x-direction =") # write(6,955) nx
   println("position of the knots ") # write(6,960)
   for i in 1:spline.nx println(spline.tx[i]) end
   println("") # write(6,965) (tx(i),i=1,nx)
   println("total number of knots in the y-direction =") # write(6,970) ny
   println("position of the knots ") # write(6,960)
   println("") # write(6,965) (ty(i),i=1,ny)
   for i in 1:spline.ny println(spline.ty[i]) end
   nc = (spline.nx-spline.kx-1)*(spline.ny-spline.ky-1)
   println("b-spline coefficients") # write(6,975)
   for i in 1:nc println(spline.c[i]) end
   println("") # write(6,980) (c(i),i=1,nc)
   #c  evaluation of the spline approximation.
   zz=zeros(length(xx), length(yy))
   ier = Bispev.bispev(spline,xx,yy,zz)
   println("spline evaluation on a given grid") # write(6,1000)
   for i in 1:mx println(xx[i]) end
   println("") # write(6,985) (xx(i),i=1,mx)
   println("") # write(6,990)
   for j in 1:my #do 230 j=1,my
     println("") # write(6,995) yy(j),(zz(i),i=j,121,11)
   end =#
    @label L300  # continue
end
#= 
stop
c  format statements.
900  format(i3)
905  format(1h1,i3,12h data points)
910  format(1h0,2(2x,1hi,5x,4hx(i),6x,4hy(i),6x,4hz(i),6x))
915  format(3f10.4)
920  format(1x,2(i3,3f10.4,5x))
925  format(e20.6)
930  format(1x,40hestimate of standard deviation of z(i) =,e15.6)
935  format(32h0least-squares spline of degrees,2i3)
940  format(28h0smoothing spline of degrees,2i3)
945  format(20h smoothing factor s=,f9.0)
950  format(1x,23hsum squared residuals =,e15.6,5x,11herror flag=,i3)
955  format(1x,42htotal number of knots in the x-direction =,i3)
960  format(1x,22hposition of the knots )
965  format(5x,10f7.3)
970  format(1x,42htotal number of knots in the y-direction =,i3)
975  format(23h0b-spline coefficients )
980  format(5x,8f9.4)
985  format(1h0,1hx,2x,11f7.1)
990  format(3x,1hy)
995  format(1x,f4.1,11f7.3)
1000 format(1h0,33hspline evaluation on a given grid)
 end =#