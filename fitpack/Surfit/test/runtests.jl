
using Surfit
#include("../src/Bispev.jl")
using Surfit.Bispev
using Test

#println(names(Surfit;imported=true, all=true))
#println(names(Bispev;imported=true, all=true))

include("logger.jl")

#=
include("fortranresults.jl")

@testset "fortranresults" begin
    LogLevel(Logging.Debug)

    sections = readResults("fortran/surfout.txt")
    @test length(sections) == 6
end
=#

"""test function for surface fitting"""
function zzz(x, y)
    r = x*x + y*y
    z = cos(r) * exp(-r)
    return z
end

function showIER(ier)
    showIER(stdout, ier)
end

function showIER(io, ier)
    print(io, "IER=$ier: ")
    if ier == 0  println(io, "normal return. the spline returned has a residual sum of squares fp") 
    elseif ier == -1 println(io, "normal return. the spline returned is an interpolating spline (fp=0).")
    elseif ier == -2 println(io, "normal return. the spline returned is the weighted least-squares polynomial of degrees kx and ky")
    elseif ier < 2 println(io, "warning. the coefficients of the spline returned have been computed as the minimal norm")
    elseif ier == 1 println(io, "error. the required storage space exceeds the available storage space")
    elseif ier == 2 println(io, "error. a theoretically impossible result was found ")
    elseif ier== 3 println(io, "error. the maximal number of iterations was reached")
    elseif ier == 4 println(io, "error. no more knots can be added")
    elseif ier == 5 println(io, "error. no more knots can be added")
    elseif ier == 10 println(io, "error.error. on entry, the input data are controlled on validity the following restrictions must be satisfied.")
    end
end
 
@testset "mnsurf" begin
    
    include("mnsurf.jl")
end
 
@testset "Surfit" begin

    iopt = 1 #Surfit.START_MINIMAL_KNOTS
    x = [
        0.0,
        0.0,
        0.0,
        0.0,
        1.1,
        1.1,
        1.1,
        1.1,
        2.2,
        2.2,
        2.2,
        2.2,
        3.1,
        3.1,
        3.1,
        3.1,
        2.7, 
        1.5
    ]
    y = [
        0.0,
        0.1,
        2.0,
        3.0,
        0.0,
        0.1,
        2.0,
        3.0,
        0.0,
        0.1,
        2.0,
        3.0,
        0.0,
        0.1,
        2.0,
        3.0,
        1.1,
        1.9
    ]
    z = [zzz(xx,yy) for (xx,yy) in zip(x,y)] 
    w = ones(Float64, length(z))
    m = length(x)
    #tx = zeros(Float64,10)
    #ty = zeros(Float64,10)
    kx = 3
    ky = 3
    s = 1000.0
    eps = 0.0000000001
    println("calling surfit")
    nxest= 2 * (kx + 1)
    nyest = 2 * (ky + 1)

    # calculate lwrk1
    u = nxest-kx-1 
    v = nyest-ky-1 
    km = max(kx,ky)+1
    ne = max(nxest,nyest)
    bx = kx*v+ky+1
    by = ky*u+kx+1
    if(bx <= by) 
        b1 = bx 
        b2 = b1+v-ky
    end 
    if(bx > by) 
        b1 = by, 
        b2 = b1+u-kx
    end
    lwrk1 = u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+2
    
   
    #wrk1 = zeros(Float64, lwrk1)

    # calcualte lwrk2
    lwrk2 = u*v*(b2+1)+b2
    wrk2 = zeros(Float64, lwrk2)

    # iwwrk
    kwrk = m + (nxest-2*kx-1)*(nyest-2*ky-1) + 1

    iwrk = zeros(Int, kwrk)
    nx = 0
    ny = 0
    tx = zeros(Float64, 20)
    ty = zeros(Float64, 20)
    fp = 0.0
    c = zeros(Float64, (length(tx)-kx-1)*(length(ty)-ky-1))
    nmax = min(length(tx), length(ty))

  
    spline = Surfit.Spline(
        w,
        0.0, 3.1, 0.0, 3.0,
        kx, ky, nx, tx, ny, ty, c, fp,""  
    )
    
    #spline = Surfit.prepareSpline(x, y, z, w, 3, 3, 20, 20)
    
    basicLogger = MinLevelLogger(FileLogger("debug.log"), Logging.Info)
    with_logger(basicLogger) do
        result =  Surfit.surfit(iopt, x,y,z,spline,  s, eps)
        showIER(result)
        if result > 0
            @error spline.errMsg
        end
        @test result <= 0
    end
    @info "spline results"
    show(stdout, spline)

    @info "Validate"

    x_test = [ 0.0, 1.1, 2.2, 3.1]
    y_test = [  0.0,  0.1, 2.0, 3.0, 3.2]


    z_result = zeros(length(x_test), length(y_test))
    ier, z1 = Bispev.bispev(spline,
        x_test, y_test, z_result) 
        
    z_expected = zeros(length(x_test), length(y_test))
    for r in 1:length(x_test)
        for c in 1:length(y_test)
            z_expected[r,c] = zzz(x_test[r], y_test[c])
        end
    end
    
    @info "z_result" z_result
    @info "z_expected" z_expected

    @test isapprox(z_result, z_expected; rtol=1e-2)
        
end
