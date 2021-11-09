module Types

export lwork1, lwork2, Spline

print("types.jl")
import Base: show
using Printf


mutable struct lwork1
    fp0::Float64
    fpint
    coord
    f
    ff
    a
    q
    bx
    by
    spx
    spy
    h
    index
    nummer
end

mutable struct lwork2
    aa::Matrix
    ff::Vector
    h::Vector
end

mutable struct Spline
    x::Vector
    y::Vector
    z::Vector
    w::Vector
    xb::Float64
    xe::Float64
    yb::Float64
    ye::Float64
    kx::Int
    ky::Int
    
    nx::Int
    tx::Vector
    ny::Int
    ty::Vector
    c::Vector
    fp::Float64
    errMsg::String
end


function show(io::IO, spl::Spline)
    io_compact = IOContext(io, :compact => true)
    xknots = spl.tx[1:spl.nx]
    yknots = spl.ty[1:spl.ny]
    println(io_compact, """Spline (xknots=$(_reallycompact(xknots)),
            yknots=$(_reallycompact(yknots)),
            kx=$(spl.kx), ky=$(spl.ky),  residual=$(spl.fp))""")
    nr_y_lines = spl.ny - spl.ky - 1
    nr_x_cols = spl.nx - spl.kx - 1
    for row in 1:nr_y_lines
        offset = (row-1)* nr_x_cols + 1
        println(io_compact, "$(_reallycompact(spl.c[offset:offset + nr_x_cols]))" )
    end
    println(io_compact, """  fp=$(spl.fp)
     errMsg = '$(spl.errMsg)'
    """)
end

function _reallycompact(a::Vector;min_elems=9, format="%.3f")
    io = IOBuffer()
    #format="%.3f"
    io_compact = IOContext(io, :compact => true)
    if length(a) <= min_elems
        write(io, "[")
        for elem in a
            @printf(io_compact,  " % .3f" , elem)
        end
        write(io, "]")
    else
        half = (min_elems+1)÷2
        write(io, "[")
        for i in 1:half
            @printf(io_compact, "% .3f ", a[i])
            write(io, ",")
        end
 
        write(io, " \u2026 ")
        for i in -half:0
            col_nr = length(a)+i
            write(io, ",")
            @printf(io_compact, " % .3f", a[col_nr]) 
        end
        write(io, "]")
        write(io, " ($(length(a)) elements)")
    end
    seekstart(io)
    return read(io, String)
end

function showMatrix(m::Matrix, width::Int, height::Int; min_rows=9, min_elems=9, format= "%.3f", io=stdout)
    io_compact = IOContext(io, :compact => true)
    if height <= min_rows
        for r in 1:height
            println(io_compact, """$(@sprintf("%2u",r)) $(_reallycompact(m[r,1:width], min_elems=min_elems,format=format))""")
        end
    else
        half = (min_rows + 1) ÷ 2
        for r in 1:half
            println(io_compact, """$(@sprintf("%2u",r)) $(_reallycompact(m[r,1:width], min_elems=min_elems,format=format))""")
        end
        println(io_compact, """[...]""")
        for r in -half:0
            row_nr = height + r
            println(io_compact, """$(@sprintf("%2u",row_nr)) $(_reallycompact(m[row_nr,1:width]))""")
        end
    end
end

#=
function show(spline::Spline2D; io=stdout)
    io_compact = IOContext(io, :compact => true)
    
    println(io, "Spline 2D kx:$(spline.kx), ky $(spline.ky), fp: $(spline.fp), fp0: $(spline.fp0), ier: $(spline.ier) ")
    println(io_compact, """tx=$(_reallycompact(spline.tx))""")
    println(io_compact, """ty=$(_reallycompact(spline.ty))""")
    cols = length(spline.tx)-spline.kx - 1
    rows =  length(spline.ty)-spline.ky - 1
    # (rows, cols) = size(spline.c)
    cc = reshape(spline.c[1:rows*cols], rows, cols)
    if rows < 10
        for r in 1:rows
            println(io_compact, """$r $(_reallycompact(cc[r,:]))""")
        end
    else
        for r in 1:4
            println(io_compact, """$r $(_reallycompact(cc[r,:]))""")
        end
        for r in -4:0
            println(io_compact, """$r $(_reallycompact(cc[end+r,:]))""")
        end
    end
end
=#
end