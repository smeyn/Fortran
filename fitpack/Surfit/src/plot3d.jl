using Plots


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

x,y,z, delta  =readData("./test/data/dasurf.txt")   

plot(scatter3d(x,y,z))