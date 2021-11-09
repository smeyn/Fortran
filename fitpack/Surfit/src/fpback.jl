"""
calculates the solution of the system of equations a*c = z 
with a = n x n upper triangular matrix of bandwidth k.

Args:
- a matrix containing the upper triangular matrix
- z vector 
- n int: actual size of matrix data in a
- c vector : modified
- k : bandwidth
"""
function fpback!(a, z, n, k, c)
    #@assert length(c) >= n  "fpback c too small. $(length(c)) n=$n"
    #@assert length(z) >= n  "fpback z too small.  $(length(z)) n=$n"
    #@assert size(a)[1] >= n  "fpback a too small a: $(size(a))  n=$n"

    k1 = k - 1
    c[n] = z[n] / a[n,1]
    i = n - 1
    if i != 0 
        for j in 2:n 
            store = z[i]
            i1 = k1
            if (j <= k1) i1 = j - 1 end
            m = i
            for l in 1:i1  
                m = m + 1
                store = store - c[m] * a[i,l + 1]
            end  
            c[i] = store / a[i,1]
            i = i - 1
        end  
    end 
    #@debug "fpback a($n,$n)" a[1:n, 1:n]
    #@debug "fpback z($n)" z[1:n]
    #@debug "fpback c"  c
    return
end
