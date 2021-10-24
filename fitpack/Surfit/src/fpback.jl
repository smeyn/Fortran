function fpback(a, z, n, k, c, nest)
# c  subroutine fpback calculates the solution of the system of
# c  equations a*c = z with a a n x n upper triangular matrix
# c  of bandwidth k.
# c  ..
# c  ..scalar arguments..
#      integer n,k,nest
# c  ..array arguments..
#      real a(nest,k),z(n),c(n)
# c  ..local scalars..
#      real store
#      integer i,i1,j,k1,l,m
# c  ..
@assert length(c) >= n  "fpback c too small. $(length(c)) n=$n"
@assert length(z) >= n  "fpback z too small.  $(length(z)) n=$n"
@assert size(a)[1] >= n  "fpback a too small a: $(size(a))  n=$n"

    k1 = k - 1
    c[n] = z[n] / a[n,1]
    i = n - 1
    if (i == 0) @goto L30 end
    for j in 2:n # do 20 j=2,n
        store = z[i]
        i1 = k1
        if (j <= k1) i1 = j - 1 end
        m = i
        for l in 1:i1 # do 10 l=1,i1
            m = m + 1
            store = store - c[m] * a[i,l + 1]
        end # 10    continue
        c[i] = store / a[i,1]
        i = i - 1
    end # 20  continue
    @label L30  
    #@debug "fpback a($n,$n)" a[1:n, 1:n]
    #@debug "fpback z($n)" z[1:n]
    #@debug "fpback c"  c
    return
end
