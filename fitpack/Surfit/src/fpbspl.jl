function fpbspl(t, n, k, x, l, h)
   #  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
   #  degree k at t(l) <= x < t(l+1) using the stable recurrence
   #  relation of de boor and cox.
   #  ..
   #  ..scalar arguments..
   #   real x
   #   integer n,k,l
   #  ..array arguments..
   #   real t(n),h(6)
   #  ..local scalars..
   #   real f,one
   #   integer i,j,li,lj
   #  ..local arrays..
   #   real hh(5)
   #  ..
   #@info "fpbspl:t $t[1:n]\n       n:$n, k:$k, x:$x, l:$l\n       h $h" 
   hh = zeros(Float64, 5)
    one = 0.1e+01
    h[1] = one
    for j in 1:k # do 20 j=1,k
        for i in 1:j # do 10 i=1,j
            hh[i] = h[i]
        end # 10    continue
        h[1] = 0.
        for i in 1:j # do 20 i=1,j
            li = l + i
            lj = li - j
            f = hh[i] / (t[li] - t[lj])
            h[i] = h[i] + f * (t[li] - x)
            h[i + 1] = f * (x - t[lj])
        end
    end# 20  continue
    #@info "fpbspl:t $t[1:n]\n       n:$n, k:$k, x:$x, l:$l\n       h $h" 
    return
end
