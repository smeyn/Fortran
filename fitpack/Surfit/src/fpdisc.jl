function fpdisc(t, n, k2, b, nest)
   #  subroutine fpdisc calculates the discontinuity jumps of the kth
   #  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
   #  ..scalar arguments..
   #   integer n,k2,nest
   #  ..array arguments..
   #   real t(n),b(nest,k2)
   #  ..local scalars..
   #   real an,fac,prod
   #   integer i,ik,j,jk,k,k1,l,lj,lk,lmk,lp,nk1,nrint
   #  ..local array..
   #   real h(12)
   #  ..
   h=zeros(12)
    k1 = k2 - 1
    k = k1 - 1
    nk1 = n - k1
    nrint = nk1 - k
    an = nrint
    fac = an / (t[nk1 + 1] - t[k1])
    for l in k2:nk1 # do 40 l=k2,nk1
        lmk = l - k1
        for j in 1:k1 # do 10 j=1,k1
            ik = j + k1
            lj = l + j
            lk = lj - k2
            h[j] = t[l] - t[lk]
            h[ik] = t[l] - t[lj]
        end# 10    continue
        lp = lmk
        for j in 1:k2 # do 30 j=1,k2
            jk = j
            prod = h[j]
            for i in 1:k # do 20 i=1,k
                jk = jk + 1
                prod = prod * h[jk] * fac
            end # 20      continue
            lk = lp + k1
            b[lmk,j] = (t[lk] - t[lp]) / prod
            lp = lp + 1
        end # 30    continue
    end # 40  continue
    return
end
