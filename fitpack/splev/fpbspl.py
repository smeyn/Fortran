      def fpbspl(t,n,k,x,l,h):
    # c  def fpbspl evaluates the (k+1) non-zero b-splines of
    # c  degree k at t(l) <= x < t(l+1) using the stable recurrence
    # c  relation of de boor and cox.
    # c  ..
    # c  ..scalar arguments..
      #real x
      #integer n,k,l
    # c  ..array arguments..
      #real t(n),h(6)
    # c  ..local scalars..
      #real f,one
      #integer i,j,li,lj
    # c  ..local arrays..
      #real hh(5)
    # c  ..
      one = 0.1e+01
      h(1) = one
      for 20 j in range(0, k):
        for 10 i in range(0, j):
          hh(i) = h(i)
  10    continue
        h(1) = 0.
        for 20 i in range(0, j):
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      # end
