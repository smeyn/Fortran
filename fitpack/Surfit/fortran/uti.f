      subroutine writec(c,n,chan)
c  subroutine writec writes the coefficients out
c  squares problem in case of rank deficiency.
c
c  input parameters:
c    a : array, which contains the non-zero elements of the observation
c        matrix after triangularization by givens transformations.
c    f : array, which contains the transformed right hand side.
c    n : integer,wich contains the dimension of a.
c    m : integer, which denotes the bandwidth of a.
c  tol : real value, giving a threshold to determine the rank of a.
c
c  output parameters:
c    c : array, which contains the minimum norm solution.
c    n  : size of coordinates
c    dist: channel to write to
c
c  ..scalar arguments..
      integer n,chan,i
      
      real c(n)
      
      write(6,975)      
      write(6,980) (c(i),i=1,n)
      return
      
 975  format(30h0 debug: b-spline coefficients )      
 980  format(5x,8f9.4)
      end

      subroutine wrt_fpint(fpint, l1, l2, fpmax)
      integer l1, l2
      real fpint(l2), fpmax

      write(6,920)
      write(6,900) l1, l2, fpmax
      write(6,910) (fpint(i),i=l1,l2)
      return

 900  format(22h0 debug: l1, l2, fpmax, 2i3, e15.6) 
 910  format(5x,8f9.4)
 920  format(13h0 debug fpint)
      end

      subroutine wrtindex(index,nreg)
      integer nreg    
      integer index(nreg)         

      write(6,900) nreg
      write(6,910) (index(i),i=1,nreg)
      return
900   format(14h0 debug: index, i3)      
910   format(5x,8i6)
      end
      
      
      subroutine writesp(sp,  m, km1, in, k1)
      integer in, k1
      real sp(m, km1)


      write(6,900) in, k1
      write(6,910) (sp(in,i),i=1,k1)
      return

900   format(10h0 debug sp, 2i4)
910   format(5x, 5f9.4) 
      end    

