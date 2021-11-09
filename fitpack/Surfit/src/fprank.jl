"""   
subroutine fprank finds the minimum norm solution of a least-
squares problem in case of rank deficiency.

args:
- a : array, which contains the non-zero elements of the observation
         matrix after triangularization by givens transformations.
- f : array, which contains the transformed right hand side.
- n : integer,wich contains the dimension of a.
- m : integer, which denotes the bandwidth of a.
- tol : real value, giving a threshold to determine the rank of a.
- c : vector (modified), which contains the minimum norm solution.
- wrk: :lwork2, scratchpad
returns:
- sq : real value, giving the contribution of reducing the rank
        to the sum of squared residuals.
- rank : integer, which contains the rank of matrix a.
"""
function fprank!(a, f, n, m,  tol, c::Vector,  wrk::lwork2)# aa, ff, h)
   
    @debug "fprank"
    i=0; ii=0; ij=0; i1=0; i2=0; j=0; jj=0; j1=0; j2=0; j3=0; k=0; kk=0; m1=0; nl=0
    m1 = m - 1
    v_cos = 0.0; v_sin = 0.0; fac=0.0 ; piv = 0.0; yi=0.0
    #  the rank deficiency nl is considered to be the number of sufficient
    #  small diagonal elements of a.
    nl = 0
    sq = 0.
    #@debug "n=$n, m=$m, f:$(f[1:n])"
    for i in 1:n # do 90 i=1,n
        if (a[i,1] > tol) continue end # @goto L90 end
        #  if a sufficient small diagonal element is found, we put it to
        #  zero. the remainder of the row corresponding to that zero diagonal
        #  element is then rotated into triangle by givens rotations .
        #  the rank deficiency is increased by one.
        nl = nl + 1
        if (i == n) continue end # @goto L90 end
        yi = f[i]
        @debug "yi=$yi"
        for j in 1:m1 # do 10 j=1,m1
            wrk.h[j] = a[i,j + 1]
        end # 10    continue
        wrk.h[m] = 0.
        i1 = i + 1
        for ii in i1:n # do 60 ii=i1,n            
            i2 = min(n - ii, m1)
            piv = wrk.h[1]
            #if (piv == 0.) @goto L30 end
            if piv != 0.
                (v_cos, v_sin, a[ii,1]) = fpgivs(piv, a[ii,1])
                yi, f[ii] = fprota(v_cos, v_sin, yi, f[ii])
                @debug "yi=$yi"
                #@debug yi
                if (i2 == 0) break end# @goto L70 end
                for j in 1:i2 # do 20 j=1,i2
                    j1 = j + 1
                    wrk.h[j1], a[ii,j1] = fprota(v_cos, v_sin, wrk.h[j1], a[ii,j1])
                    wrk.h[j] = wrk.h[j1]
                end # 20      continue
                #@goto L50
            else
                #@label L30
                if (i2 == 0) break end #@goto L70 end
                for j in 1:i2 # do 40 j=1,i2
                    wrk.h[j] = wrk.h[j + 1]
                end # 40      continue
            end # @label L50
            wrk.h[i2 + 1] = 0.         
        end # 60    continue
       
        #  add to the sum of squared residuals the contribution of deleting
        #  the row with small diagonal element.
        #@label L70
        sq = sq + yi^2
        @debug sq, yi
        #@label L90
    end # 90  continue
    #  rank denotes the rank of a.
    rank = n - nl
    #=
    #  let b denote the (rank*n) upper trapezoidal matrix which can be
    #  obtained from the (n*n) upper triangular matrix a by deleting
    #  the rows and interchanging the columns corresponding to a zero
    #  diagonal element. if this matrix is factorized using givens
    #  transformations as  b = (r) (u)  where
    #    r is a (rank*rank) upper triangular matrix,
    #    u is a (rank*n) orthonormal matrix
    #  then the minimal least-squares solution c is given by c = b' v,
    #  where v is the solution of the system  (r) (r)' v = g  and
    #  g denotes the vector obtained from the old right hand side f, by
    #  removing the elements corresponding to a zero diagonal element of a.
    #  initialization.
    =#
    fill!(wrk.aa, 0.)
    #=
    for i in 1:rank # do 100 i=1,rank
        for j in 1:m# do 100 j=1,m
            wrk.aa[i,j] = 0.
        end 
    end 
    =#
    #  form in aa the upper triangular matrix obtained from a by
    #  removing rows and columns with zero diagonal elements. form in ff
    #  the new right hand side by removing the elements of the old right
    #  hand side corresponding to a deleted row.
    ii = 0
    for i in 1:n # do 120 i=1,n
        if (a[i,1] <= tol) continue end #@goto L120 end
        ii = ii + 1
        wrk.ff[ii] = f[i]
        wrk.aa[ii,1] = a[i,1]
        jj = ii
        kk = 1
        j = i
        j1 = min(j - 1, m1)
        if (j1 != 0) #@goto L120
            for k in 1:j1 # do 110 k=1,j1
                j = j - 1
                #if (a[j,1] <= tol) @goto L110 end
                if a[j,1] > tol
                    kk = kk + 1
                    jj = jj - 1
                    wrk.aa[jj,kk] = a[j,k + 1]
                end #@label L110#    continue
            end
            #@label L120 # continue
        end
    end
   #  form successively in h the columns of a with a zero diagonal element.
    ii = 0
    for i in 1:n # do 200 i=1,n
        ii = ii + 1
        if (a[i,1] > tol) continue end#@goto L200 end
        ii = ii - 1
        if (ii == 0) continue end # @goto L200 end
        jj = 1
        j = i
        j1 = min(j - 1, m1)
        for k in 1:j1# do 130 k=1,j1
            j = j - 1
            #if (a[j,1] <= tol) @goto L130 end
            if a[j,1] > tol
                wrk.h[jj] = a[j,k + 1]
                jj = jj + 1          
            end #    @label L130  #  continue
        end
        for kk in jj:m# do 140 kk=jj,m
            wrk.h[kk] = 0.
        end # 140    continue
        #  rotate this column into wrk.aa by givens transformations.
        jj = ii
        #@debug "before jj=$jj"
        for i1 in 1:ii # do 190 i1=1,ii
            j1 = min(jj - 1, m1)
            piv = wrk.h[1]
            #if (piv != 0.) @goto L160 end
            if piv == 0
                if (j1 == 0) break end # @goto L200 end
                for j2 in 1:j1 # do 150 j2=1,j1
                    j3 = j2 + 1
                    wrk.h[j2] = wrk.h[j3]
                end # 150      continue
                #@goto L180
            else # @label L160 
                v_cos, v_sin, wrk.aa[jj,1] = fpgivs(piv, wrk.aa[jj,1])
                if (j1 == 0) break end #@goto L200 end
                @debug "in loop jj=$jj"
                kk = jj
                for j2 in 1:j1 # do 170 j2=1,j1
                    j3 = j2 + 1
                    kk = kk - 1
                    wrk.h[j3], wrk.aa[kk,j3] = fprota(v_cos, v_sin, wrk.h[j3], wrk.aa[kk,j3])
                    wrk.h[j2] = wrk.h[j3]
                end # 170      continue
            end #@label L180      
            jj = jj - 1
            #@debug "in loop jj=$jj"
            wrk.h[j3] = 0.
        end # 190    continue
        #@label L200
    end # 200  continue
   #  solve the system (aa) (f1) = ff
    wrk.ff[rank] = wrk.ff[rank] / wrk.aa[rank,1]
    i = rank - 1
    if i != 0 #if (i == 0) @goto L230 end
        for j in 2:rank # do 220 j=2,rank
            store = wrk.ff[i]
            i1 = min(j - 1, m1)
            k = i
            for ii in 1:i1# do 210 ii=1,i1
                k = k + 1
                stor1 = wrk.ff[k]
                stor2 = wrk.aa[i,ii + 1]
                store = store - stor1 * stor2
            end # 210    continue
            stor1 = wrk.aa[i,1]
            wrk.ff[i] = store / stor1
            i = i - 1
        end # 220  continue
        #  solve the system  (aa)' (f2) = f1
    end #@label L230  
    wrk.ff[1] = wrk.ff[1] / wrk.aa[1,1]
    if rank != 1# if (rank == 1) @goto L260 end
        for j in 2:rank # do 250 j=2,rank
            store = wrk.ff[j]
            i1 = min(j - 1, m1)
            k = j
            for ii in 1:i1 # do 240 ii=1,i1
                k = k - 1
                stor1 = wrk.ff[k]
                stor2 = wrk.aa[k,ii + 1]
                store = store - stor1 * stor2
            end # 240    continue
            stor1 = wrk.aa[j,1]
            wrk.ff[j] = store / stor1
        end # 250  continue
        #  premultiply f2 by the transpoze of a.
    end # @label L260
    k = 0
    for i in 1:n# do 280 i=1,n
        store = 0.
        if (a[i,1] > tol) k = k + 1 end
        j1 = min(i, m)
        kk = k
        ij = i + 1
        for j in 1:j1 # do 270 j=1,j1
            ij = ij - 1
            if (a[ij,1] > tol) #  if (a[ij,1] <= tol) 
                #@debug "a[$ij, 1] ($(a[ij,1]) <= $tol"
                #@goto L270 end
                stor1 = a[ij,j]
                stor2 = wrk.ff[kk]
                #@debug "stor1: $stor1, stor2: $stor2"
                store = store + stor1 * stor2
                kk = kk - 1
                #@label L270
            end
        end # 270    continue
        #@debug "fprank: c[$i]=$store"
        c[i] = store
        # @debug "c[$i] = $store"
    end # 280  continue
   #  add to the sum of squared residuals the contribution of putting
   #  to zero the small diagonal elements of matrix (a).
    stor3 = 0.
    for i in 1:n # do 310 i=1,n
        if (a[i,1] > tol) continue end # @goto L310 end
        store = f[i]
        i1 = min(n - i, m1)
        if i1 != 0 #  if (i1 == 0) @goto L300 end
            for j in 1:i1 # do 290 j=1,i1
                ij = i + j
                stor1 = c[ij]
                stor2 = a[i,j + 1]
                store = store - stor1 * stor2
            end # 290    continue
        end #     @label L300
        fac = a[i,1] * c[i]
        stor1 = a[i,1]
        stor2 = c[i]
        stor1 = stor1 * stor2
        stor3 = stor3 + stor1 * (stor1 - store - store)
        # @label L310
    end # 310  continue
    fac = stor3
    sq = sq + fac
    return sq, rank
end

