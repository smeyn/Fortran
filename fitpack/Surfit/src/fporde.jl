"""
   fporde sorts the data points (x(i),y(i)),i=1,2,...,m
   according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong
   to. for each panel a stack is constructed  containing the numbers
   of data points lying inside; index(j),j=1,2,...,nreg points to the
   first data point in the jth panel while nummer(i),i=1,2,...,m gives
   the number of the next data point in the panel.

args:
- x: x coordinates
- y: y coordinates
- m : length of x and y
- kx: degree of spline in x direction
- ky : degree of spline in y direction
- tx: knot vector in x direction
- nx: number of knots in tx
- ty: knot vector in y direction
- ny: number of knots in ty
- nummer: modified, contains the number of points in panel
- index: modified, points to first data point in the panel
"""   
function fporde!(x::Vector, y::Vector, m, kx, ky, tx, nx, ty, ny, nummer, index, nreg)

    kx1 = kx + 1
    ky1 = ky + 1
    nk1x = nx - kx1
    nk1y = ny - ky1
    nyy = nk1y - ky
    fill!(index, 0)
    #for i in 1:nreg # do 10 i=1,nreg
    #    index[i] = 0
    #end # 10  continue
    for im in 1:m # do 60 im=1,m
        xi = x[im]
        yi = y[im]
        l = kx1
        l1 = l + 1
        #@label L20 
        while (xi >= tx[l1]) && (l != nk1x) # if (xi < tx[l1] || l == nk1x) @goto L30 end
            l = l1
            l1 = l + 1
        end #   @goto L20
        #@label L30
        k = ky1
        k1 = k + 1
        #@label L40
        while (yi >= ty[k1]) && (k != nk1y)  # if (yi < ty[k1] || k == nk1y) @goto L50 end
            k = k1
            k1 = k + 1
            #@goto L40
        end #@label L50    
        num = (l - kx1) * nyy + k - ky
        nummer[im] = index[num]
        index[num] = im
    end # 60  continue
    return
end

