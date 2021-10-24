function fporde(x, y, m, kx, ky, tx, nx, ty, ny, nummer, index, nreg)
   #  subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m
   #  according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong
   #  to. for each panel a stack is constructed  containing the numbers
   #  of data points lying inside; index(j),j=1,2,...,nreg points to the
   #  first data point in the jth panel while nummer(i),i=1,2,...,m gives
   #  the number of the next data point in the panel.
   #  ..
   #  ..scalar arguments..
   #   integer m,kx,ky,nx,ny,nreg
   #  ..array arguments..
   #   real x(m),y(m),tx(nx),ty(ny)
   #   integer nummer(m),index(nreg)
   #  ..local scalars..
   #   real xi,yi
   #   integer i,im,k,kx1,ky1,k1,l,l1,nk1x,nk1y,num,nyy
   #  ..
    kx1 = kx + 1
    ky1 = ky + 1
    nk1x = nx - kx1
    nk1y = ny - ky1
    nyy = nk1y - ky
    for i in 1:nreg # do 10 i=1,nreg
        index[i] = 0
    end # 10  continue
    for im in 1:m # do 60 im=1,m
        xi = x[im]
        yi = y[im]
        l = kx1
        l1 = l + 1
        @label L20 
        if (xi < tx[l1] || l == nk1x) @goto L30 end
        l = l1
        l1 = l + 1
        @goto L20
        @label L30
        k = ky1
        k1 = k + 1
        @label L40
        if (yi < ty[k1] || k == nk1y) @goto L50 end
        k = k1
        k1 = k + 1
        @goto L40
        @label L50    
        num = (l - kx1) * nyy + k - ky
        nummer[im] = index[num]
        index[num] = im
    end # 60  continue
    return
end

