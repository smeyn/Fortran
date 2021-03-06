module Surfit



include("types.jl") 
using .Types

include("fpsurf.jl")
include("fprati.jl")
include("fpback.jl")
include("fpbspl.jl")
include("fpgivs.jl")
include("fpdisc.jl")
include("fporde.jl")
include("fprank.jl")
include("fprota.jl")
include("Bispev.jl")



function surfit(iopt, spline::Spline,  s, eps;io=nothing)
#= 
       c given the set of data points (x(i),y(i),z(i)) and the set of positive
c numbers w(i),i=1,...,m, subroutine surfit determines a smooth bivar-
c iate spline approximation s(x,y) of degrees kx and ky on the rect-
c angle xb <= x <= xe, yb <= y <= ye.
c if iopt = -1 surfit calculates the weighted least-squares spline
c according to a given set of knots.
c if iopt >= 0 the total numbers nx and ny of these knots and their
c position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
c ally by the routine. the smoothness of s(x,y) is then achieved by
c minimalizing the discontinuity jumps in the derivatives of s(x,y)
c across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
c the amounth of smoothness is determined by the condition that f(p) =
c sum ((w(i)*(z(i)-s(x(i),y(i))))**2) be <= s, with s a given non-neg-
c ative constant, called the smoothing factor.
c the fit is given in the b-spline representation (b-spline coefficients
c c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
c uated by means of subroutine bispev.
c
c calling sequence:
   #      call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
   #     *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
c
c parameters:
   #   iopt  : integer flag. on entry iopt must specify whether a weighted
   #           least-squares spline (iopt=-1) or a smoothing spline (iopt=0
   #           or 1) must be determined.
   #           if iopt=0 the routine will start with an initial set of knots
   #           tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
   #           1,...,ky+1. if iopt=1 the routine will continue with the set
   #           of knots found at the last call of the routine.
   #           attention: a call with iopt=1 must always be immediately pre-
   #                      ceded by another call with iopt=1 or iopt=0.
   #           unchanged on exit.
   #   m     : integer. on entry m must specify the number of data points.
   #           m >= (kx+1)*(ky+1). unchanged on exit.
   #   x     : real array of dimension at least (m).
   #   y     : real array of dimension at least (m).
   #   z     : real array of dimension at least (m).
   #           before entry, x(i),y(i),z(i) must be set to the co-ordinates
   #           of the i-th data point, for i=1,...,m. the order of the data
   #           points is immaterial. unchanged on exit.
   #   w     : real array of dimension at least (m). before entry, w(i) must
   #           be set to the i-th value in the set of weights. the w(i) must
   #           be strictly positive. unchanged on exit.
   #   xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
   #   yb,ye   aries of the rectangular approximation domain.
   #           xb<=x(i)<=xe,yb<=y(i)<=ye,i=1,...,m. unchanged on exit.
   #   kx,ky : integer values. on entry kx and ky must specify the degrees
   #           of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
   #           (kx=ky=3) splines. unchanged on exit.
   #   s     : real. on entry (in case iopt>=0) s must specify the smoothing
   #           factor. s >=0. unchanged on exit.
   #           for advice on the choice of s see further comments
   #   nxest : integer. unchanged on exit.
   #   nyest : integer. unchanged on exit.
   #           on entry, nxest and nyest must specify an upper bound for the
   #           number of knots required in the x- and y-directions respect.
   #           these numbers will also determine the storage space needed by
   #           the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
   #           in most practical situation nxest = kx+1+sqrt(m/2), nyest =
   #           ky+1+sqrt(m/2) will be sufficient. see also further comments.
   #   nmax  : integer. on entry nmax must specify the actual dimension of
   #           the arrays tx and ty. nmax >= nxest, nmax >=nyest.
   #           unchanged on exit.
   #   eps   : real.
   #           on entry, eps must specify a threshold for determining the
   #           effective rank of an over-determined linear system of equat-
   #           ions. 0 < eps < 1.  if the number of decimal digits in the
   #           computer representation of a real number is q, then 10**(-q)
   #           is a suitable value for eps in most practical applications.
   #           unchanged on exit.
   #   nx    : integer.
   #           unless ier=10 (in case iopt >=0), nx will contain the total
   #           number of knots with respect to the x-variable, of the spline
   #           approximation returned. if the computation mode iopt=1 is
   #           used, the value of nx should be left unchanged between sub-
   #           sequent calls.
   #           in case iopt=-1, the value of nx should be specified on entry
   #   tx    : real array of dimension nmax.
   #           on succesful exit, this array will contain the knots of the
   #           spline with respect to the x-variable, i.e. the position of
   #           the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
   #           position of the additional knots tx(1)=...=tx(kx+1)=xb and
   #           tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
   #           if the computation mode iopt=1 is used, the values of tx(1),
   #           ...,tx(nx) should be left unchanged between subsequent calls.
   #           if the computation mode iopt=-1 is used, the values tx(kx+2),
   #           ...tx(nx-kx-1) must be supplied by the user, before entry.
   #           see also the restrictions (ier=10).
   #   ny    : integer.
   #           unless ier=10 (in case iopt >=0), ny will contain the total
   #           number of knots with respect to the y-variable, of the spline
   #           approximation returned. if the computation mode iopt=1 is
   #           used, the value of ny should be left unchanged between sub-
   #           sequent calls.
   #           in case iopt=-1, the value of ny should be specified on entry
   #   ty    : real array of dimension nmax.
   #           on succesful exit, this array will contain the knots of the
   #           spline with respect to the y-variable, i.e. the position of
   #           the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
   #           position of the additional knots ty(1)=...=ty(ky+1)=yb and
   #           ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
   #           if the computation mode iopt=1 is used, the values of ty(1),
   #           ...,ty(ny) should be left unchanged between subsequent calls.
   #           if the computation mode iopt=-1 is used, the values ty(ky+2),
   #           ...ty(ny-ky-1) must be supplied by the user, before entry.
   #           see also the restrictions (ier=10).
   #   c      : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
   #           on succesful exit, c contains the coefficients of the spline
   #           approximation s(x,y)
   #   fp    : real. unless ier=10, fp contains the weighted sum of
   #           squared residuals of the spline approximation returned.
   #   wrk1  : real array of dimension (lwrk1). used as workspace.
   #           if the computation mode iopt=1 is used the value of wrk1(1)
   #           should be left unchanged between subsequent calls.
   #           on exit wrk1(2),wrk1(3),...,wrk1(1+(nx-kx-1)*(ny-ky-1)) will
   #           contain the values d(i)/max(d(i)),i=1,...,(nx-kx-1)*(ny-ky-1)
   #           with d(i) the i-th diagonal element of the reduced triangular
   #           matrix for calculating the b-spline coefficients. it includes
   #           those elements whose square is less than eps,which are treat-
   #           ed as 0 in the case of presumed rank deficiency (ier<-2).
   #   lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
   #           the array wrk1 as declared in the calling (sub)program.
   #           lwrk1 must not be too small. let
   #             u = nxest-kx-1, v = nyest-ky-1, km = max(kx,ky)+1,
   #             ne = max(nxest,nyest), bx = kx*v+ky+1, by = ky*u+kx+1,
   #             if(bx <= by) b1 = bx, b2 = b1+v-ky
   #             if(bx > by) b1 = by, b2 = b1+u-kx  then
   #           lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
   #   wrk2  : real array of dimension (lwrk2). used as workspace, but
   #           only in the case a rank deficient system is encountered.
   #   lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
   #           the array wrk2 as declared in the calling (sub)program.
   #           lwrk2 > 0 . a save upper boundfor lwrk2 = u*v*(b2+1)+b2
   #           where u,v and b2 are as above. if there are enough data
   #           points, scattered uniformly over the approximation domain
   #           and if the smoothing factor s is not too small, there is a
   #           good chance that this extra workspace is not needed. a lot
   #           of memory might therefore be saved by setting lwrk2=1.
   #           (see also ier > 10)
   #   iwrk  : integer array of dimension (kwrk). used as workspace.
   #   kwrk  : integer. on entry kwrk must specify the actual dimension of
   #           the array iwrk as declared in the calling (sub)program.
   #           kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1).
   #   ier   : integer. unless the routine detects an error, ier contains a
   #           non-positive value on exit, i.e.
   #    ier=0  : normal return. the spline returned has a residual sum of
   #             squares fp such that abs(fp-s)/s <= tol with tol a relat-
   #             ive tolerance set to 0.001 by the program.
   #    ier=-1 : normal return. the spline returned is an interpolating
   #             spline (fp=0).
   #    ier=-2 : normal return. the spline returned is the weighted least-
   #             squares polynomial of degrees kx and ky. in this extreme
   #             case fp gives the upper bound for the smoothing factor s.
   #    ier<-2 : warning. the coefficients of the spline returned have been
   #             computed as the minimal norm least-squares solution of a
   #             (numerically) rank deficient system. (-ier) gives the rank.
   #             especially if the rank deficiency which can be computed as
   #             (nx-kx-1)*(ny-ky-1)+ier, is large the results may be inac-
   #             curate. they could also seriously depend on the value of
   #             eps.
   #    ier=1  : error. the required storage space exceeds the available
   #             storage space, as specified by the parameters nxest and
   #             nyest.
   #             probably causes : nxest or nyest too small. if these param-
   #             eters are already large, it may also indicate that s is
   #             too small
   #             the approximation returned is the weighted least-squares
   #             spline according to the current set of knots.
   #             the parameter fp gives the corresponding weighted sum of
   #             squared residuals (fp>s).
   #    ier=2  : error. a theoretically impossible result was found during
   #             the iteration proces for finding a smoothing spline with
   #             fp = s. probably causes : s too small or badly chosen eps.
   #             there is an approximation returned but the corresponding
   #             weighted sum of squared residuals does not satisfy the
   #             condition abs(fp-s)/s < tol.
   #    ier=3  : error. the maximal number of iterations maxit (set to 20
   #             by the program) allowed for finding a smoothing spline
   #             with fp=s has been reached. probably causes : s too small
   #             there is an approximation returned but the corresponding
   #             weighted sum of squared residuals does not satisfy the
   #             condition abs(fp-s)/s < tol.
   #    ier=4  : error. no more knots can be added because the number of
   #             b-spline coefficients (nx-kx-1)*(ny-ky-1) already exceeds
   #             the number of data points m.
   #             probably causes : either s or m too small.
   #             the approximation returned is the weighted least-squares
   #             spline according to the current set of knots.
   #             the parameter fp gives the corresponding weighted sum of
   #             squared residuals (fp>s).
   #    ier=5  : error. no more knots can be added because the additional
   #             knot would (quasi) coincide with an old one.
   #             probably causes : s too small or too large a weight to an
   #             inaccurate data point.
   #             the approximation returned is the weighted least-squares
   #             spline according to the current set of knots.
   #             the parameter fp gives the corresponding weighted sum of
   #             squared residuals (fp>s).
   #    ier=10 : error. on entry, the input data are controlled on validity
   #             the following restrictions must be satisfied.
   #             -1<=iopt<=1, 1<=kx,ky<=5, m>=(kx+1)*(ky+1), nxest>=2*kx+2,
   #             nyest>=2*ky+2, 0<eps<1, nmax>=nxest, nmax>=nyest,
   #             xb<=x(i)<=xe, yb<=y(i)<=ye, w(i)>0, i=1,...,m
   #             lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
   #             kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1)
   #             if iopt=-1: 2*kx+2<=nx<=nxest
   #                         xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
   #                         2*ky+2<=ny<=nyest
   #                         yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
   #             if iopt>=0: s>=0
   #             if one of these conditions is found to be violated,control
   #             is immediately repassed to the calling program. in that
   #             case there is no approximation returned.
   #    ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
   #             space for computing the minimal least-squares solution of
   #             a rank deficient system of linear equations. ier gives the
   #             requested value for lwrk2. there is no approximation re-
   #             turned but, having saved the information contained in nx,
   #             ny,tx,ty,wrk1, and having adjusted the value of lwrk2 and
   #             the dimension of the array wrk2 accordingly, the user can
   #             continue at the point the program was left, by calling
   #             surfit with iopt=1.
c
c further comments:
   #   by means of the parameter s, the user can control the tradeoff
   #    between closeness of fit and smoothness of fit of the approximation.
   #    if s is too large, the spline will be too smooth and signal will be
   #    lost ; if s is too small the spline will pick up too much noise. in
   #    the extreme cases the program will return an interpolating spline if
   #    s=0 and the weighted least-squares polynomial (degrees kx,ky)if s is
   #    very large. between these extremes, a properly chosen s will result
   #    in a good compromise between closeness of fit and smoothness of fit.
   #    to decide whether an approximation, corresponding to a certain s is
   #    satisfactory the user is highly recommended to inspect the fits
   #    graphically.
   #    recommended values for s depend on the weights w(i). if these are
   #    taken as 1/d(i) with d(i) an estimate of the standard deviation of
   #    z(i), a good s-value should be found in the range (m-sqrt(2*m),m+
   #    sqrt(2*m)). if nothing is known about the statistical error in z(i)
   #    each w(i) can be set equal to one and s determined by trial and
   #    error, taking account of the comments above. the best is then to
   #    start with a very large value of s ( to determine the least-squares
   #    polynomial and the corresponding upper bound fp0 for s) and then to
   #    progressively decrease the value of s ( say by a factor 10 in the
   #    beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
   #    approximation shows more detail) to obtain closer fits.
   #    to choose s very small is strongly discouraged. this considerably
   #    increases computation time and memory requirements. it may also
   #    cause rank-deficiency (ier<-2) and endager numerical stability.
   #    to economize the search for a good s-value the program provides with
   #    different modes of computation. at the first call of the routine, or
   #    whenever he wants to restart with the initial set of knots the user
   #    must set iopt=0.
   #    if iopt=1 the program will continue with the set of knots found at
   #    the last call of the routine. this will save a lot of computation
   #    time if surfit is called repeatedly for different values of s.
   #    the number of knots of the spline returned and their location will
   #    depend on the value of s and on the complexity of the shape of the
   #    function underlying the data. if the computation mode iopt=1
   #    is used, the knots returned may also depend on the s-values at
   #    previous calls (if these were smaller). therefore, if after a number
   #    of trials with different s-values and iopt=1, the user can finally
   #    accept a fit as satisfactory, it may be worthwhile for him to call
   #    surfit once more with the selected value for s but now with iopt=0.
   #    indeed, surfit may then return an approximation of the same quality
   #    of fit but with fewer knots and therefore better if data reduction
   #    is also an important objective for the user.
   #    the number of knots may also depend on the upper bounds nxest and
   #    nyest. indeed, if at a certain stage in surfit the number of knots
   #    in one direction (say nx) has reached the value of its upper bound
   #    (nxest), then from that moment on all subsequent knots are added
   #    in the other (y) direction. this may indicate that the value of
   #    nxest is too small. on the other hand, it gives the user the option
   #    of limiting the number of knots the routine locates in any direction
   #    for example, by setting nxest=2*kx+2 (the lowest allowable value for
   #    nxest), the user can indicate that he wants an approximation which
   #    is a simple polynomial of degree kx in the variable x.
c
   #   other subroutines required:
   #     fpback,fpbspl,fpsurf,fpdisc,fpgivs,fprank,fprati,fprota,fporde
c
   #   references:
   #    dierckx p. : an algorithm for surface fitting with spline functions
   #                 ima j. numer. anal. 1 (1981) 267-283.
   #    dierckx p. : an algorithm for surface fitting with spline functions
   #                 report tw50, dept. computer science,k.u.leuven, 1980.
   #    dierckx p. : curve and surface fitting with splines, monographs on
   #                 numerical analysis, oxford university press, 1993.
c
   #   author:
   #     p.dierckx
   #     dept. computer science, k.u. leuven
   #     celestijnenlaan 200a, b-3001 heverlee, belgium.
   #     e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
   #   creation date : may 1979
   #   latest update : march 1987
c
   #   ..
   #   ..scalar arguments..
      real xb,xe,yb,ye,s,eps,fp
      integer iopt,m,kx,ky,nxest,nyest,nmax,nx,ny,lwrk1,lwrk2,kwrk,ier
   #   ..array arguments..
      real x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),
     * c((nxest-kx-1)*(nyest-ky-1)),wrk1(lwrk1),wrk2(lwrk2)
      integer iwrk(kwrk)
   #   ..local scalars..
      real tol
      integer i,ib1,ib3,jb1,ki,kmax,km1,km2,kn,kwest,kx1,ky1,la,lbx,
     * lby,lco,lf,lff,lfp,lh,lq,lsx,lsy,lwest,maxit,ncest,nest,nek,
     * nminx,nminy,nmx,nmy,nreg,nrint,nxk,nyk
   #   ..function references..
      integer max0
   #   ..subroutine references..
   #     fpsurf
   #   .. =#
    @info "surfit start"
#  we set up the parameters tol and maxit.
    maxit = 20
    tol = 0.1e-02
   #   before starting computations a data check is made. if the input data
   #   are invalid,control is immediately repassed to the calling program.
    ier = 10
    @info "checking eps"
    if (eps <= 0. || eps >= 1.)
        @warn "eps out of range: $eps"
        @goto L70 end
   @info "checking kx"        
   if (spline.kx <= 0 || spline.kx > 5) 
        @warn "kx out of range: $(spline.kx)"
        spline.errMsg = "kx out of range: $(spline.kx)"
        @goto L70 end
    kx1 = spline.kx + 1
    @info "checking ky"        
    if (spline.ky <= 0 || spline.ky > 5) 
        @warn "ky out of range: $(spline.ky)"
        spline.errMsg = "ky out of range: $(spline.ky)"
        @goto L70 end
    @info "first 3 checks"    
    ky1 = spline.ky + 1
    kmax = max(spline.kx, spline.ky)
    km1 = kmax + 1
    km2 = km1 + 1
    @info "checking IOPT"
    if (iopt < (-1) || iopt > 1) 
        @warn "iopt invalid: $iopt"
        spline.errMsg = "iopt invalid: $iopt"
        @goto L70 end

    m = min(length(spline.x), length(spline.y), length(spline.z))
   if (m < (kx1 * ky1)) 
      @warn "m < kx1 *ky1  ($m < $kx1 * $ky1)"
      spline.errMsg = "m < kx1 *ky1  ($m < $kx1 * $ky1)"
      @goto L70 
   end
    nminx = 2 * kx1
    nxest =  length(spline.tx)
    nyest =  length(spline.ty)
    if length(spline.c) <   (nxest - spline.kx - 1) * (nyest - spline.ky - 1)
        @warn "c is too small $(length(spline.c)). must be at least $( (nxest - spline.kx - 1) * (nyest - spline.ky - 1))."
        spline.errMsg = "c is too small $(length(spline.c)). must be at least $( (nxest - spline.kx - 1) * (nyest - spline.ky - 1))."
        @goto L70 
    end
    nmax = max(nxest, nyest)  
    if (nxest < nminx || nxest > nmax)
        @warn "nxest out of range $nminx <= $(nxest) <= $nmax"
        spline.errMsg = "nxest out of range $nminx <= $(nxest) <= $nmax"
        @goto L70 end
    nminy = 2 * ky1
    
    if (nyest < nminy || nyest > nmax) 
        @warn "nyest out of range $nminy <= $(nyest) <= $nmax"
        spline.errMsg = "nyest out of range $nminy <= $(nyest) <= $nmax"
        @goto L70 end
      
    nest = max(nxest, nyest)
    nxk = nxest - kx1
    nyk = nyest - ky1
    ncest = nxk * nyk
    nmx = nxest - nminx + 1
    nmy = nyest - nminy + 1
    nrint = nmx + nmy
    nreg = nmx * nmy
    ib1 = spline.kx * nyk + ky1
    jb1 = spline.ky * nxk + kx1
    ib3 = kx1 * nyk + 1
    if (ib1 <= jb1) @goto L10 end
    ib1 = jb1
    ib3 = ky1 * nxk + 1
    @label  L10
    @info "Surfit L10"
    # lwest = ncest * (2 + ib1 + ib3) + 2 * (nrint + nest * km2 + m * km1) + ib3
    kwest = m + nreg
    # if (lwrk1 < lwest || kwrk < kwest) @goto L70 end
    if (spline.xb >= spline.xe || spline.yb >= spline.ye)
      spline.errMsg = "xb >= xe || yb >= ye. ($spline.xb) >= $(spline.xe) || $(spline.yb) >= $(spline.ye)"
       @goto L70 
   end
    for i in 1:m # do 20 i=1,m
        if (spline.w[i] <= 0.) 
         spline.errMsg = "(spline.w[$i] $(spline.w[i]) <= 0."
         @goto L70 end
        if (spline.x[i] < spline.xb || spline.x[i] > spline.xe) 
         spline.errMsg = "x[$i]>xb || x[$i]>xe  $(spline.x[i]) < $(spline.xb) || $(spline.x[i]) > $(spline.xe)"
         @goto L70 end
        if (spline.y[i] < spline.yb || spline.y[i] > spline.ye) 
         spline.errMsg = "y[$i]>yb || y[$i]>ye  $(spline.y[i]) < $(spline.yb) || $(spline.y[i]) > $(spline.xye)"
         @goto L70 end
    end #  @label  L20  continue
    @info "surfit L20"
    if (iopt >= 0) @goto L50 end
    if (spline.nx < nminx || spline.nx > nxest) 
      spline.errMsg = "(spline.nx < nminx || spline.nx > nxest) ($(spline.nx) < $(nminx) || $(spline.nx) > $nxest) "
      @goto L70 end
    nxk = spline.nx - kx1
    spline.tx[kx1] = spline.xb
    spline.tx[nxk + 1] = spline.xe
    for i in kx1:nxk # do 30 i=kx1,nxk
        if (spline.tx[i + 1] <= spline.tx[i]) 
         spline.errMsg = "spline.tx[$(i + 1)] ($(spline.tx[i + 1])) <= spline.tx[$i] ($(spline.tx[i]))"
         @goto L70 end
    end #    @label  L30  continue
    if (spline.ny < nminy || spline.ny > nyest) 
      spline.errMsg = "(spline.ny < nminy || spline.ny > nyest) ($(spline.ny) < $(nminy) || $(spline.ny) > $nyest) "
      @goto L70 end
    nyk = spline.ny - ky1
    spline.ty[ky1] = spline.yb
    spline.ty[nyk + 1] = spline.ye
    for i in ky1:nyk# do 40 i=ky1,nyk
        if (spline.ty[i + 1] <= spline.ty[i]) 
         spline.errMsg = "spline.ty[$(i + 1)] ($(spline.ty[i + 1])) <= spline.ty[$i] ($(spline.ty[i]))"
         
         @goto L70 end
        @label  L40  # continue
        @goto L60
        @label  L50  
        if (s < 0.) 
         spline.errMsg = "s ($s) < 0.0"
         @goto L70 end
        @label  L60  
        ier = 0
   #   we partition the working space and determine the spline approximation
        kn = 1
        ki = kn + m
        lq = 2
        la = lq + ncest * ib3
        lf = la + ncest * ib1
        lff = lf + ncest
        lfp = lff + ncest
        lco = lfp + nrint
        lh = lco + nrint
        lbx = lh + ib3
        nek = nest * km2
        lby = lbx + nek
        lsx = lby + nek
        lsy = lsx + m * km1
        wrk1 = lwork1(
         0.0,  
         # fpint = 
         zeros(Float64,  nrint),
         # coord = 
         zeros(Float64,  nrint),
         # f = 
         zeros(Float64,  ncest),
         # ff = 
         zeros(Float64,  ncest),
         # a = 
         zeros(Float64, ncest, ib1),
         # q = 
         zeros(Float64, ncest, ib3),
         # bx = 
         zeros(Float64, nest, km2),
         # by = 
         zeros(Float64, nest, km2),
         # spx =  
         zeros(Float64, m, km1),
         # spy = 
         zeros(Float64, m, km1),
         # h = 
         zeros(Float64, ib3),
         # index = 
         zeros(Int, nreg),
         # nummer =  
         zeros(Int,  m)
         )
         #    n : integer,wich contains the dimension of q (ncest, ib3).

         #    m : integer, which denotes the bandwidth of q.
        wrk2 = lwork2(
            zeros(Float64, ncest, ib3), # aa(n, m)
            zeros(Float64, ncest),  # ff(n)
            zeros(Float64, ib3)  # h(m)
         )
   
        ier = fpsurf(iopt,m,spline ,s, nxest, nyest, eps,tol,maxit,nest,km1,km2,ib1,ib3,ncest,nrint,
         nreg, wrk1,wrk2; io=io)
        @label  L70 
        @info "Surfit L70. ier = $ier"
        return ier
    end
    return ier
end      
end