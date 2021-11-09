"""
given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
gives the value of p such that the rational interpolating function
of the form r(p) = (u*p+v)/(p+w) equals zero at p.
"""
function fprati(p1, f1, p2, f2, p3, f3)
    if p3 <= 0
        #if (p3 > 0.) @goto L10 end
        # value of p in case p3 = infinity.
        p = (p1 * (f1 - f3) * f2 - p2 * (f2 - f3) * f1) / ((f1 - f2) * f3)
    else
        #@goto L20
        # value of p in case p3 ^= infinity.
        #@label  L10
        h1 = f1 * (f2 - f3)
        h2 = f2 * (f3 - f1)
        h3 = f3 * (f1 - f2)
        p = -(p1 * p2 * h3 + p2 * p3 * h1 + p3 * p1 * h2) / (p1 * h1 + p2 * h2 + p3 * h3)
    end
    # adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
    #@label L20  
    if f2 >0.0
        #if (f2 < 0.) @goto L30 end
        p1 = p2
        f1 = f2
    else
        #@goto L40
        #@label L30
        p3 = p2
        f3 = f2
        #@label L40 
    end
    #fprati = p
    return p, p1, p2, p3, f1, f2, f3
end
