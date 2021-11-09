
"""
applies a givens rotation to a and b.
"""
function fprota(vcos, vsin, a, b)
   
    stor1 = a
    stor2 = b
    b = vcos * stor2 + vsin * stor1
    a = vcos * stor1 - vsin * stor2
    return a, b
end
