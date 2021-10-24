function fprota(vcos, vsin, a, b)
   #  subroutine fprota applies a givens rotation to a and b.
   #  ..
   #  ..scalar arguments..
      # real cos,sin,a,b
# c ..local scalars..
 #     real stor1,stor2
   #  ..
    stor1 = a
    stor2 = b
    b = vcos * stor2 + vsin * stor1
    a = vcos * stor1 - vsin * stor2
    return a, b
end
