"""
calculates the parameters of a givens transformation .
args:
piv: float, the pivot value, first non zero element in vector 1
ww: float, the element i  vector 2 that needs to become zero

returns:
- vcos: the cosine value for the ,atrix rotations
- vsin : the sine value for the matrix rotation

referntce:  Diercks p 56
"""
function fpgivs(piv, ww)
  #@debug "fpgivs(piv=$piv, ww=$ww)" 
    one = 0.1e+01
    oldww = ww
    @assert ! isnan(piv)
    store = abs(piv)
    if (store >= ww) 
      dd = store * sqrt(one + (ww / piv)^2) 
    else
      dd = ww * sqrt(one + (piv / ww)^2) 
    end
    v_cos = ww / dd
    v_sin = piv / dd
    ww = dd
    #@debug "fpgivs returns ww=$oldww->$ww, v_cos=$v_cos, v_sin = $v_sin"
    return  v_cos, v_sin, ww
end
    