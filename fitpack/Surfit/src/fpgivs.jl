function fpgivs(piv, ww)
# c  subroutine fpgivs calculates the parameters of a givens
# c  transformation .
# c  ..
# c  ..scalar arguments..
      # real piv,ww,cos,sin
# c  ..local scalars..
#      real dd,one,store
# c  ..function references..
      # real abs,sqrt
# c  ..
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
    