module LatticeGFH90

import Nemo: FiniteField

export toto, FiniteField

function toto(x)
    return x+1
end

function FiniteField(y::Int)
    return y^2
end

end
