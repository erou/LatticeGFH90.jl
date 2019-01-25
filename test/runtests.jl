using LatticeGFH90, Test

include("embeddings-test.jl")

function test_all()
    test_embeddings()
end

test_all()
