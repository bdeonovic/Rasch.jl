using Rasch
using Optim
using Test
using DelimitedFiles
using LinearAlgebra

@testset "Rasch.jl" begin
    X = readdlm(joinpath(dirname(pathof(Rasch)), "../test/lsat.txt"));
    n, m = size(X)
    r = vec(mapslices(sum, X, dims=2))
    s = Float64.(vec(mapslices(sum, X, dims=1)))
    f = [sum(r.==j) for j in 1:m]
    t = Float64.(vec(mapslices(sum, X .* r, dims=1)))
    q = Float64.(vec(mapslices(sum, X .* r .^2, dims=1)))

    res0, delta_hat0, πs0 =  haberman(m,reshape(s,(m,1)),f,zeros(m-1),Optim.Newton(), Optim.Options(f_tol=1e-6, g_tol=0))
    res1, delta_hat1, πs1 =  haberman(m,hcat(s,t),f,zeros(2m-2),Optim.Newton(), Optim.Options(f_tol=1e-6, g_tol=0))
    res2, delta_hat2, πs2 =  haberman(m,hcat(s,t,q),f,zeros(3m-3),Optim.Newton(), Optim.Options(f_tol=1e-6, g_tol=0)).
end
