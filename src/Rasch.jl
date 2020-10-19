module Rasch
    using Optim
    using LineSearches

    import ElementarySymmetricFunctions: poisbin_sum_taub!, poisbin_sum_taub_dervs_2!

    function logsumexp(x::AbstractVector{T}) where T <: Real
        u, idx = findmax(x)
        u + log1p(sum(exp.(x[1:end .!= idx] .- u)))
    end
    function logsumexpfast(x::AbstractVector{T}) where T <: Real
        u = maximum(x)
        u + log(sum(exp.(x .- u)))
    end

    include("haberman.jl")
    export haberman

end
