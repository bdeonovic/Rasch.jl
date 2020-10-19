function haberman(m::D, s::AbstractMatrix{T}, f::AbstractVector{D}, 
                  β0::AbstractVector{T} = zeros(T, size(s)[2]*(m-1)),
                  method::Optim.SecondOrderOptimizer = Optim.Newton(), 
                  options=Optim.Options()) where {T <: Real, D <: Integer}

  Q = size(s)[2]
  δ = zeros(T, m, Q)
  γ = zeros(T, m+1, m+1)
  η = zeros(T, m, m, m+1)

  pr = zeros(T, m, m+1)

  πs = zeros(T, m, m+1)
  πd = zeros(T, m, m, m+1)

  tmp = zeros(T, m-1)
  last_β = rand(T, Q*(m-1))
  function calculate_common!(x, last_x)
    if x != last_x
      copyto!(last_x,x)
      @inbounds for q in 1:Q
        δ[1,q] = 1 
        δ[2:end,q] .= exp.(x[((q-1)*(m-1)+1):(q*(m-1))])
      end

      @inbounds for r in 0:m
        pr[1,r+1] = 0.5
        tmp .= 0.0
        for q in 1:Q
          @views tmp .+= r^(q-1) .* x[((q-1)*(m-1)+1):(q*(m-1))] 
        end
        @views pr[2:end,r+1] .= 1 ./ (1 .+ exp.(-tmp))
        @views poisbin_sum_taub!(γ[:,r+1], pr[:,r+1])
        @views poisbin_sum_taub_dervs_2!(η, pr[:,r+1])

        for i in 1:m
          πs[i,r+1] = r == 0 ? 0.0 : pr[i,r+1] * η[i,i,r] / γ[r+1,r+1] 
          for j in 1:m
            if r == 0
              πd[i,j,r+1] = 0.0
            elseif r == 1
              πd[i,j,r+1] = -pr[i,r+1] * pr[j,r+1] * f[r] * η[i,i,r] * η[j,j,r] / γ[r+1,r+1]^2 
            else
              πd[i,j,r+1] = pr[i,r+1] * pr[j,r+1] * f[r] * (η[i,j,r-1] / γ[r+1,r+1] - η[i,i,r] * η[j,j,r] / γ[r+1,r+1]^2)
            end
          end
        end
      end
    end
  end

  function neglogLC(β)
    calculate_common!(β, last_β)
    result = 0.0
    @inbounds for q in 1:Q
      @views result -= s[:,q]'log.(δ[:,q])
    end
    @inbounds for r in 1:m
      result += f[r] * (log(γ[r+1,r+1]) - log(γ[1,r+1]))
    end
    return result
  end

  function neg_g!(storage, β)
    calculate_common!(β, last_β)
    @inbounds for i in 1:(m-1)
      for q in 1:Q
        storage[i + (q-1)*(m-1)] = -s[i+1,q]
      end
      for r in 1:m
        A = f[r] * πs[i+1,r+1]
        for q in 1:Q
          storage[i + (q-1)*(m-1)] += r^(q-1) * A
        end
      end
    end
  end

  function neg_h!(storage, β)
    calculate_common!(β, last_β)
    storage .= 0.0
    @inbounds for r in 1:m
      for i in 1:(m-1)
        A = f[r] * (πs[i+1,r+1] - πs[i+1,r+1]^2) 
        for q1 in 1:Q
          for q2 in 1:Q
            storage[i + (q1-1)*(m-1),i + (q2-1)*(m-1)] += r^(q1-1 + q2-1) * A
          end
        end
        for j in (i+1):(m-1)
          B = πd[i+1,j+1,r+1]
          for q1 in 1:Q
            for q2 in 1:Q
              storage[i + (q1-1)*(m-1),j + (q2-1)*(m-1)] += r^(q1-1 + q2-1) * B
              storage[j + (q2-1)*(m-1),i + (q1-1)*(m-1)] += r^(q1-1 + q2-1) * B
            end
          end
        end
      end
    end
  end

  df = TwiceDifferentiable(neglogLC, neg_g!, neg_h!, β0)
  res = Optim.optimize(df, β0, method, options)
  delta_hat = zeros(m, Q)

  @inbounds for q in 1:Q
    delta_hat[1,q] = 1.0
    delta_hat[2:m,q] .= exp.(res.minimizer[((q-1)*(m-1)+1):(q*(m-1))])
  end

  return (res, delta_hat, πs)
end