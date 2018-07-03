function kl2(dist::MvNormal,
            samps::Array{Float64,2},window = 3)
    μ1 = params(dist)[1]
    μ2 = vec(mean(samps,2))
    d = length(μ1)
    d2 = length(μ2)
    Σ1 = params(dist)[2] * eye(d)
    Σ2 = cov(samps,2)

    print(size(μ1),size(μ2),size(Σ1),size(Σ2))
    nwind = div(d,window)
    subdivs = zeros(nwind)
    subtr = zeros(nwind)
    subdif = zeros(nwind)
    sublog = zeros(nwind)
    for w = 1:nwind
        r = ((w-1) * window + 1):(w * window)
        diff = (μ2[r] - μ1[r])
        subtr[w] = (trace(inv(Σ2[r,r]) * Σ1[r,r]) - window) / 2
        subdif[w] = diff' * Σ2[r,r] * diff / 2
        d2, d1 = det(Σ2[r,r]), det(Σ1[r,r])
        if (d2 <= 0)
            print("kl2 non-positive determinant", d1, " ", d2, " ")
            print(Σ1[r,r])
            print(Σ2[r,r])
            d2 = d1
        end
        sublog[w] = (log(d2/d1)) / 2
        subdivs[w] = subtr[w] + subdif[w] + sublog[w]
    end
    #print("KL parts:",mean(subtr),",",mean(subdif),",",mean(sublog),",",mean(subdivs))
    [mean(subdivs),mean(subtr),mean(subdif),mean(sublog)]
end
