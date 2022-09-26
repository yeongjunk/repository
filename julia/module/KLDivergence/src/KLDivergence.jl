module KLDivergence 

export KL1

function KL1(ψ_1::AbstractVector{F}, ψ_2::AbstractVector{F}) where F
    @assert length(ψ_1) == length(ψ_2)
    kl = 0.0
    for i in 1:length(ψ_1)
        kl += abs2(ψ_1[i]) * log(abs(ψ_1[i]/ψ_2[i]))
    end
    return kl
end

function KL1(ψs::Matrix{F}) where F
    kls = Array{Float64}(undef, size(ψs, 2)-1) 
    for i in 1:size(ψs, 2)-1
        @views kls[i] = KL1(ψs[:, i], ψs[:, i+1])
    end

    return kls
end

end # module
