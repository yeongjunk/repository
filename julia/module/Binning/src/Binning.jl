module Binning
using StatsBase
export binning, binning_unsorted, binning_sorted, binned_mean, binned_histogram, binned_histogram!, binned_hist_coarse 

mutable struct Hist{F,E}
    weights::Array{F}
    edges::E
    function Hist{F,E}(weights::Vector{F}, edges::AbstractArray) where {F,E}
        new{F,E}(weights, edges) 
    end
end

function normalize!(h::Hist)
    normalize!(h.weights)
end


@doc """
Bin x of given bin edge given by x_edges
"""
@inline function binning(x::T, x_edges::AbstractArray{T}) where T <: Number
    L = 1 #leftmost edge
    R = length(x_edges) #rightmost edge
    #j = 1;
    if x_edges[1] <= x < x_edges[end]
        while L+1 < R
            m = (L+R)÷2
            # println("$(j),$(L), $(m), $(R)")
            if x  >= x_edges[m]
                 L = m
            else
                R = m
            end
            #j += 1
       end
        return L
    else
        return 0
    end
end

@doc """
Binning sorted array.
"""
function binning_sorted(x::AbstractArray{T}, x_edges::AbstractArray{T}) where T <: Number
    lbl = Array{Int64}(undef, length(x))
    if x[end] < x_edges[1] || x[1] > x_edges[end]
        fill!(lbl, 0)
        return lbl
    else
        j = 1
        while x[j] < x_edges[1]
            j += 1
        end
        lbl[1:j-1] .= 0

        for i in 1:length(x_edges)-1
            while j <= length(x) && x_edges[i] <= x[j] < x_edges[i+1]
                lbl[j] = i
                j += 1
            end
        end
        if j <= length(x)
            lbl[j:end] .= 0
        end
        return lbl
    end
end

function binning_unsorted(x::AbstractArray, x_edges::AbstractArray)
    lbl = binning.(x, Ref(x_edges))
end

function binning(x::AbstractArray, x_edges::AbstractArray)
    if issorted(x)
        lbl = binning_sorted(x, x_edges)
    else
        lbl = binning_unsorted(x, x_edges)
    end
    return lbl
end
"""
Compute the x-resolved mean, standard deviations, and numbers(histogram) of y = f(x), using the edge of the bins x_edges
"""
function binned_mean(x::AbstractArray, y::AbstractArray{F}, x_edges::AbstractArray) where F
    n = length(x_edges)
    lbl = binning(x, x_edges)
    y_mean = Array{F}(undef, n-1)
    y_std = similar(y_mean)
    y_num = Array{Int64}(undef, n-1)
    for i in 1:n-1
        idx = findall(x -> x == i, lbl)            
        y_mean[i] = mean(y[idx]) 
        y_num[i] = length(y[idx])
        y_std[i] = std(y[idx])
    end
    return y_mean, y_std, y_num
end

function binned_histogram(x::AbstractArray, x_edges::AbstractArray)
    y = zeros(Int64, length(x_edges)- 1)
    binned_histogram!(y, x, x_edges)
    return y
end

"""
Update the histogram y of given edges of sorted array
"""
function binned_histogram!(y::AbstractArray, x::AbstractArray, x_edges::AbstractArray)
    @assert length(y) == length(x_edges) - 1 
    lbl = binning(x, x_edges)
    filter!(t -> t != 0, lbl)
    prev_ind = 1
    for i in 1:length(x_edges) - 1
        ind = searchsortedlast(lbl, i)
        y[i] += (ind - prev_ind + 1)
        prev_ind = ind + 1
    end
end

"""
Convert the histogram y, binned over x_edges, into coarser histogram with x_edges_2. This is only possible if x_edges_2 ⊆ x_edges
"""    
function binned_hist_coarse(y, x_edges, x_edges_2)
    @assert x_edges_2 ⊆ x_edges 
    n = length(x_edges_2)
    y_coarse = zeros(eltype(y), n-1) 
    for i in 1:n-1 
        idx = findall(t -> (t >= x_edges_2[i]) && (x_edges_2[i+1] > t), x_edges)
        for j in 1:length(idx)
            y_coarse[i] += y[idx[j]]
        end
    end
    return y_coarse
end
end
