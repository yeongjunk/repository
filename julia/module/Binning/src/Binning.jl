module Binning

    export binning, binning_unsorted, binning_sorted

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

    function binning_unsorted(x::AbstractArray, x_edges)
        lbl = binning.(x, Ref(x_edges))
    end

    function binning(x::AbstractArray, x_edges)
        if issorted(x)
            lbl = binning_sorted(x, x_edges)
        else
            lbl = binning_unsorted(x, x_edges)
        end
        return lbl
    end
end
