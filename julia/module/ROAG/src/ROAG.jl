# Level spacing statistics

module ROAG
    export levelspacing!, roag!

    @doc """
    Adjacent level spacing of E. Overwrites the output to E.
    """
    function levelspacing!(E::AbstractArray{Float64, 1}; sorted = true)
        if !sorted
            sort!(E)
        end
        E[1:end-1] .= diff(E)
        pop!(E)
    end


    @doc """
    Ratio of adjacent gaps(roag) of eigenvalues. Overwrites the output to E.
    """
    function roag!(E::AbstractArray{Float64, 1}; sorted = true)
        levelspacing!(E; sorted = sorted)
        for α = 1:(length(E)-1)
            E[α] = min(E[α],E[α+1])/max(E[α],E[α+1])
        end
        pop!(E)
    end
end
