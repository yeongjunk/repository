module LatticeDraw
using LaTeXStrings
using Lattice
using Plots

export get_coordinate, lattice_wave_data_2d, lattice_wave_plot_2d, isconnected, lattice_ham_plot_2d
function get_coordinate(ltc::Lattice2D)
    ms = Int64[]
    ns = Int64[]
    for i in 1:ltc.N*ltc.M*ltc.U
        m,n, _ = site.(Ref(ltc), i)
        push!(ms, m)
        push!(ns, n)
    end
    return hcat(ms, ns)
end

function lattice_wave_data_2d(ltc::Lattice2D, psi)
    mn =  get_coordinate(ltc)
    df = DataFrame(x = mn[:,1], y=mn[:,2], psi_abs = abs.(psi))
    return df
end

function lattice_wave_plot_2d(ltc::Lattice2D, psi::Vector)
    df = latticedata_2d(ltc, psi)
    p = scatter(df.x, df.y, marker_z = df.psi_abs, 
        msw = 0, c = :lajolla, frame = :box, label = false, 
        markersize = 5.5, colorbar_title = L"\psi", size = (500, 500), thickness_scaling = 1.2, format = :png, dpi = 300)
    xlabel!(L"x")
    ylabel!(L"y")
    return p
end

function isconnected(x; tol = 1E-12)
    return x > 1E-12
end

function lattice_ham_plot_2d(ltc, H; tol = 1E-12)
    mn = get_coordinate(ltc::Lattice2D) 
    p = plot(legend = false, format = :png);

    for i in 1:size(H, 1)
        for j in i:size(H, 1)
            if isconnected(H[i,j], tol=tol)
                c1 = site(ltc, i)[1:end-1]
                c2 = site(ltc, j)[1:end-1]
                plot!(p, c1, c2)
            end
        end
    end
    scatter!(p, mn[:, 1], mn[:, 2])
    return p
end
end
