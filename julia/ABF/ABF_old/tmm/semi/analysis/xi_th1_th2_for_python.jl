using MAT

function py_form(x,y)
    data = Array{Float64}(undef, length(θI), 2)
    for i in 1:length(x)
            data[i,1] = x[i]
            data[i,2] = y[i]
    end
    return data
end

function py_form(x,y,z)
    data = Array{Float64}(undef, length(θI), 2)
    for i in 1:length(x)
            data[i,1] = x[i]
            data[i,2] = y[i]
    end
    return data
end

# function py_form2(x, y, z)
#     data = Array{Float64}(undef, length(θI)*length(θII), 3)
#     for i in 1:length(x)
#         for j in 1:length(y)
#             data[(i-1)*length(θII)+j,1] = x[i]
#             data[(i-1)*length(θII)+j,2] = y[j]
#             data[(i-1)*length(θII)+j,3] = z[i,j]
#         end
#     end
#     return data
# end


include("xi_abf2.jl")

#data save directory
dir_weak = "/Users/pcs/data/ABF/rawdata/nu2-tm-semi-detangled-weak-log/"
dir_full = "/Users/pcs/data/ABF/rawdata/nu2-tm-semi-detangled/"
#read all the data and the parameter
ξ_weak, params_weak = ξ_read(dir_weak*"/config-weak-xi_all.mat")
(E, W, θI, θII) =params_expand(params_weak)

weak_xi_th1_th2_W2 = py_form2(θI, θII,ξ_weak[1,3,:,:])
weak_xi_th1_th2_W4 = py_form2(θI, θII,ξ_weak[1,2,:,:])
weak_xi_th1_th2_W6 = py_form2(θI, θII,ξ_weak[1,1,:,:])

weak_xi_th_W2 = py_form(θI, [ξ_weak[1,3,i,i] for i in 1:size(ξ_weak,3)])
weak_xi_th_W4 = py_form(θI, [ξ_weak[1,3,i,i] for i in 1:size(ξ_weak,3)])
weak_xi_th_W6 = py_form(θI, [ξ_weak[1,3,i,i] for i in 1:size(ξ_weak,3)])



ξ, params = ξ_read(dir_full*"/config-xi_all.mat")
(E, W, θI, θII) =params_expand(params)


xi_th1_th2_1 = py_form2(E, W, ξ[:,:,1,1])
xi_th1_th2_2 = py_form2(E, W, ξ[:,:,9,9])
xi_th1_th2_3 = py_form2(E, W, ξ[:,:,17,17])
xi_th1_th2_4 = py_form2(E, W, ξ[:,:,25,25])


matdata = Dict();
matdata["weak_xi_th1_th2_W2"] = weak_xi_th1_th2_W2
matdata["weak_xi_th1_th2_W4"] = weak_xi_th1_th2_W4
matdata["weak_xi_th1_th2_W6"] = weak_xi_th1_th2_W6

matdata["weak_xi_th_W2"] = weak_xi_th_W2
matdata["weak_xi_th_W4"] = weak_xi_th_W4
matdata["weak_xi_th_W6"] = weak_xi_th_W6

matdata["xi_th1_th2_1"] = xi_th1_th2_1
matdata["xi_th1_th2_2"] = xi_th1_th2_2
matdata["xi_th1_th2_3"] = xi_th1_th2_3
matdata["xi_th1_th2_4"] = xi_th1_th2_4


fn_out = dir_weak*"/py_form.mat"
matwrite(fn_out, matdata, compress=true)
 # write compressed MATLAB file
