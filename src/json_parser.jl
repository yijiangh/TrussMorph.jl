using TrussMorph
tm = TrussMorph
import JSON
using Dates
using Printf

function parse_truss_json(file_path::String; parse_morph=false)
    data = Dict()
    open(file_path, "r") do f
        data_txt = read(f, String)
        data = JSON.parse(data_txt)
    end
    dim = data["dimension"]
    n_nodes = data["node_num"]
    n_elements = data["element_num"]

    # get material properties
    # pressure: kN/cm^2 -> kN/m^2
    # density: kN/m^3
    mp = tm.MaterialProperties(data["material_properties"]["material_name"],
                            data["material_properties"]["youngs_modulus"] * 1e4,
                            data["material_properties"]["shear_modulus"] * 1e4,
                            data["material_properties"]["poisson_ratio"],
                            data["material_properties"]["density"])

    X = Matrix{Float64}(undef, n_nodes,2)
    T = Matrix{Int}(undef, n_elements,2)
    fix_node_ids = []

    # get node coord
    for i=1:n_nodes
        X[i,:] = hcat(data["node_list"][i]["point"]["X"],
                      data["node_list"][i]["point"]["Y"])
                      # data["node_list"][i]["point"]["Z"])
        if 1 == data["node_list"][i]["is_grounded"]
            push!(fix_node_ids, i)
        end
    end

    # get fixities
    fix_dof = [1,2,6]
    if 2 != dim
        fix_dof = 1:1:6
    end

    S = Matrix{Int}(undef, length(fix_node_ids),length(fix_dof)+1)
    for i=1:length(fix_node_ids)
        S[i,1] = fix_node_ids[i]
        S[i,2:end] = data["node_list"][fix_node_ids[i]]["fixities"][fix_dof]'
    end

    # get element node ids
    for i=1:n_elements
        T[i,:] = (data["element_list"][i]["end_node_ids"] .+ 1)'
    end

    morph_data = Dict()
    if parse_morph
        morph_data = data["morph_data"]
        # @show morph_data
    end

    return Truss(X, T, S, mp), morph_data
end

function parse_load_json(file_path::String, node_dof::Int)
    data = Dict()
    open(file_path, "r") do f
        data_txt = read(f, String)
        data = JSON.parse(data_txt)
    end

    dim = data["dimension"]
    n_load_nodes = length(data["point_load_list"])
    @assert(dim == 2)

    Load = zeros(n_load_nodes, 1+node_dof)
    for i=1:n_load_nodes
        Load[i,1] = data["point_load_list"][i]["applied_node_id"] + 1
        if 2 == dim
            Load[i,2] = data["point_load_list"][i]["Fx"]
            Load[i,3] = data["point_load_list"][i]["Fy"]
            if 3 == node_dof
                Load[i,4] = data["point_load_list"][i]["Mz"]
            end
        end
    end

    @assert(n_load_nodes > 0)
    return Load
    # TODO: include_self_weight
end

function save_morph_path_json(morph_path::Array{Matrix{Float64}}, file_dir::String, st_file_name::String, end_file_name::String, t0::TrussMorph.Truss, sE::Vector{Float64}, wE::Vector{Float64}, totE::Vector{Float64}, sumE::Float64, parm_smooth::Float64, parm_weight::Float64)
    pure_st_file_name = SubString(st_file_name, 1:length(st_file_name)-length(".json"))
    pure_end_file_name = SubString(end_file_name, 1:length(end_file_name)-length(".json"))
    result_file_name = pure_st_file_name * "-" *  pure_end_file_name
    f_file_dir = joinpath(file_dir, result_file_name)
    if !ispath(f_file_dir)
        mkpath(f_file_dir)
    end

    # same material
    mp_data = Dict()
    mp_data["material_name"] = t0.mp.name
    mp_data["youngs_modulus"] = t0.mp.E * 1e-4
    mp_data["youngs_modulus_unit"] = "kN/cm2"
    mp_data["shear_modulus"] = t0.mp.G * 1e-4
    mp_data["shear_modulus_unit"] = "kN/cm2"
    mp_data["tensile_yeild_stress"] = "N/A"
    mp_data["tensile_yeild_stress_unit"] = "kN/cm2"
    mp_data["density"] = t0.mp.ρ
    mp_data["density_unit"] = "kN/m3"
    mp_data["poisson_ratio"] = t0.mp.μ

    mp_data["radius"] = "N/A"
    mp_data["radius_unit"] = "centimeter"

    # same topology
    topo_data = Dict[]
    for e=1:size(t0.T,1)
        e_data = Dict()
        e_data["end_node_ids"] = t0.T[e,:]
        e_data["element_id"] = e-1
        e_data["layer_id"] = 0
        push!(topo_data, e_data)
    end

    # morph data
    morph_data = Dict()
    morph_data["smooth_parameter"] = parm_smooth
    morph_data["weight_parameter"] = parm_weight
    morph_data["time_in_path"] = string(0,"/",0)
    morph_data["weight_energy"] = 0.0
    morph_data["smoothness_energy"] = 0.0
    morph_data["tot_energy"] = 0.0

    plen = length(morph_path)
    for i=1:plen
        data = Dict()
        data["model_name"] = result_file_name * "_mp" * string(i,"-",plen)
        data["model_type"]= "2D_frame"
        data["unit"] = "meter"
        data["generate_time"] = string(Dates.now())
        data["dimension"] = size(morph_path[i],2)

        i_morph_data = morph_data #deepcopy
        i_morph_data["time_in_path"] = string(i,"/",plen)
        i_morph_data["smoothness_energy"] = sE[i]
        i_morph_data["weight_energy"] = wE[i]
        i_morph_data["tot_energy"] = totE[i]
        i_morph_data["sum_energy"] = sumE
        data["morph_data"] = i_morph_data

        data["node_num"] = size(morph_path[i],1)
        data["element_num"] = size(t0.T,1)

        data["material_properties"] = mp_data
        data["node_list"] = Dict[]
        for j=1:size(morph_path[i],1)
            pt_data = Dict()
            pt_data["point"] = Dict()
            pt_data["point"]["X"] = morph_path[i][j,1]
            pt_data["point"]["Y"] = morph_path[i][j,2]
            pt_data["node_id"] = j-1

            pt_fix = findall(x->x==j, t0.S[:,1])
            pt_data["is_grounded"] = !isempty(pt_fix)
            if !isempty(pt_fix)
                @assert(length(pt_fix) == 1)
                pt_data["fixities"] = ones(Int, 6)
                if 2 == size(morph_path[i],1)
                    pt_data["fixities"][1] = t0.S[pt_fix,1]
                    pt_data["fixities"][2] = t0.S[pt_fix,2]
                    pt_data["fixities"][6] = t0.S[pt_fix,3]
                end
            else
                pt_data["fixities"] = []
            end
            push!(data["node_list"], pt_data)
        end
        data["element_list"] = topo_data

        # write
        stringdata = JSON.json(data)
        result_json_name = result_file_name * "_" * string(@sprintf("%03d",i),"-",plen) * ".json"
        open(joinpath(f_file_dir,result_json_name), "w") do f
                write(f, stringdata)
        end
    end
end
