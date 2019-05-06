import JSON

function parse_truss_json(file_path::String)
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
    mp = MaterialProperties(data["material_properties"]["youngs_modulus"] * 1e4,
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

    return Truss(X, T, S, mp)
end

function parse_load_json(file_path::String)
    data = Dict()
    open(file_path, "r") do f
        data_txt = read(f, String)
        data = JSON.parse(data_txt)
    end

    dim = data["dimension"]
    n_load_nodes = length(data["point_load_list"])

    if 2 == dim
        node_dof = 3
    else
        node_dof = 6
    end

    Load = zeros(n_load_nodes, 1+node_dof)
    for i=1:n_load_nodes
        Load[i,1] = data["point_load_list"][i]["applied_node_id"] + 1
        if 2 == dim
            Load[i,2] = data["point_load_list"][i]["Fx"]
            Load[i,3] = data["point_load_list"][i]["Fy"]
            Load[i,4] = data["point_load_list"][i]["Mz"]
        end
    end

    @assert(n_load_nodes > 0)
    return Load
    # TODO: include_self_weight
end
