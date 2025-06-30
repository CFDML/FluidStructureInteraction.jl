### ablation function
function update_ablation_exist!(chunk::Peridynamics.AbstractBodyChunk, Δt::Float64)
    new_ablation_point_idx = Vector{Int}()
    for i in eachindex(chunk.ch.loc_points)
        if chunk.storage.exist[i] == 1
            Δd = cal_ablation_rate(chunk.storage.temperature[i])
            chunk.storage.ablation_deep[i] += Δd * Δt
    
            if chunk.storage.ablation_deep[i] >= chunk.paramsetup.dx
                chunk.storage.exist[i] = 0
                push!(new_ablation_point_idx, i)
            end
        end
    end
    return new_ablation_point_idx
end

function cal_ablation_rate(temperature::Float64) # any form you prefer! :)
    tk = 1500.0
    if temperature <= tk 
        Δd = 0.0
    else 
        Δd = (temperature - tk) / tk * Δx * 100
    end
    return Δd
end

function update_ab_bcs_add!(chunk::Peridynamics.AbstractBodyChunk, new_ablation_point_idx::Vector{Int})
    if !isempty(new_ablation_point_idx)
        position = chunk.system.position[:, 1:chunk.ch.n_loc_points]
        #position = chunk.storage.position[:, 1:chunk.ch.n_loc_points]
        n_radius = 1.1 * chunk.paramsetup.δ
        new_add_bcs_idx = find_neighbors_indices(position, new_ablation_point_idx, n_radius)

        for id in new_add_bcs_idx
            θx = position[1, id]/sqrt(position[1, id]^2 + position[2, id]^2 + position[3, id]^2)
            chunk.storage.hsource[1, id] = q_const * abs(θx) * chunk.storage.exist[id]
            chunk.storage.b_ext[1, id] = f_const * abs(θx) * chunk.storage.exist[id]
            id_ch = chunk.ch.loc_points[id]
            q_matrix[1, id_ch] = q_const * abs(θx) * chunk.storage.exist[id]
            f_matrix[1, id_ch] = f_const * abs(θx) * chunk.storage.exist[id]
        end
        for id in new_ablation_point_idx
            chunk.storage.hsource[1, id] = q_const * chunk.storage.exist[id]
            chunk.storage.b_ext[1, id] = f_const * chunk.storage.exist[id]
            chunk.storage.velocity_half[1, id] = 0.0
            chunk.storage.velocity_half[2, id] = 0.0

            id_ch = chunk.ch.loc_points[id]
            q_matrix[1, id_ch] = q_const * chunk.storage.exist[id]
            f_matrix[1, id_ch] = f_const * chunk.storage.exist[id]
        end

    end
end

function find_neighbors_indices(position::Matrix{Float64}, new_ablation_point_idx::Vector{Int}, radius::Float64)
    nhs = GridNeighborhoodSearch{3}(search_radius=radius, n_points=size(position, 2))
    initialize_grid!(nhs, position)

    all_neighbors_indices = Int[]

    for point_id in new_ablation_point_idx
        foreach_neighbor(position, position, nhs, point_id; search_radius=radius) do i, j, _, L
            if i != j
                push!(all_neighbors_indices, j)
            end
        end
    end
    return all_neighbors_indices
end

function update_mech_state!(chunk::Peridynamics.AbstractBodyChunk)
    for i in eachindex(chunk.ch.loc_points)
        if chunk.storage.exist[i] < 1
            chunk.storage.velocity[:, i] .= 0.0
            chunk.storage.velocity_half[:, i] .= 0.0
            chunk.storage.velocity_half_old[:, i] .= 0.0
            chunk.storage.displacement[:, i] .= 0.0
        end
    end
end


