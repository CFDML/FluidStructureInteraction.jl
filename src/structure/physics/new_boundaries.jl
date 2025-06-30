## Here we have the boundaries module for thermal and coupling
## We introduced a new way to apply the boundary condition from matrix

#apply_boundary_conditions
function hsource_bc!(f::F, b::Body, name::Symbol) where {F<:Function}
    Peridynamics.check_if_set_is_defined(b.point_sets, name)
    type = Peridynamics.check_boundary_condition_function(f)
    if type === :sdbc
        sdbc = Peridynamics.SingleDimBC(f, :hsource, name, 0x1)
        Peridynamics.add_boundary_condition!(b, b.single_dim_bcs, sdbc)
    elseif type === :pdsdbc
        pdsdbc = Peridynamics.PosDepSingleDimBC(f, :hsource, name, 0x1)
        Peridynamics.add_boundary_condition!(b, b.posdep_single_dim_bcs, pdsdbc)
    end
    return nothing
end

function hsource_databc!(body::Peridynamics.AbstractBody, data::Matrix, name::Symbol,
                              dimspec::Vector{T} = [1]) where {T<:Union{Integer,Symbol}}
    Peridynamics.check_if_set_is_defined(body.point_sets, name)
    Peridynamics.check_databc_data_dimensions(body, data, dimspec)
    dims = Peridynamics.get_dims(dimspec)
    databc = Peridynamics.DataBC(data, :hsource, name, dims)
    Peridynamics.add_boundary_condition!(body, body.data_bcs, databc)
    return nothing
end


function temperature_bc!(f::F, b::Body, name::Symbol) where {F<:Function}
    Peridynamics.check_if_set_is_defined(b.point_sets, name)
    type = Peridynamics.check_boundary_condition_function(f)
    if type === :sdbc
        sdbc = Peridynamics.SingleDimBC(f, :temperature, name, 0x1)
        Peridynamics.add_boundary_condition!(b, b.single_dim_bcs, sdbc)
    elseif type === :pdsdbc
        pdsdbc = Peridynamics.PosDepSingleDimBC(f, :temperature, name, 0x1)
        Peridynamics.add_boundary_condition!(b, b.posdep_single_dim_bcs, pdsdbc)
    end
    return nothing
end

function temperature_databc!(body::Peridynamics.AbstractBody, data::Matrix, name::Symbol,
                              dimspec::Vector{T} = [1]) where {T<:Union{Integer,Symbol}}
    Peridynamics.check_if_set_is_defined(body.point_sets, name)
    Peridynamics.check_databc_data_dimensions(body, data, dimspec)
    dims = Peridynamics.get_dims(dimspec)
    databc = Peridynamics.DataBC(data, :temperature, name, dims)
    Peridynamics.add_boundary_condition!(body, body.data_bcs, databc)
    return nothing
end


function second_bcs!(chunk::Peridynamics.BodyChunk, conv::Vector{Int}, radi::Vector{Int})
    if !isempty(conv)
        (; paramsetup, storage) = chunk
        for i in conv
            temp = storage.temperature[1, i]
            storage.hsource[i] = paramsetup.h * (paramsetup.tem∞ - temp) / paramsetup.dx
        end
    end
        if !isempty(radi)
        for i in conv
            temp = storage.temperature[1, i]
            storage.hsource[i] = paramsetup.hσ * paramsetup.hϵ * (temp^4 - paramsetup.tem∞^4) / paramsetup.dx
        end
    end
end

function find_sec_bcs_points(dh::Peridynamics.AbstractDataHandler)
    bc_points_convection = [Int[] for _ in 1:dh.n_chunks]
    bc_points_radiation = [Int[] for _ in 1:dh.n_chunks]
    for key in keys(dh.chunks[1].psets)
        if occursin("conv", string(key))
            for ic in eachindex(dh.chunks)
                bc_points_convection[ic] = dh.chunks[ic].psets[key]
            end
        end
        
        if occursin("radi", string(key))
            for ic in eachindex(dh.chunks)
                bc_points_radiation[ic] = dh.chunks[ic].psets[key]
            end
        end
    end
    return bc_points_convection, bc_points_radiation
end

function temperature_ic!(b::Body, name::Symbol, value::Real)
    Peridynamics.check_if_set_is_defined(b.point_sets, name)
    sdic = Peridynamics.SingleDimIC(convert(Float64, value), :temperature, name, 0x1)
    Peridynamics.add_initial_condition!(b, b.single_dim_ics, sdic)
    return nothing
end

function temperature_ic!(f::Matrix, b::Body)
    vic = Peridynamics.VIC(f, :temperature)
    Peridynamics.add_initial_condition!(b, b.v_ics, vic)
    return nothing
end

