function solve_therm!(dh::Peridynamics.AbstractDataHandler, vv::Thermstep,
                  options::Peridynamics.AbstractJobOptions)
    Peridynamics.export_reference_results(dh, options)
    Δt = vv.Δt
    t_factors = modify_t(dh)
    conv, radi = find_sec_bcs_points(dh)
    if mpi_isroot()
        p = Progress(vv.n_steps; dt=1, desc="TIME INTEGRATION LOOP", color=:normal,
                     barlen=40, enabled=Peridynamics.progress_bars())
    end
    for n in 1:vv.n_steps
        th_timestep!(dh, options, Δt, n, t_factors, conv, radi)
        Peridynamics.mpi_isroot() && next!(p)
    end
    Peridynamics.mpi_isroot() && Peridynamics.finish!(p)
    return dh
end
#=
function th_timestep!(dh::Peridynamics.AbstractThreadsBodyDataHandler, options::Peridynamics.AbstractJobOptions,
                          Δt::Float64, n::Int, t_factors::Vector{Vector{Float64}},
                          conv::Vector{Vector{Int}}, radi::Vector{Vector{Int}})
    t = n * Δt
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        Peridynamics.apply_boundary_conditions!(chunk, t)
        second_bcs!(chunk, conv[chunk_id], radi[chunk_id])
        update_temperature!(chunk, Δt)
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_loc_to_halo!(dh, chunk_id)
        calc_pflux!(dh.chunks[chunk_id], t_factors[chunk_id]) 
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_halo_to_loc!(dh, chunk_id)
        Peridynamics.export_results(dh, options, chunk_id, n, t)
    end
    #
    return nothing
end


function verlet_timestep!(dh::Peridynamics.AbstractMPIBodyDataHandler, options::Peridynamics.AbstractJobOptions,
                          Δt::Float64, Δt½::Float64, n::Int)
    t = n * Δt
    chunk = dh.chunk
    @timeit_debug TO "update_vel_half!" update_vel_half!(chunk, Δt½)
    @timeit_debug TO "apply_boundary_conditions!" apply_boundary_conditions!(chunk, t)
    @timeit_debug TO "update_disp_and_pos!" update_disp_and_pos!(chunk, Δt)
    @timeit_debug TO "exchange_loc_to_halo!" exchange_loc_to_halo!(dh)
    @timeit_debug TO "calc_force_density!" calc_force_density!(chunk)
    @timeit_debug TO "exchange_halo_to_loc!" exchange_halo_to_loc!(dh)
    @timeit_debug TO "calc_damage!" calc_damage!(chunk)
    @timeit_debug TO "update_acc_and_vel!" update_acc_and_vel!(chunk, Δt½)
    @timeit_debug TO "export_results" export_results(dh, options, n, t)
    return nothing
end


function Peridynamics.req_point_data_fields_timesolver(::Type{Thermstep})
    fields = (:position, :temperature, :pflux, :hsource)
    return fields
end

function Peridynamics.req_data_fields_timesolver(::Type{Thermstep})
    return ()
end
=#

function fsi_log_timesolver(options::Peridynamics.AbstractJobOptions, fv::Flowstep, vv::Thermstep)
    msg = "VELOCITY VERLET TIME FUILD SOLVER\n"
    msg *= Peridynamics.msg_qty("number of time steps", fv.n_steps)
    msg *= Peridynamics.msg_qty("time step size", fv.Δt)
    msg *= Peridynamics.msg_qty("time step safety factor", fv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", fv.end_time)
    msg = "VELOCITY VERLET TIME STRUCTURE SOLVER\n"
    msg *= Peridynamics.msg_qty("number of time steps", vv.n_steps)
    msg *= Peridynamics.msg_qty("time step size", vv.Δt)
    msg *= Peridynamics.msg_qty("time step safety factor", vv.safety_factor)
    msg *= Peridynamics.msg_qty("simulation time", vv.end_time)
    Peridynamics.log_it(options, msg)
    return nothing
end

