const SUBMIT_KWARGS = (:quiet,)

function x_submit(job::Xjob, bcst::Bcstruct, geo::Post2D, mode::String; output=nothing, output_vars = nothing, ref = nothing, kwargs...)
    o = Dict{Symbol,Any}(kwargs)
    Peridynamics.check_kwargs(o, SUBMIT_KWARGS)
    quiet = Peridynamics.get_submit_options(o)
    Peridynamics.set_quiet!(quiet)
    
    if mode == "T" # thermal diffusion in structure
        if Peridynamics.mpi_run()
            ret = submit_mpi_therm(job)
        else 
            ret = submit_threads_therm(job, bcst.hsource, geo, nthreads(); output = output, output_vars = output_vars, ref = ref)
        end
        return ret
    else
        error("Other physical field is actively under development")
    end
end

function submit_threads_therm(job::Xjob, hsource_bc::Matrix{Float64}, geo::Post2D, n_chunks::Int; output = nothing, output_vars = nothing, ref = nothing)

    simulation_duration = @elapsed begin
        logo_init_logs(job.options)
        flow_log_spatial_setup(job.options, job.flow_setup)
        Peridynamics.log_spatial_setup(job.options, job.spatial_setup)
        Peridynamics.log_create_data_handler_start()
        dh = Peridynamics.threads_data_handler(job.spatial_setup, job.s_time_solver, n_chunks)
        Peridynamics.init_time_solver!(job.s_time_solver, dh)
        Peridynamics.initialize!(dh, job.s_time_solver)
        Peridynamics.log_create_data_handler_end()
        Peridynamics.log_data_handler(job.options, dh)
        fsi_log_timesolver(job.options, job.f_time_solver, job.s_time_solver)
        solve_therm!(dh, job, geo, job.options, hsource_bc; output = output, output_vars = output_vars, ref = ref)
    end
    Peridynamics.log_simulation_duration(job.options, simulation_duration)
    return dh
end

function solve_therm!(dh::Peridynamics.AbstractDataHandler, job::Xjob, geo::Post2D, options::Peridynamics.AbstractJobOptions, hsource_bc::Matrix{Float64}; output = nothing, output_vars = nothing, ref = nothing)
    ks = job.flow_setup
    Δt = min(job.s_time_solver.Δt, job.f_time_solver.Δt)
    sys_n = job.s_time_solver.n_steps
    
    ctr, a1face, a2face = init_fvm(job.flow_setup; structarray = true) 

    ps = job.flow_setup.ps
    if isnothing(ref)
        ib = IBM2D(ps, geo.pos, geo.bc_edge, ReferenceVariables(;length = 1.0, temperature = 1.0, velocity = 1.0, density = 1.0, pressure = 1.0))
    else
        ib = IBM2D(ps, geo.pos, geo.bc_edge, ref)
    end

    t_factors = modify_t(dh)
    conv, radi = find_sec_bcs_points(dh)    
    
    pro = Progress(sys_n; dt=1, desc="solve...", color=:normal, barlen=20,
                 enabled=Peridynamics.progress_bars())    
    for n in 1:sys_n
        cur_t = n * Δt   
        timestep_pd!(dh, options, Δt, n, t_factors, conv, radi)

        pd_out = merge_chunk(dh)
        temperature = pd_out[size(pd_out, 1)-1,:] .+ job.spatial_setup.point_params[1].rft #[pd_out[5, hull_indices[iter]] for iter in bi2pd]  
        update_ghost!(ctr, ps, ks.gas, ib, temperature)

        ibm_step!(job.flow_setup, ib, ctr, a1face, a2face, n, Δt/ib.rv.time, options)

        update_hsource_on_edge!(ps, ctr, geo.pos, geo.area, temperature, hsource_bc, ib)

        output_to_main!(Base.@locals, options, n, output, output_vars)

        next!(pro)     
    end
    finish!(pro)

end

function x_submit(job::Job, mode::String; kwargs...)
    o = Dict{Symbol,Any}(kwargs)
    Peridynamics.check_kwargs(o, SUBMIT_KWARGS)
    quiet = Peridynamics.get_submit_options(o)
    Peridynamics.set_quiet!(quiet)
    
    if mode == "T" # thermal diffusion in structure
        if Peridynamics.mpi_run()
            ret = T_submit_mpi(job)
        else 
            ret = T_submit_threads(job, nthreads())
        end
        return ret
    else
        error("Other physical field is actively under development")
    end
end

function timestep_pd!(dh::Peridynamics.AbstractThreadsBodyDataHandler, options::Peridynamics.AbstractJobOptions,
                          Δt::Float64, n::Int, t_factors::Vector{Vector{Float64}},
                          conv::Vector{Vector{Int}}, radi::Vector{Vector{Int}})
    t = n * Δt
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        Peridynamics.apply_boundary_conditions!(chunk, t)
        second_bcs!(chunk, conv[chunk_id], radi[chunk_id])
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_loc_to_halo!(dh, chunk_id)
        calc_pflux!(dh.chunks[chunk_id], t_factors[chunk_id]) 
    end

    @threads :static for chunk_id in eachindex(dh.chunks)
        Peridynamics.exchange_halo_to_loc!(dh, chunk_id)
        chunk = dh.chunks[chunk_id]        
        update_temperature!(chunk, Δt)
        Peridynamics.apply_boundary_conditions!(chunk, t)
        Peridynamics.export_results(dh, options, chunk_id, n, t)
    end
    #
    return nothing
end

function merge_chunk(C::Peridynamics.AbstractThreadsBodyDataHandler)
    n = C.chunks[end].system.chunk_handler.loc_points[end]
    pd_out = fill(0.0, 8, n)
    for i in 1:C.n_chunks
        for j in C.chunks[i].system.chunk_handler.loc_points
            m = C.chunks[i].system.chunk_handler.localizer[j]
            pd_out[1:3, j] += C.chunks[i].storage.position[1:3, m]
            pd_out[4:6, j] += C.chunks[i].storage.velocity[1:3, m]
            pd_out[7, j] += C.chunks[i].storage.temperature[1, m]
            pd_out[8, j] += C.chunks[i].storage.damage[m]
        end
    end
    return pd_out
end
