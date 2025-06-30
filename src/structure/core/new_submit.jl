const SUBMIT_KWARGS = (:quiet,)

function cou_submit(job::Job, mode::String; kwargs...)
    o = Dict{Symbol,Any}(kwargs)
    Peridynamics.check_kwargs(o, SUBMIT_KWARGS)
    quiet = Peridynamics.get_submit_options(o)
    Peridynamics.set_quiet!(quiet)

    if Peridynamics.mpi_run()
        ret = cou_submit_mpi(job, mode)
    else
        ret = cou_submit_threads(job, mode, nthreads())
    end
    return ret
end

function get_submit_options(o::Dict{Symbol,Any})
    local quiet::Bool
    if haskey(o, :quiet)
        quiet = Bool(o[:quiet])
    else
        quiet = false
    end
    return quiet
end

function cou_submit_threads(job::Job, mode::String, n_chunks::Int)
    simulation_duration = @elapsed begin
        Peridynamics.init_logs(job.options)
        Peridynamics.log_spatial_setup(job.options, job.spatial_setup)
        Peridynamics.log_create_data_handler_start()
        dh = Peridynamics.threads_data_handler(job.spatial_setup, job.time_solver, n_chunks)
        Peridynamics.init_time_solver!(job.time_solver, dh)
        Peridynamics.initialize!(dh, job.time_solver)
        Peridynamics.log_create_data_handler_end()
        Peridynamics.log_data_handler(job.options, dh)
        Peridynamics.log_timesolver(job.options, job.time_solver)
        if mode == "cou"
            cou_solve!(dh, job.time_solver, job.options)
        elseif mode == "therm"
            th_solve!(dh, job.time_solver, job.options)
        else
            Peridynamics.solve!(dh, job.time_solver, job.options)
        end
    end
    Peridynamics.log_simulation_duration(job.options, simulation_duration)
    return dh
end


