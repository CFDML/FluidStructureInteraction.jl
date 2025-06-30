
function logo_init_logs(options::Peridynamics.AbstractJobOptions)
    Peridynamics.mpi_isroot() || return nothing
    Peridynamics.print_log(fluidstructureinteraction_banner())
    Peridynamics.print_log(Peridynamics.get_run_info())
    Peridynamics.set_progress_bars!()
    options.export_allowed || return nothing
    mkpath(options.vtk)
    logo_init_logfile(options)
    return nothing
end

function logo_init_logfile(options::Peridynamics.AbstractJobOptions)
    open(options.logfile, "w+") do io
        write(io, get_logfile_head_fsi())
        write(io,fluidstructureinteraction_banner(color=false))
        write(io, Peridynamics.get_run_info())
    end
    return nothing
end

function get_logfile_head_fsi()
    msg = "LOGFILE CREATED ON "
    msg *= Dates.format(Dates.now(), "yyyy-mm-dd, HH:MM:SS")
    msg *= "\n"
    msg *= "\n"
    return msg
end

function fluidstructureinteraction_banner(; color::Bool=true, indentation::Int=10)
    indent = indentation > 0 ? " "^indentation : ""
    if color
        fluid_color = "\e[38;2;65;105;225m"
        structure_color = "\e[38;2;50;205;150m"
        interaction_color = "\e[38;2;220;20;60m"
        bold = "\e[1m"
        reset = "\e[0m"
    else
        fluid_color = structure_color = interaction_color = bold = reset = ""
    end

    msg  = indent * "$(fluid_color)    ______ __        _      __ $(reset)$(structure_color)_____  __                      __                   $(reset)\n"
    msg *= indent * "$(fluid_color)   / ____// /__  __ (_)____/ /$(reset)$(structure_color)/ ___/ / /_ _____ __  __ _____ / /_ __  __ _____ ___ $(reset)\n"
    msg *= indent * "$(fluid_color)  / /_   / // / / // // __  /$(reset)$(structure_color) \\__ \\ / __// ___// / / // ___// __// / / // ___// _ \\$(reset)\n"
    msg *= indent * "$(fluid_color) / __/  / // /_/ // // /_/ /$(reset)$(structure_color) ___/ // /_ / /   / /_/ // /__ / /_ / /_/ // /   /  __/$(reset)\n"
    msg *= indent * "$(fluid_color)/_/    /_/ \\__,_//_/ \\__,_/$(reset)$(structure_color) /____/ \\__//_/    \\__,_/ \\___/ \\__/ \\__,_//_/    \\___/ $(reset)\n"
    msg *= indent * "\n"
    msg *= indent * "$(interaction_color)         ____        __                            __   _                          $(reset)\n"
    msg *= indent * "$(interaction_color)        /  _/____   / /_ ___   _____ ____ _ _____ / /_ (_)____   ____              $(reset)\n"
    msg *= indent * "$(interaction_color)        / / / __ \\ / __// _ \\ / ___// __ `// ___// __// // __ \\ / __ \\             $(reset)\n"
    msg *= indent * "$(interaction_color)      _/ / / / / // /_ /  __// /   / /_/ // /__ / /_ / // /_/ // / / /             $(reset)\n"
    msg *= indent * "$(interaction_color)     /___//_/ /_/ \\__/ \\___//_/    \\__,_/ \\___/ \\__//_/ \\____//_/ /_/             $(reset)\n"
    
    msg *= indent * "\n"
    msg *= indent * "$(fluid_color)Fluid$(reset)$(structure_color)Structure$(reset)$(interaction_color)Interaction$(reset): Fluid-Structure Multiphysics Coupling Framework\n"
    msg *= indent * "Copyright © 2025 Shiwei Hu\n\n"
    
    return msg
end
