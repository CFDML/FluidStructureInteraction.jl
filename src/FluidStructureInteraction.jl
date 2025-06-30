module FluidStructureInteraction

using Peridynamics
using KitBase
using Base.Threads
using TimerOutputs
using ProgressMeter
using PointNeighbors
using LinearAlgebra
using DelimitedFiles
using Glob
using WriteVTK
using Dates

abstract type AbstractFlowTimeSolver end

# Pre processing
export Post2D, Post3D

# Coupling with DSMC_sparta
export glob, write_sprata_bc_file_2d, read_bc_from_sparta, write_sparta_files, writedlm, readdlm

# New Material models
export BBTMaterial


# New Discretization
export Peridynamics, KitBase, hsource_bc!, hsource_databc!, temperature_ic!, temperature_bc!, temperature_databc!, second_bcs!, find_sec_bcs_points 

# Running simulations
export Thermstep, Flowstep, Xjob, Job, cou_submit, x_submit, IBM2D, Bcstruct

include("IBM/post_2d.jl")
include("IBM/post_3d.jl")
# include("coupling_solvers/pre_dsmc.jl")

include("structure/physics/modify.jl")
include("structure/physics/new_boundaries.jl")
include("structure/time_solvers/thermstep.jl")
include("structure/physics/bond_based_thermal_diffusion.jl")
#include("structure/core/new_submit.jl")


# include("IBM/ibm2d.jl")
include("IBM/flow2pd.jl")
include("IBM/ibm_2d_struct_flow2pd.jl")
include("IBM/ibm_2d_struct_pd2flow.jl")
include("IBM/ibm_2d.jl")
include("IBM/ibm_2d_flow2pd.jl")
include("IBM/ibm_2d_pd2flow.jl")
include("IBM/similarity_transform.jl")

include("fluid/FlowTimesolver.jl")
include("fluid/Evolution.jl")
include("fluid/Advance.jl")
include("fluid/output_vtk.jl")

#coupling_solvers
include("coupling_solvers/logs.jl") 
include("coupling_solvers/x_job.jl")
include("coupling_solvers/x_submit.jl") # central of jobs
end




