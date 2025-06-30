struct BBTMaterial{Correction} <: Peridynamics.AbstractBondSystemMaterial{Correction} end

BBTMaterial() = BBTMaterial{NoCorrection}()

struct BBTPointParameters <: Peridynamics.AbstractPointParameters
    δ::Float64
    rho::Float64
    E::Float64
    nu::Float64
    G::Float64
    K::Float64
    λ::Float64
    μ::Float64
    Gc::Float64
    εc::Float64
    bc::Float64
    kc::Float64 # thermal conductivity
    kp::Float64 # microconductivity
    aph::Float64 # thermal expansion
    cv::Float64 # specific heat capacity
    rft::Float64 # Reference temperature
    h::Float64 # convective heat transfer coefficient,
    hσ::Float64 # Stefan-Boltzman constant,
    hϵ::Float64 # emissivity
    tem∞::Float64 # temperature of the surrounding medium
end

@inline function Peridynamics.allowed_material_kwargs(::BBTMaterial)
    return (Peridynamics.DEFAULT_POINT_KWARGS..., :kc, :aph, :cv, :rft, :h, :hσ, :hϵ, :tem∞, :D)
end

function BBTPointParameters(::BBTMaterial, p::Dict{Symbol,Any})
    δ = Peridynamics.get_horizon(p)
    rho = Peridynamics.get_density(p)
    if haskey(p, :nu)
        msg = "keyword `nu` is not allowed for BBTMaterial!\n"
        msg *= "Bond-based peridynamics has a limitation on the possion ratio.\n"
        msg *= "Therefore, when using BBTMaterial, `nu` is hardcoded as 1/3 for 2D && 1/4 for 3D.\n"
        throw(ArgumentError(msg))
    end
    
    if haskey(p, :D)
        p[:nu] = 1/3
    else 
        p[:nu] = 1/4    
    end
   
    E, nu, G, K, λ, μ = Peridynamics.get_elastic_params(p)
    Gc, εc = Peridynamics.get_frac_params(p, δ, K)
    kc, aph, cv, rft, h, hσ, hϵ, tem∞ = get_thermal_params(p, δ)

    if haskey(p, :D) #2D 
        bc = 9 * E / (π * 1 * δ^3) # bond constant
        kp = 6 * kc / (π * 1 * δ^3) # microcndicitvity constant
    else #3D
        bc = 12 * E / (π * δ^4) 
        kp = 6 * kc / (π * δ^4) 
    end   
    return BBTPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc, kc, kp, aph, cv, rft, h, hσ, hϵ, tem∞)
end

function get_thermal_params(p::Dict{Symbol,Any}, δ)
    haskey(p, :kc) || throw(UndefKeywordError(:kc))
    haskey(p, :aph) || throw(UndefKeywordError(:aph))
    haskey(p, :cv) || throw(UndefKeywordError(:cv))
    haskey(p, :rft) || throw(UndefKeywordError(:rft))
    haskey(p, :h) || throw(UndefKeywordError(:h))  
    haskey(p, :hσ) || throw(UndefKeywordError(:hσ))  
    haskey(p, :hϵ) || throw(UndefKeywordError(:hϵ))  
    haskey(p, :tem∞) || throw(UndefKeywordError(:tem∞))  

    kc::Float64 = float(p[:kc])
    kc ≤ 0 && throw(ArgumentError("`kc` should be larger than zero!\n"))
    aph::Float64 = float(p[:aph])
    aph ≤ 0 && throw(ArgumentError("`aph` should be larger than zero!\n"))
    cv::Float64 = float(p[:cv])
    cv ≤ 0 && throw(ArgumentError("`cv` should be larger than zero!\n"))
    rft::Float64 = float(p[:rft])
    h::Float64 = float(p[:h])
    h ≤ 0 && throw(ArgumentError("`h` should be larger than zero!\n"))
    hσ::Float64 = float(p[:hσ])
    hσ ≤ 0 && throw(ArgumentError("`hσ` should be larger than zero!\n"))
    hϵ::Float64 = float(p[:hϵ])
    hϵ ≤ 0 && throw(ArgumentError("`hϵ` should be larger than zero!\n"))
    tem∞::Float64 = float(p[:tem∞])
        
    return kc, aph, cv, rft, h, hσ, hϵ, tem∞
end

@Peridynamics.params BBTMaterial BBTPointParameters


# BBT&therm store
struct BBTthermtStorage <: Peridynamics.AbstractStorage
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    acceleration::Matrix{Float64}
    b_int::Matrix{Float64}
    b_ext::Matrix{Float64}
    damage::Vector{Float64}
    bond_active::Vector{Bool}
    n_active_bonds::Vector{Int}
    temperature::Matrix{Float64}
    pflux::Matrix{Float64}
    hsource::Matrix{Float64}
end

function BBTthermtStorage(::BBTMaterial, ::Thermstep,
     system::Peridynamics.BondSystem, ch)
    n_loc_points = length(ch.loc_points)
    n_halo = length(ch.halo_points)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, n_loc_points)
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    temperature = zeros(1, n_loc_points+n_halo)
    pflux = zeros(1, n_loc_points)
    hsource = zeros(1, n_loc_points)
    s = BBTthermtStorage(position, displacement, velocity, velocity_half, acceleration,
                        b_int, b_ext, damage, bond_active, n_active_bonds, temperature, pflux, hsource)
    return s
end

@Peridynamics.storage BBTMaterial Thermstep BBTthermtStorage

Peridynamics.point_data_field(s::BBTthermtStorage, ::Val{:temperature}) = getfield(s, :temperature)
Peridynamics.point_data_field(s::BBTthermtStorage, ::Val{:pflux}) = getfield(s, :pflux)
Peridynamics.point_data_field(s::BBTthermtStorage, ::Val{:hsource}) = getfield(s, :hsource)

@Peridynamics.loc_to_halo_fields BBTthermtStorage :temperature

const BBTStorage = Union{BBTthermtStorage}

# Thermal 
@inline function update_temperature!(b::Peridynamics.AbstractBodyChunk, Δt::Float64)
    for point_id in Peridynamics.each_point_idx(b)
        param = Peridynamics.get_params(b, point_id)
        k = Δt / (param.rho * param.cv)
        _update_temperature!(b.storage.temperature, b.storage.pflux, b.storage.hsource, k, point_id)
    end
    return nothing
end

@inline function _update_temperature!(temperature, pflux, hsource, k, i)
    temperature[1, i] += (pflux[1, i] + hsource[1, i]) * k
    return nothing
end

function calc_pflux!(chunk::Peridynamics.AbstractBodyChunk, mbd_t::Vector{Float64})
    (; system, mat, paramsetup, storage) = chunk
    storage.pflux .= 0.0
    storage.n_active_bonds .= 0
    for point_id in eachindex(chunk.ch.loc_points)
        pflux_point!(storage, system, mat, paramsetup, point_id, mbd_t)
    end
    return nothing
end

function pflux_point!(storage::BBTStorage, system::Peridynamics.BondSystem, 
    ::BBTMaterial, param::BBTPointParameters, i::Int, mbd_t::Vector{Float64}) 

    for bond_id in system.bond_ids[i]
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        
        mof_th = mbd_t[bond_id]
        
        Δtem = storage.temperature[1, j] - storage.temperature[1, i]
        storage.pflux[1, i] += Peridynamics.bond_failure(storage, bond_id) * (param.kp * Δtem / L) * 
        system.volume[j] * mof_th
    end
    return nothing
end

