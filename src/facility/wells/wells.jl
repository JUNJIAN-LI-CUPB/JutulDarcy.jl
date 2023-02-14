export WellGrid, MultiSegmentWell
export TotalMassFlux, PotentialDropBalanceWell, SegmentWellBoreFrictionHB

export InjectorControl, ProducerControl, SinglePhaseRateTarget, BottomHolePressureTarget

export Perforations
export MixedWellSegmentFlow
export segment_pressure_drop


abstract type WellPotentialFlowDiscretization <: PotentialFlowDiscretization end

"""
Two point approximation with flux for wells
"""
struct MixedWellSegmentFlow <: WellPotentialFlowDiscretization end

abstract type WellGrid <: PorousMediumGrid
    # Wells are not porous themselves per se, but they are discretizing
    # part of a porous medium.
end

function Base.show(io::IO, w::WellGrid)
    if w isa SimpleWell
        nseg = 0
        nn = 1
        n = "SimpleWell"
    else
        nseg = size(w.neighborship, 2)
        nn = length(w.volumes)
        n = "MultiSegmentWell"
    end
    print(io, "$n [$(w.name)] ($(nn) nodes, $(nseg) segments, $(length(w.perforations.reservoir)) perforations)")
end

# Total velocity in each well segment
struct TotalMassFlux <: ScalarVariable
    scale
    max_abs
    max_rel
    function TotalMassFlux(;scale = 3600*24, max_abs = nothing, max_rel = nothing)
        new(scale, max_abs, max_rel)
    end
end

relative_increment_limit(tmf::TotalMassFlux) = tmf.max_rel
absolute_increment_limit(tmf::TotalMassFlux) = tmf.max_abs

associated_entity(::TotalMassFlux) = Faces()
Jutul.variable_scale(t::TotalMassFlux) = t.scale

default_surface_cond() = (p = 101325.0, T = 288.15) # Pa and deg. K from ISO 13443:1996 for natural gas

struct SimpleWell{SC, P} <: WellGrid where {SC, P}
    perforations::P
    surface::SC
    name::Symbol
end

function SimpleWell(reservoir_cells; name = :Well, surface_conditions = default_surface_cond(), kwarg...)
    nr = length(reservoir_cells)
    WI, gdz = common_well_setup(nr; kwarg...)
    perf = (self = ones(Int64, nr), reservoir = vec(reservoir_cells), WI = WI, gdz = gdz)
    return SimpleWell(perf, surface_conditions, name)
end
struct MultiSegmentWell{V, P, N, A, C, SC, S} <: WellGrid
    volumes::V          # One per cell
    perforations::P     # (self -> local cells, reservoir -> reservoir cells, WI -> connection factor)
    neighborship::N     # Well cell connectivity
    top::A              # "Top" node where scalar well quantities live
    centers::C          # Coordinate centers of nodes
    surface::SC         # p, T at surface
    name::Symbol        # Symbol that names the well
    segment_models::S   # Segment pressure drop model for each segment
    function MultiSegmentWell(reservoir_cells, volumes::AbstractVector, centers;
                                                        N = nothing,
                                                        name = :Well,
                                                        perforation_cells = nothing,
                                                        segment_models = nothing,
                                                        reference_depth = 0,
                                                        dz = nothing,
                                                        surface_conditions = default_surface_cond(),
                                                        accumulator_volume = mean(volumes),
                                                        kwarg...)
        nv = length(volumes)
        nc = nv + 1
        reservoir_cells = vec(reservoir_cells)
        nr = length(reservoir_cells)
        if isnothing(N)
            @debug "No connectivity. Assuming nicely ordered linear well."
            N = vcat((1:nv)', (2:nc)')
        elseif maximum(N) == nv
            N = vcat([1, 2], N+1)
        end
        nseg = size(N, 2)
        @assert size(N, 1) == 2
        @assert size(centers, 1) == 3
        volumes = vcat([accumulator_volume], volumes)
        ext_centers = hcat([centers[1:2, 1]..., reference_depth], centers)
        @assert length(volumes) == size(ext_centers, 2)
        if !isnothing(reservoir_cells) && isnothing(perforation_cells)
            @assert length(reservoir_cells) == nv "If no perforation cells are given, we must 1->1 correspondence between well volumes and reservoir cells."
            perforation_cells = collect(2:nc)
        end
        perforation_cells = vec(perforation_cells)

        if isnothing(segment_models)
            Δp = SegmentWellBoreFrictionHB(1.0, 1e-4, 0.1)
            segment_models = repeat([Δp], nseg)
        else
            segment_models::AbstractVector
            @assert length(segment_models) == nseg
        end
        if isnothing(dz)
            dz = centers[3, :] - reference_depth
        end
        @assert length(perforation_cells) == nr
        WI, gdz = common_well_setup(nr; dz = dz, kwarg...)
        perf = (self = perforation_cells, reservoir = reservoir_cells, WI = WI, gdz = gdz)
        accumulator = (reference_depth = reference_depth, )
        new{typeof(volumes), typeof(perf), typeof(N), typeof(accumulator), typeof(ext_centers), typeof(surface_conditions), typeof(segment_models)}(volumes, perf, N, accumulator, ext_centers, surface_conditions, name, segment_models)
    end
end

function common_well_setup(nr; dz = nothing, WI = nothing, gravity = gravity_constant)
    if isnothing(dz)
        @warn "dz not provided for well. Assuming no gravity."
        gdz = zeros(nr)
    else
        @assert length(dz) == nr  "Must have one connection drop dz per perforated cell"
        gdz = dz*gravity
    end
    if isnothing(WI)
        @warn "No well indices provided. Using 1e-12."
        WI = repeat(1e-12, nr)
    else
        @assert length(WI) == nr  "Must have one well index per perforated cell"
    end
    return (WI, gdz)
end

export setup_well, setup_vertical_well
function setup_well(g, K, reservoir_cells::AbstractVector;
                                        reference_depth = nothing, 
                                        geometry = tpfv_geometry(g),
                                        skin = 0.0,
                                        Kh = nothing,
                                        radius = 0.1,
                                        simple_well = false,
                                        dir = :z,
                                        kwarg...)
    n = length(reservoir_cells)
    # Make sure these are cell indices
    reservoir_cells = map(i -> cell_index(g, i), reservoir_cells)
    centers = geometry.cell_centroids[:, reservoir_cells]
    if isnothing(reference_depth)
        reference_depth = centers[3, 1]
    end
    volumes = zeros(n)
    WI = zeros(n)
    dz = zeros(n)
    for (i, c) in enumerate(reservoir_cells)
        if K isa AbstractVector
            k = K[c]
        else
            k = K[:, c]
        end
        WI[i] = compute_peaceman_index(g, k, radius, c, skin = skin, Kh = Kh)
        center = vec(centers[:, i])
        dz[i] = center[3] - reference_depth
        if dir isa Symbol
            d = dir
        else
            d = dir[i]
        end
        Δ = cell_dims(g, c)
        d_index = findfirst(isequal(d), [:x, :y, :z])
        h = Δ[d_index]
        volumes[i] = h*π*radius^2
    end
    if simple_well
        W = SimpleWell(reservoir_cells, WI = WI, dz = dz, reference_depth = reference_depth, kwarg...)
    else
        W = MultiSegmentWell(reservoir_cells, volumes, centers; WI = WI, dz = dz, reference_depth = reference_depth, kwarg...)
    end
    return W
end

function Jutul.plot_primitives(mesh::MultiSegmentWell, plot_type; kwarg...)
    # By default, no plotting is supported
    if plot_type == :lines
        centers = mesh.centers
        if ndims(centers) == 3
            # Some bug somewhere in MRST parser.
            centers = dropdims(centers, dims = 3)
        end
        pts = collect(centers')
        top = [pts[1, 1] pts[1, 2] pts[1, 3] - 10.0]
        pts = vcat(top, pts)
        @. pts[:, 3] *= -1

        function cell_mapper(x::AbstractVector)
            return vcat(x[1], x)
        end

        function cell_mapper(x::AbstractMatrix)
            return hcat(x[:, 1], x)
        end
        mapper = (Cells = x -> cell_mapper(x), )
        out = (points = pts, mapper = mapper, top_text = String(mesh.name), marker_size = 20)
    else
        out = nothing
    end
    return out
end

function setup_vertical_well(g, K, i, j; heel = 1, toe = grid_dims_ijk(g)[3], kwarg...)
    @assert heel <= toe
    @assert heel > 0
    @assert toe > 0
    k_range = heel:toe
    n = length(k_range)
    @assert n > 0
    reservoir_cells = zeros(Int64, n)
    for (ix, k) in enumerate(k_range)
        reservoir_cells[ix] = cell_index(g, (i, j, k))
    end
    return setup_well(g, K, reservoir_cells; kwarg...)
end

include("mswells_equations.jl")



function fluid_volume(grid::WellGrid)
    return grid.volumes
end

# Well segments

function get_neighborship(::SimpleWell)
    # No interior connections.
    return zeros(Int64, 2, 0)
end

function number_of_cells(W::SimpleWell)
    return 1
end

function number_of_cells(W::WellGrid)
    return length(W.volumes)
end

function declare_entities(W::WellGrid)
    c = (entity = Cells(),         count = number_of_cells(W))
    f = (entity = Faces(),         count = number_of_faces(W))
    p = (entity = Perforations(),  count = length(W.perforations.self))
    return [c, f, p]
end

Base.@propagate_inbounds function well_perforation_flux!(out, w, sys::Union{ImmiscibleSystem, SinglePhaseSystem}, state_res, state_well, rhoS, conn)
    rc = conn.reservoir
    wc = conn.well
    # Reservoir quantities
    ρ = state_res.PhaseMassDensities
    # Extra mobility needed
    kr = state_res.RelativePermeabilities
    μ = state_res.PhaseViscosities
    is_ms = w isa MultiSegmentWell
    nph = size(ρ, 1)
    if is_ms
        # Well quantities
        ρ_w = state_well.PhaseMassDensities
        # Saturation instead of mobility - use total mobility form
        s_w = state_well.Saturations
        λ_t = 0
        for ph in 1:nph
            λ_t += kr[ph, rc]/μ[ph, rc]
        end
    else
        ρλ_t = 0
        for ph in 1:nph
            ρλ_t += ρ[ph, rc]*kr[ph, rc]/μ[ph, rc]
        end
        X = state_well.MassFractions
    end
    for ph in 1:nph
        # ψ is pressure difference from reservoir to well. If it is negative, we are injecting into the reservoir.
        ψ = perforation_phase_potential_difference(conn, state_res, state_well, ph)
        if ψ < 0
            # Injection
            if is_ms
                out[ph] = s_w[ph, wc]*ρ_w[ph, wc]*ψ*λ_t
            else
                out[ph] = X[ph]*ψ*ρλ_t
            end
        else
            # Production
            λ = kr[ph, rc]/μ[ph, rc]
            out[ph] = λ*ρ[ph, rc]*ψ
        end
    end
    return out
end

const WellDomain = DiscretizedDomain{<:WellGrid}
include("mswells.jl")
include("stdwells.jl")

include("stdwells_equations.jl")

# Some utilities
function mix_by_mass(masses, total::T, values) where T
    v = zero(T)
    @inbounds for i in eachindex(masses)
        v += masses[i]*values[i]
    end
    return v/total
end

function mix_by_saturations(s, values)
    v = zero(eltype(s))
    @inbounds for i in eachindex(s)
        v += s[i]*values[i]
    end
    return v
end

function mix_by_saturations(s::Real, values)
    return s*values[]
end

function setup_forces(model::SimulationModel{D, S}; mask = nothing) where {D <: DiscretizedDomain{G}, S<:MultiPhaseSystem} where G<:WellGrid
    mask::Union{Nothing, PerforationMask}
    return (mask = mask,)
end

function apply_perforation_mask!(storage::NamedTuple, mask::AbstractVector)
    function mask_row!(M::AbstractMatrix, m, ix)
        for i in axes(M, 1)
            M[i, ix] *= m
        end
    end
    function mask_row!(M::AbstractVector, m, ix)
        M[ix] *= m
    end
    for (k, s) in pairs(storage)
        if k == :numeric
            continue
        end
        v = s.entries
        for i in 1:Jutul.number_of_entities(s)
            mask_value = mask[i]
            for j in Jutul.vrange(s, i)
                mask_row!(v, mask_value, j)
            end
        end
    end
end

function flash_wellstream_at_surface(well_model, well_state, rhoS)
    fsys = flow_system(well_model.system)
    return flash_wellstream_at_surface(well_model, fsys, well_state, rhoS)
end

function flash_wellstream_at_surface(well_model, system::ImmiscibleSystem, well_state, rhoS)
    if haskey(well_state, :MassFractions)
        X = well_state.MassFractions
    else
        X = well_state.TotalMasses[:, 1]
    end
    vol = X./rhoS
    volfrac = vol./sum(vol)
    return (rhoS, volfrac)
end

function flash_wellstream_at_surface(well_model, system::SinglePhaseSystem, well_state, rhoS)
    return (rhoS, [1.0])
end