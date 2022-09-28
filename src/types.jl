
abstract type MultiPhaseSystem <: JutulSystem end
abstract type MultiComponentSystem <: MultiPhaseSystem end
const DarcyFlowModel = SimulationModel{<:Any, <:MultiPhaseSystem, <:Any, <:Any}

abstract type CompositionalSystem <: MultiComponentSystem end
const CompositionalModel = SimulationModel{D, S, F, C} where {D, S<:CompositionalSystem, F, C}

abstract type BlackOilSystem <: MultiComponentSystem end

abstract type PhaseVariables <: GroupedVariables end
abstract type ComponentVariable <: GroupedVariables end

struct MultiPhaseCompositionalSystemLV{E, T, O, R} <: CompositionalSystem where T<:Tuple
    phases::T
    components::Vector{String}
    equation_of_state::E
    rho_ref::R
    function MultiPhaseCompositionalSystemLV(equation_of_state, phases = (LiquidPhase(), VaporPhase()); reference_densities = ones(length(phases)), other_name = "Water")
        c = copy(equation_of_state.mixture.component_names)
        phases = tuple(phases...)
        T = typeof(phases)
        nph = length(phases)
        @assert nph == 2 || nph == 3
        reference_densities = tuple(reference_densities...)
        @assert length(reference_densities) == nph
        if nph == 3
            other = only(filter(x -> !(isa(x, LiquidPhase) || isa(x, VaporPhase)), phases))
            O = typeof(other)
            push!(c, other_name)
        else
            O = Nothing
        end
        only(findall(isequal(LiquidPhase()), phases))
        only(findall(isequal(VaporPhase()), phases))
        new{typeof(equation_of_state), T, O, typeof(reference_densities)}(phases, c, equation_of_state, reference_densities)
    end
end
const LVCompositionalModel = SimulationModel{D, S, F, C} where {D, S<:MultiPhaseCompositionalSystemLV{<:Any, <:Any, <:Any}, F, C}

export StandardBlackOilSystem
struct StandardBlackOilSystem{D, W, R, F, T, P} <: BlackOilSystem
    saturation_table::D
    rho_ref::R
    phase_indices::T
    phases::P
    function StandardBlackOilSystem(sat; phases = (AqueousPhase(), LiquidPhase(), VaporPhase()), reference_densities = [786.507, 1037.84, 0.969758], formulation::Symbol = :varswitch)
        phases = tuple(phases...)
        nph = length(phases)
        if nph == 2 && length(reference_densities) == 3
            reference_densities = reference_densities[2:3]
        end
        reference_densities = tuple(reference_densities...)
        @assert LiquidPhase() in phases
        @assert VaporPhase() in phases
        @assert nph == 2 || nph == 3
        @assert length(reference_densities) == nph
        phase_ind = zeros(Int64, nph)
        has_water = nph == 3
        if has_water
            phase_ind[1] = findfirst(isequal(AqueousPhase()), phases)
            offset = 1
        else
            offset = 0
        end
        phase_ind[1 + offset] = findfirst(isequal(LiquidPhase()), phases)
        phase_ind[2 + offset] = findfirst(isequal(VaporPhase()), phases)
        phase_ind = tuple(phase_ind...)
        @assert formulation == :varswitch || formulation == :zg
        new{typeof(sat), has_water, typeof(reference_densities), formulation, typeof(phase_ind), typeof(phases)}(sat, reference_densities, phase_ind, phases)
    end
end

function Base.show(io::IO, d::StandardBlackOilSystem)
    print(io, "StandardBlackOilSystem with $(d.phases)")
end


const BlackOilVariableSwitchingSystem = StandardBlackOilSystem{<:Any, <:Any, <:Any, :varswitch, <:Any, <:Any}
const BlackOilGasFractionSystem = StandardBlackOilSystem{<:Any, <:Any, <:Any, :zg, <:Any, <:Any}

const BlackOilModelVariableSwitching = SimulationModel{<:Any, BlackOilVariableSwitchingSystem, <:Any, <:Any}
const BlackOilModelGasFraction       = SimulationModel{<:Any, BlackOilGasFractionSystem,        <:Any, <:Any}
const StandardBlackOilModel          = SimulationModel{<:Any, <:StandardBlackOilSystem, <:Any, <:Any}
const StandardBlackOilModelWithWater = SimulationModel{<:Any, <:StandardBlackOilSystem{<:Any, true, <:Any, <:Any, <:Any, <:Any}, <:Any, <:Any}

struct ImmiscibleSystem{T, F} <: MultiPhaseSystem where {T<:Tuple, F<:NTuple}
    phases::T
    rho_ref::F
    function ImmiscibleSystem(phases; reference_densities = ones(length(phases)))
        phases = tuple(phases...)
        reference_densities = tuple(reference_densities...)
        new{typeof(phases), typeof(reference_densities)}(phases, reference_densities)
    end
end
Base.show(io::IO, t::ImmiscibleSystem) = print(io, "ImmiscibleSystem with $(join([typeof(p) for p in t.phases], ", "))")


struct SinglePhaseSystem{P, F} <: MultiPhaseSystem where {P, F<:AbstractFloat}
    phase::P
    rho_ref::F
    function SinglePhaseSystem(phase = LiquidPhase(); reference_density = 1.0)
        return new{typeof(phase), typeof(reference_density)}(phase, reference_density)
    end
end