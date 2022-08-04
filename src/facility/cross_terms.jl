struct ReservoirFromWellCT{T<:AbstractVector, I<:AbstractVector} <: Jutul.AdditiveCrossTerm
    WI::T
    reservoir_cells::I
    well_cells::I
end

Jutul.symmetry(::ReservoirFromWellCT) = Jutul.CTSkewSymmetry()

function update_cross_term_in_entity!(out, i,
    state_t, state0_t,
    state_s, state0_s, 
    model_t, model_s,
    param_t, param_s,
    ct::ReservoirFromWellCT, eq, dt, ldisc = local_discretization(ct, i))
    # Unpack properties
    reservoir_cell = ct.reservoir_cells[i]
    well_cell = ct.well_cells[i]
    WI = ct.WI[i]
    rhoS = param_s[:reference_densities]

    p_well = state_s.Pressure
    p_res = state_t.Pressure
    # Todo: Fix conn -> cell pressure drop
    ρgdz = 0
    dp = -WI*(p_well[well_cell] - p_res[reservoir_cell] + ρgdz)
    # Call smaller interface that is easy to specialize
    well_perforation_flux!(out, model_t.system, state_t, state_s, rhoS, dp, reservoir_cell, well_cell)
end

Jutul.cross_term_entities(ct::ReservoirFromWellCT, eq::ConservationLaw, model) = ct.reservoir_cells
Jutul.cross_term_entities_source(ct::ReservoirFromWellCT, eq::ConservationLaw, model) = ct.well_cells

# Well influence on facility
struct FacilityFromWellCT <: Jutul.AdditiveCrossTerm
    well::Symbol
end

well_top_node() = 1

Jutul.cross_term_entities(ct::FacilityFromWellCT, eq::ControlEquationWell, model) = get_well_position(model.domain, ct.well)

function Jutul.prepare_cross_term_in_entity!(i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    param_f, param_w,
    ct::FacilityFromWellCT, eq, dt, ldisc = local_discretization(ct, i))
    well_symbol = ct.well
    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)

    target = ctrl.target
    q_t = facility_surface_mass_rate_for_well(facility, well_symbol, state_facility)
    if !isa(target, DisabledTarget)
        cfg = state_facility.WellGroupConfiguration
        limits = current_limits(cfg, well_symbol)
        if !isnothing(limits)
            rhoS = param_w[:reference_densities]
            rhoS, S = flash_wellstream_at_surface(well, state_well, rhoS)
            rhoS = tuple(rhoS...)
            apply_well_limit!(cfg, target, well, state_well, well_symbol, rhoS, S, value(q_t), limits)    
        end
    end
end

function update_cross_term_in_entity!(out, i,
    state_facility, state0_facility,
    state_well, state0_well,
    facility, well,
    param_f, param_w,
    ct::FacilityFromWellCT, eq, dt, ldisc = local_discretization(ct, i))

    well_symbol = ct.well
    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)

    target = ctrl.target
    q_t = facility_surface_mass_rate_for_well(facility, well_symbol, state_facility)
    if isa(target, DisabledTarget)
        # Early return - no cross term needed.
        t = q_t + 0*bottom_hole_pressure(state_well)
        t_num = 0.0
    else
        cfg = state_facility.WellGroupConfiguration
        need_rates = isa(ctrl, ProducerControl) && !isa(target, BottomHolePressureTarget)
        rhoS = param_w[:reference_densities]
        if need_rates
            rhoS, S = flash_wellstream_at_surface(well, state_well, rhoS)
        else
            S = nothing
        end
        rhoS = tuple(rhoS...)
        t = well_target(ctrl, target, well, state_well, rhoS, S)
        if rate_weighted(target)
            t *= q_t
        end
        t_num = target.value
    end
    t += 0*bottom_hole_pressure(state_well) + 0*q_t
    scale = target_scaling(target)
    eq = (t - t_num)/scale
    out[] = eq
end

# Facility influence on well
struct WellFromFacilityCT <: Jutul.AdditiveCrossTerm
    well::Symbol
end

Jutul.cross_term_entities(ct::WellFromFacilityCT, eq::ConservationLaw, model) = well_top_node()

function update_cross_term_in_entity!(out, i,
    state_well, state0_well,
    state_facility, state0_facility,
    well, facility,
    param_w, param_f,
    ct::WellFromFacilityCT, eq, dt, ldisc = local_discretization(ct, i))

    well_symbol = ct.well
    pos = get_well_position(facility.domain, well_symbol)

    cfg = state_facility.WellGroupConfiguration
    ctrl = operating_control(cfg, well_symbol)
    qT = state_facility.TotalSurfaceMassRate[pos] 
    # Hack for sparsity detection
    qT += 0*bottom_hole_pressure(state_well)

    if isa(ctrl, InjectorControl)
        if value(qT) < 0
            @warn "Injector $well_symbol is producing?"
        end
        mix = ctrl.injection_mixture
        nmix = length(mix)
        ncomp = number_of_components(well.system)
        @assert nmix == ncomp "Injection composition length ($nmix) must match number of components ($ncomp)."
    else
        if value(qT) > 0
            @warn "Producer $well_symbol is injecting?"
        end
        masses = state_well.TotalMasses[:, well_top_node()]
        mass = sum(masses)
        mix = masses./mass
    end
    for i in eachindex(out)
        out[i] = -mix[i]*qT
    end
end

