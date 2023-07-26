function current_unit_system(deck)
    rs = deck["RUNSPEC"]
    choices = ["METRIC", "SI", "LAB", "FIELD"]
    found = false
    sys = missing
    for choice in choices
        if haskey(rs, choice)
            @assert ismissing(sys) "Cannot have multiple unit keywords (encountered $choice and $sys)"
            sys = choice
        end
    end
    if ismissing(sys)
        @warn "Units not found. Assuming field units."
        out = :field
    else
        out = Symbol(lowercase(sys))
    end
    return out
end

Base.@kwdef struct DeckUnitSystem{S, T}
    length::T = 1.0
    area::T = 1.0
    time::T = 1.0
    density::T = 1.0
    pressure::T = 1.0
    mol::T = 1.0
    mass::T = 1.0
    concentration::T = 1.0
    compressibility::T = 1.0
    viscosity::T = 1.0
    surface_tension::T = 1.0
    jsurface_tension::T = 1.0
    permeability::T = 1.0
    liquid_volume_surface::T = 1.0
    liquid_volume_reservoir::T = 1.0
    gas_volume_surface::T = 1.0
    gas_volume_reservoir::T = 1.0
    volume::T = 1.0
    transmissibility::T = 1.0
    rock_conductivity::T = 1.0
    volume_heat_capacity::T = 1.0
    mass_heat_capacity::T = 1.0
    relative_temperature::Symbol = :Celsius
    absolute_temperature::Symbol = :Kelvin
end

function DeckUnitSystem(sys::Symbol, T = Float64)
    #meter, day, kilogram, bar = si_units(:meter, :day, :kilogram, :bar)
    u = Jutul.all_units()
    m = u[:meter]
    K = u[:kelvin]
    day = u[:day]
    centi = u[:centi]
    kilogram = u[:kilogram]
    ft = u[:feet]
    psi = u[:psi]
    pound = u[:pound]
    kilo = u[:kilo]
    stb = u[:stb]
    rankine = u[:rankine]
    btu = u[:btu]

    # Commons
    cP = u[:centi]*u[:poise]
    mD = u[:milli]*u[:darcy]
    if sys == :metric
        len = m
        kJ = u[:kilo]*u[:joule]
        volume = m^3
        area = m^2
        time = day
        pressure = u[:bar]
        mol = u[:kilo]
        mass = kilogram
        viscosity = cP
        surface_tension = u[:newton]/m
        jsurface_tension = u[:dyne]/(centi*m)
        permeability = mD
        liquid_volume_surface = volume
        liquid_volume_reservoir = volume
        gas_volume_surface = volume
        gas_volume_reservoir = volume
        rock_conductivity = kJ/(m*day*K)
        volume_heat_capacity = kJ/(volume*K)
        mass_heat_capacity = kJ/(mass*K)
        relative_temperature = :Celsius
        absolute_temperature = :Kelvin
    elseif sys == :field
        len = ft
        area = ft^2
        time = day
        pressure = psi
        mol = pound*kilo
        mass = pound
        viscosity = cP
        surface_tension = u[:lbf]/u[:inch]
        jsurface_tension = u[:dyne]/(centi*m)
        permeability = mD
        liquid_volume_surface = stb
        liquid_volume_reservoir = stb
        gas_volume_surface = kilo*ft^3
        gas_volume_reservoir = stb
        volume = ft^3
        rock_conductivity = btu / (ft*day*rankine)
        volume_heat_capacity = btu / (ft^3*rankine)
        mass_heat_capacity = btu / (pound*rankine)
        relative_temperature = :Fahrenheit
        absolute_temperature = :Rankine
    elseif sys == :lab
        error("Not implemented")
    else
        @assert sys == :si
        len = 1.0
        area = 1.0
        time = 1.0
        pressure = 1.0
        mol = 1.0
        mass = 1.0
        viscosity = 1.0
        surface_tension = 1.0
        jsurface_tension = 1.0
        permeability = 1.0
        liquid_volume_surface = 1.0
        liquid_volume_reservoir = 1.0
        gas_volume_surface = 1.0
        gas_volume_reservoir = 1.0
        volume = 1.0
        rock_conductivity = 1.0
        volume_heat_capacity = 1.0
        mass_heat_capacity = 1.0
        relative_temperature = :Celsius
        absolute_temperature = :Kelvin
    end
    density = mass/volume
    concentration = mass/volume
    compressibility = 1.0/pressure
    transmissibility = viscosity * volume / (time * pressure)
    return DeckUnitSystem{sys, T}(
        length = len,
        area = area,
        time = time,
        density = density,
        pressure = pressure,
        mol = mol,
        mass = mass,
        concentration = concentration,
        compressibility = compressibility,
        viscosity = viscosity,
        surface_tension = surface_tension,
        jsurface_tension = jsurface_tension,
        permeability = permeability,
        liquid_volume_surface = liquid_volume_surface,
        liquid_volume_reservoir = liquid_volume_reservoir,
        gas_volume_surface = gas_volume_surface,
        gas_volume_reservoir = gas_volume_reservoir,
        volume = volume,
        transmissibility = transmissibility,
        rock_conductivity = rock_conductivity,
        volume_heat_capacity = volume_heat_capacity,
        mass_heat_capacity = mass_heat_capacity,
        relative_temperature = relative_temperature,
        absolute_temperature = absolute_temperature
    )
end

function swap_unit_system!(x::AbstractArray, systems, k)
    for i in eachindex(x)
        x[i] = swap_unit_system(x[i], systems, k)
    end
    return x
end

function swap_unit_system(v, systems::Nothing, k)
    return v
end

function swap_unit_system(val, systems, k)
    return swap_unit_system(val, systems, Val(k))
end


function swap_unit_system(val, systems::NamedTuple, ::Val{k}) where k
    (; to, from) = systems
    to_unit = getproperty(to, k)
    from_unit = getproperty(from, k)
    val_si = convert_to_si(val, from_unit)
    val_final = convert_from_si(val_si, to_unit)
    return val_final
end