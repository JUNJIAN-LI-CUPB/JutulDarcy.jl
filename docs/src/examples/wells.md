# Example with wells
This example demonstrates how to set up a 3D domain with a layered permeability field, define wells and solve a simple production-injection schedule. We begin by loading the `Jutul` package that contains generic features like grids and linear solvers and the `JutulDarcy` package itself.

## Preliminaries
```@example intro_wells
using JutulDarcy, Jutul
```
`JutulDarcy` uses SI units internally. It is therefore convenient to define a few constants at the start of the script to have more managable numbers later on.
```@example intro_wells
bar = 1e5
day = 3600*24.0
Darcy = 9.869232667160130e-13
```
!!! info "Note about units"
    JutulDarcy does currently not make us of conversion factors or explicit units can in principle use any consistent unit system. Some default scaling of variables assume that the magnitude pressures and velocities roughly match that of strict SI (e.g. Pascals and cubic meters per second). These scaling factors are primarily used when iterative linear solvers are used.


## Defining a porous medium
We start by defining the static part of our simulation problem -- the porous medium itself.
### Defining the grid
The first step is to create a grid for our simulation domain. We make a tiny 5 by 5 grid with 4 layers that discretizes a physical domain of 2000 by 1500 by 50 meters.
```@example intro_wells
nx = ny = 5
nz = 4
dims = (nx, ny, nz)
g = CartesianMesh(dims, (2000.0, 1500.0, 50.0))
```
### Adding properties and making a domain
The grid by itself does not fully specify a porous medium. For that we need to specify the permeability in each cell and the porosity. Permeability, often denoted by a positive-definite tensor K, describes the relationship between a pressure gradient and the flow through the medium. Porosity is a dimensionless number between 0 and 1 that describes how much of the porous medium is void space where fluids can be present. The assumption of Darcy flow becomes less reasonable for high porosity values and the flow equations break down at zero porosity. A porosity of 0.2 is then a safe choice.

Jutul uses the `DataDomain` type to store a domain/grid together with data. For porous media simulation, `JutulDarcy` includes a convenience function `reservoir_domain` that contains defaults for permeability and porosity. We specify the permeability per-cell with varying values per layer in the vertical direction and a single porosity value for all cells that the function will expand for us. From the output, we can see that basic geometry primitives are also automatically added:
```@example intro_wells
nlayer = nx*ny # Cells in each layer
K = vcat(
    repeat([0.65], nlayer),
    repeat([0.3], nlayer),
    repeat([0.5], nlayer),
    repeat([0.2], nlayer)
    )*Darcy
# Set up domain with given permeability and porosity
domain = reservoir_domain(g, permeability = K, porosity = 0.2)
```

## Defining wells
Now that we have a porous medium with all static properties set up, it is time to introduce some driving forces. Jutul assumes no-flow boundary conditions on all boundary faces unless otherwise specified so we can go ahead and add wells to the model. 

### A vertical producer well
We will define two wells: A first well is named "Producer" and is a vertical well positioned at `(1, 1)`. By default, the [`vertical_well`](@ref) function perforates all layers in the model.

```@example intro_wells
## Set up a vertical well in the first corner, perforated in all layers
Prod = setup_vertical_well(domain, 1, 1, name = :Producer);
```

### A single-perforation injector
We also define an injector by [`setup_well`](@ref). This function allows us to pass a vector of either cell indices or tuples of logical indices that the well trajectory will follow.

```@example intro_wells
## Set up an injector in the upper left corner
Inj = setup_well(domain, [(nx, ny, 1)], name = :Injector);
```
## Choosing a fluid system

```@example intro_wells
## Set up a two-phase immiscible system and define a density secondary variable
phases = (LiquidPhase(), VaporPhase())
rhoLS = 1000.0
rhoGS = 100.0
rhoS = [rhoLS, rhoGS]
sys = ImmiscibleSystem(phases, reference_densities = rhoS)
```

### Creating the model
```@example intro_wells
## Set up a reservoir model that contains the reservoir, wells and a facility that controls the wells
model, parameters = setup_reservoir_model(domain, sys, wells = [Inj, Prod])
model
```

### Replace the density function with our custom version
```@example intro_wells
c = [1e-6/bar, 1e-4/bar]
ρ = ConstantCompressibilityDensities(p_ref = 1*bar, density_ref = rhoS, compressibility = c)
replace_variables!(model, PhaseMassDensities = ρ)
```

```@example intro_wells
## Set up initial state
state0 = setup_reservoir_state(model, Pressure = 150*bar, Saturations = [1.0, 0.0])
```


```@example intro_wells
## Set up time-steps
dt = repeat([30.0]*day, 12*5)
## Inject a full pore-volume (at reference conditions) of gas
# We first define an injection rate
reservoir = reservoir_model(model);
pv = pore_volume(model, parameters)
inj_rate = sum(pv)/sum(dt)
```

```@example intro_wells
# We then set up a total rate target (positive value for injection)
# together with a corresponding injection control that specifies the
# mass fractions of the two components/phases for pure gas injection,
# with surface density given by the known gas density.
rate_target = TotalRateTarget(inj_rate)
I_ctrl = InjectorControl(rate_target, [0.0, 1.0], density = rhoGS)
# The producer operates at a fixed bottom hole pressure
bhp_target = BottomHolePressureTarget(50*bar)
P_ctrl = ProducerControl(bhp_target)
# Set up the controls. One control per well in the Facility.
controls = Dict()
controls[:Injector] = I_ctrl
controls[:Producer] = P_ctrl
# Set up forces for the whole model. For this example, all forces are defaulted
# (amounting to no-flow for the reservoir).
forces = setup_reservoir_forces(model, control = controls)
```

```@example intro_wells
## Finally simulate!
wd, states, reports = simulate_reservoir(state0, model, dt, parameters = parameters, forces = forces);
```