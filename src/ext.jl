export plot_well!, plot_well_results, plot_reservoir_simulation_result

function plot_well!

end

function plot_well_results

end

"""
    plot_reservoir_simulation_result(model::MultiModel, res::ReservoirSimResult; wells = true, reservoir = true)

Plot a reservoir simulation result. If `wells=true` well curves will be shown
interactively. If `reservoir=true` the reservoir quantities will be visualized
in 3D. These options can be combined.
"""
function plot_reservoir_simulation_result(model::MultiModel, res::ReservoirSimResult; wells = true, reservoir = true)
    Jutul.check_plotting_availability()
    if reservoir
        rmodel = reservoir_model(model)
        fig = plot_interactive(rmodel, res.states)
        g = physical_representation(rmodel.data_domain)
        ax = fig.current_axis[]
        for (k, m) in pairs(model.models)
            w = physical_representation(m.data_domain)
            if w isa WellDomain
                plot_well!(ax, g, w)
            end
        end
    else
        fig = nothing
    end
    if wells
        plot_well_results(res.wells, res.time, new_window = true)
    end
    if reservoir
        display(fig)
    end
    return fig
end

export plot_reservoir
function plot_reservoir(model, arg...; well_fontsize = 18, well_linewidth = 5, kwarg...)
    rmodel = reservoir_model(model)
    fig = plot_interactive(rmodel.data_domain, arg...; kwarg...)
    g = physical_representation(rmodel.data_domain)
    ax = fig.current_axis[]
    for (k, m) in pairs(model.models)
        w = physical_representation(m.data_domain)
        if w isa WellDomain
            plot_well!(ax, g, w, fontsize = well_fontsize, linewidth = well_linewidth)
        end
    end
    return fig
end

export simulate_reservoir_parray
function simulate_reservoir_parray(case, mode = :mpi; kwarg...)
    sim, cfg = setup_reservoir_simulator(case; mode = mode, kwarg...)
    return simulate!(sim, case.dt, forces = case.forces, config = cfg)
end

export setup_reservoir_simulator_parray
function setup_reservoir_simulator_parray

end
