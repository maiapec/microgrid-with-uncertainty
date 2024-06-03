using CSV
using DataFrames
using Plots

include("Util.jl")


function plot_accumulated_loss_of_load(loss_of_load_dict)
    # Get unique years, MPCs, and scenarios
    years = unique(key[3] for key in keys(loss_of_load_dict))
    mpcs = unique(key[1] for key in keys(loss_of_load_dict))
    scenarios = unique(key[2] for key in keys(loss_of_load_dict))

    # Determine the number of subplots needed
    num_years = length(years)
    num_cols = 2
    num_rows = cld(num_years, num_cols)  # ceiling division to determine rows

    # Initialize the plot with a dynamic layout
    p = Plots.plot(title="Accumulated Loss of Load Over the Year", layout=(num_rows, num_cols))

    # Loop over each year to create subplots
    for (i, year) in enumerate(years)
        subplot = Plots.plot(title="Year: $year", xlabel="Time (Hours)", ylabel="Accumulated Loss of Load")
        
        for mpc in mpcs
            for scenario in scenarios
                for (key, (loss_of_load)) in loss_of_load_dict
                    if key[3] == year && key[1] == mpc && key[2] == scenario
                        # Calculate the accumulated loss of load
                        accumulated_loss_of_load = cumsum(loss_of_load)
                        # Debug output
                        println("Plotting for Year: $year, MPC: $mpc, Scenario: $scenario")
                        println("Accumulated Loss of Load: ", accumulated_loss_of_load[end])
                        # Plot the line
                        Plots.plot!(subplot, accumulated_loss_of_load, label="MPC: $mpc, Scenario: $scenario")
                    end
                end
            end
        end
        
        # Add the subplot to the main plot
        display(Plots.plot(subplot))
    end

    return p
end

### Plotting heat map ###

using Plots

function plot_performance_heatmap(loss_of_load_dict)
    # Get unique years, MPCs, and scenarios
    years = unique(key[3] for key in keys(loss_of_load_dict))
    mpcs = unique(key[1] for key in keys(loss_of_load_dict))
    scenarios = unique(key[2] for key in keys(loss_of_load_dict))

    # Initialize a dictionary to hold total loss of load per scenario and MPC
    total_loss_of_load = Dict((mpc, scenario, year) => 0.0 for mpc in mpcs, scenario in scenarios, year in years)

    # Calculate the total loss of load for each scenario and MPC
    for (key, loss_of_load) in loss_of_load_dict
        mpc, scenario, year = key[1], key[2], key[3]
        total_loss_of_load[(mpc, scenario, year)] = loss_of_load[end]
    end

    # Assuming "baseline" is one of the scenarios, calculate the average relative loss of load
    baseline = "Baseline"
    relative_loss_of_load = Dict((mpc, scenario) => [] for mpc in mpcs, scenario in scenarios )

    for mpc in mpcs
        for scenario in scenarios
            for year in years
                baseline_loss = total_loss_of_load[(baseline, scenario, year)]
                scenario_loss = total_loss_of_load[(mpc, scenario, year)]
                relative_loss = (scenario_loss - baseline_loss) / baseline_loss
                push!(relative_loss_of_load[(mpc, scenario)], relative_loss)
            end
        end
    end

    # Calculate the average relative loss of load
    avg_relative_loss_of_load = Dict((mpc, scenario) => mean(relative_loss_of_load[(mpc, scenario)]) for mpc in mpcs, scenario in scenarios)

    # Convert the dictionary to a matrix for heatmap plotting
    relative_loss_matrix = [avg_relative_loss_of_load[(mpc, scenario)] for mpc in mpcs, scenario in scenarios]
    # Create heatmap
    heatmap(mpcs, scenarios[2:end], relative_loss_matrix, xlabel="MPC", ylabel="Scenario", zlabel="Average Relative Loss of Load", title="Performance Heatmap")
end


# Example usage to retrieve data
log_file = "mpc_run_log.csv"
loss_of_load_dict = retrieve_loss_of_load(log_file)

display(loss_of_load_dict)

# Generate and display the plot
p = plot_accumulated_loss_of_load(loss_of_load_dict)
display(p)

# p2 = plot_performance_heatmap(loss_of_load_dict)
# display(p2)