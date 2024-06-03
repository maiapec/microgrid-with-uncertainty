# Defines the scenarios of error functions to be used in the weather forecast simulations

function get_all_scenarios()
    # Define the error functions to be used in the scenarios
    begin
        # 0. Scenario where MPC is performed on ground truth

        errorfunction_dict_ground_truth = Dict()

        # 1. 3 different scenarios with medium-range smooth error on one variable only (to see which one has more impact - we probably expect temperature):
        # Note: Humidity is not used in this model because data on dew point is available.

        errorfunction_dict_only_temp = Dict(
                "Temperature" => (Add_smooth_uncertainty_to_variable, (0.95, 2.0, 7, false, false))
        )
        errorfunction_dict_only_wind= Dict(
                "Wind Speed" => (Add_smooth_uncertainty_to_variable, (0.95, 2.0, 7, false, true))
        )
        errorfunction_dict_only_dewpoint= Dict(
                "Dew Point" => (Add_smooth_uncertainty_to_variable, (0.95, 2.0, 7, false, false))
        )

        # 2. Scenario with low-range smooth error on all three variables

        errorfunction_dict_all_smooth_low = Dict(
                "Temperature" => (Add_smooth_uncertainty_to_variable, (0.95, 1.0, 5, false, false)), 
                "Wind Speed" => (Add_smooth_uncertainty_to_variable, (0.95, 1.0, 5, false, true)),
                "Dew Point" => (Add_smooth_uncertainty_to_variable, (0.95, 1.0, 5, false, false))
        )

        # 3. Scenario with medium-range smooth error on all three variables

        errorfunction_dict_all_smooth_medium = Dict(
                "Temperature" => (Add_smooth_uncertainty_to_variable, (0.95, 2.0, 7, false, false)), 
                "Wind Speed" => (Add_smooth_uncertainty_to_variable, (0.95, 2.0, 7, false, true)),
                "Dew Point" => (Add_smooth_uncertainty_to_variable, (0.95, 2.0, 7, false, false))
        )

        # 4. Scenario with high-range smooth error on all three variables

        errorfunction_dict_all_smooth_high = Dict(
                "Temperature" => (Add_smooth_uncertainty_to_variable, (0.95, 3.0, 10, false, false)), 
                "Wind Speed" => (Add_smooth_uncertainty_to_variable, (0.95, 3.0, 10, false, true)),
                "Dew Point" => (Add_smooth_uncertainty_to_variable, (0.95, 3.0, 10, false, false))
        )

        # 5. Scenario with medium-range randomly sparsed noise on all three variables

        errorfunction_dict_all_random_medium = Dict(
                "Temperature" => (Add_random_uncertain_points_to_variable, (6, 0.5, false, false)), 
                "Wind Speed" => (Add_random_uncertain_points_to_variable, (6, 0.5, false, true)),
                "Dew Point" => (Add_random_uncertain_points_to_variable, (6, 0.5, false, false))
        )
        # (low would be 5, 0.3)(high would be 7, 0.8)

        # 6. Scenario with medium-range full noise on all three variables

        errorfunction_dict_all_noisy_medium = Dict(
                "Temperature" => (Add_noise_to_variable, (0.5, 0.5, false, false)), 
                "Wind Speed" => (Add_noise_to_variable, (0.5, 0.5, false, true)),
                "Dew Point" => (Add_noise_to_variable, (0.5, 0.5, false, false))
        )
        # (low would be 0.3, 0.3)(high would be 0.7, 0.7)
    end
    # Make a dictionary of all scenarios
    begin
        weather_scenarios = Dict(
            "Ground Truth" => errorfunction_dict_ground_truth,
            "Only Smooth Temp" => errorfunction_dict_only_temp,
            "Only Smooth Wind" => errorfunction_dict_only_wind,
            "Only Smooth Dew Point" => errorfunction_dict_only_dewpoint,
            "All Smooth Low" => errorfunction_dict_all_smooth_low,
            "All Smooth Medium" => errorfunction_dict_all_smooth_medium,
            "All Smooth High" => errorfunction_dict_all_smooth_high,
            "All Random Medium" => errorfunction_dict_all_random_medium,
            "All Noisy Medium" => errorfunction_dict_all_noisy_medium
        )
    end
    return weather_scenarios
end

# Extract only a few scenarios for testing
function get_specific_scenarios(namelist)
    weather_scenarios = get_all_scenarios()
    specific_scenarios = Dict()
    for name in namelist
        # check if the name is in the dictionary
        if !haskey(weather_scenarios, name)
            println("Scenario $name was not defined. Please check in the MakeScenarios.jl file. Skipping.")
            continue
        end
        specific_scenarios[name] = weather_scenarios[name]
    end
    return specific_scenarios
end
