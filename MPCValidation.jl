# Cross-validation function

function cross_validate_mpc(mpc_functions, weather_scenarios, simulated_years, weather_years, TimeStart, NumRun, Opt_Horizon, Schedule, log_file)

    # Prepare weather variability schedule
    std_schedule = make_std_schedule()
    
    for year in simulated_years[1:end-1]
        
        current_year_data = weather_years[year-first(simulated_years)+1] # get current and next weather years
        next_year_data = weather_years[year-first(simulated_years) + 2]
        local weather = current_year_data
        local weather_next = next_year_data

        println("Loaded weather")
        
        # Run baseline scenario for each year
        println("\n-------------------------------------------")
        println("Running Baseline")
        println("-------------------------------------------")
        runtime = @elapsed begin
            weather_baseline = copy(weather)
            current_states = [0.5*BatterySize, 0.5*PCM_H_Size, 0.5*PCM_C_Size, 71]
            baseline_loss_of_load = Optimize_Baseline(weather_baseline, Schedule, current_states) # loss of load
        end
        log_run(log_file, "Baseline", "None", year, runtime, Opt_Horizon, baseline_loss_of_load)
        println("\nDone with baseline")

        # Loop through all weather forecast scenarios
        for (scenario_name, forecast_error_dict) in weather_scenarios

            # Prepare weather forecast (same to use for all mpc functions)
            begin
                # Make forecast for the whole year
                Weather_forecast_list = Make_weather_forecast_full_year(weather, weather_next, std_schedule, forecast_error_dict, Opt_Horizon)
                # Add PV
                Weather_forecast_with_PV = Make_PV_forecast_full_year(Weather_forecast_list, 
                                                noct_installed, module_height, wind_height, 
                                                module_emissivity, module_absorption, 
                                                module_surface_tilt, module_width, module_length, Berg_tilt, Opt_Horizon)
                # # Quantify error in the forecast
                # RMSE_Temp, RMSE_WS, RMSE_RH = Quantify_error_in_forecast(Weather_forecast_with_PV, weather, weather_next)
                # forecast_error = [RMSE_Temp, RMSE_WS, RMSE_RH]
            end
            
            # Loop through all MPC functions
            for (mpc_name, (mpc_func, mpc_param)) in mpc_functions

                # Unpack parameters
                if length(mpc_param) > 0
                    (hours_penalized,) = mpc_param
                    previous_actions = zeros(13, hours_penalized+1)
                end
                # Re initialize states
                current_states = [0.5*BatterySize, 0.5*PCM_H_Size, 0.5*PCM_C_Size, 71]
                
                # Take runtime of MPC function
                runtime = @elapsed begin
                    println("\n-------------------------------------------")
                    println("Running $mpc_name on $scenario_name")
                    println("-------------------------------------------")
                    Costs, AccumulatedCosts, AllStates, ActionMatrix, TotalCost = mpc_func(Weather_forecast_with_PV, Schedule, Opt_Horizon, TimeStart, NumRun, current_states, mpc_param...)
                end
                # Log only cost profile for now
                log_run(log_file, mpc_name, scenario_name, year, runtime, Opt_Horizon, Costs) 

                trace1 = scatter(
                x = 1:NumRun,  
                y=AccumulatedCosts,
                name="Optimized Storage System"
                )

                p_ll = plot([trace1], Layout(title="Accumulated Loss of Load Over Time", xaxis_title="Hours", yaxis_title="Loss of Load (kWh)"))
                display(p_ll)
                save_plot(p_ll, folder_path, "Accumulated Loss of Load $mpc_name on $scenario_name for year $year)_$today_date", "png")

                println("\nThe total cost is $TotalCost kWh.")
                println("The runtime is $runtime seconds.\n")
            end
        end
    end
end