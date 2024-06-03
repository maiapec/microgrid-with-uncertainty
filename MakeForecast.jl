begin
    using JuMP
    using CSV
    using DataFrames
    using PlotlyJS
    using Dates
    using XLSX
    using FileIO
    using Base
    using Random
    using Statistics
    using Gurobi
    using PyCall
    include("ErrorFunctions.jl")
end

# 1. Function to simulate forecast for the weather data (only interested in non-PV variables)

function Make_weather_forecast_full_year(weather::DataFrame, weather_next::DataFrame, std_schedule, errorfunction_dict, horizon=24*7)
    """ Simulates forecast for weather quantities in variables based on the given errorfunction.
        Input: weather data for the year of interest and the next one, to have forecast available for the whole year.
        errorfunction_dict is a dictionary with the error functions for each variable.
        std_schedule is the variability schedule used in the error functions.

        Returns: a list weather_forecast_list where weather_forecast_list[i] is the weather forecast DataFrames for the week starting at hour i.
        Temperatures are converted to [°F].
    """
    begin
        # For debugging purposes
        # display(weather)
        # display(weather_next)

        # Concatenate weather and the first week of the next year 
        weather_all = vcat(weather, weather_next[1:horizon, :])

        # Create a list that stores the weather forecasts dataframes to be used at each hour for the upcoming week
        weather_forecast_list = []

        for i in 1:8760
            weather_week = copy(weather_all[i:(i-1+horizon), :])

            for (var, (error_function, params)) in pairs(errorfunction_dict)
                std_schedule_var = std_schedule[i:(i-1+horizon), var]
                weather_week[!, var] = error_function(weather_week[!, var], std_schedule_var, params...)
            end
            
            # Add the forecasted week to the list
            push!(weather_forecast_list, weather_week)

        end
    end 

    return weather_forecast_list

end

# 2. Once the weather forecast has been simulated, we can calculate the PV output

function Make_PV_output(weather,
                        noct_installed, module_height, wind_height, 
                        module_emissivity, module_absorption, 
                        module_surface_tilt, module_width, module_length, Berg_tilt,
                        horizon=24*7)
    """ Computes the PV output for a given weather DataFrame.
        If used for forecasting PV from weather forecast data, this function should be called after the weather forecast has been simulated.
        Returns the weather DataFrame with the PV output added.
    """
    begin
    # Corrected function to create DateTime from row
    begin
        # Apply function to each row to create a new column with DateTime objects
        weather.datetime = [create_datetime_from_row(row) for row in eachrow(weather)] 
    end
    
    # Use PyCall to format data into Pandas
    begin
        # Import necessary Python libraries
        pd = pyimport("pandas")
    
        # Convert Julia DataFrame to Pandas DataFrame
        columns = names(weather)
        data = [getproperty(weather, col) for col in columns]
        pandas_df = pd.DataFrame(Dict(zip(columns, data)))
        
        # Ensure the datetime column is in the correct format and set as the index
        pandas_df["datetime"] = pd.to_datetime(pandas_df["datetime"])
        pandas_df = pandas_df.set_index("datetime")
    end

    # No more need to simulate forecast for temperature as it has been done in the weather forecast previously.

    # Calculate PV module cell temperature with PVLib [°C]
    begin
        # Proceed with pvlib operations
        pvlib = pyimport("pvlib")
    
        # Prepare data input for pvlib cell temperature calculation
        ambient_temperature = pandas_df["Temperature"] # [°C] independent variable
        wind_speed = pandas_df["Wind Speed"] # [m/s] independent variable
        total_irradiance = pandas_df["GHI"] # [W/m^2] independent variable
    
        # Calculate PV cell temperature using the Fuentes Model
        temp_fuentes = pvlib.temperature.fuentes(
            total_irradiance,
            ambient_temperature,
            wind_speed,
            noct_installed,
            module_height=module_height,
            wind_height=wind_height,
            emissivity=module_emissivity,
            absorption=module_absorption,
            surface_tilt=module_surface_tilt,
            module_width=module_width,
            module_length=module_length
        ) # [°C]

        # Store PV cell temperature into weather file
        temp_fuentes_array_manual = [temp_fuentes.iloc[i] for i in 1:length(temp_fuentes)]
        weather.celltemp = temp_fuentes_array_manual
    end
    
    # Calculate PV output [kW DC Output/kW Capacity]
    begin
        # Calculate PV output [kW DC Output/kW Capacity]  
        PV_output = [];
        for i = 1:(horizon)
            pv_o = calculate_PV(weather[i, :datetime], weather[i, :DNI], weather[i, :DHI], weather[i, :GHI], 180-Berg_tilt, module_surface_tilt, weather[i, :celltemp])
            push!(PV_output, pv_o)
        end
        
        # Store PV output into weather file
        weather.PV = PV_output
    end

    # Add solar time column to the weather file
    weather = add_solar_time(weather)
    
    end
    
    return weather
    
end

function Make_PV_forecast_full_year(weather_forecast_list,
                                    noct_installed, module_height, wind_height, 
                                    module_emissivity, module_absorption, 
                                    module_surface_tilt, module_width, module_length, Berg_tilt,
                                    horizon=24*7)
    """ Computes the PV output for a given weather forecast list (list of forecasts of the upcoming week for each hour of the year)
        and PV parameters.
        Returns the weather forecast list with the PV output added. No additional noise is added to the PV output currently.
    """                         
    weather_forecast_with_PV = []

    println("\nMaking forecast for the year")
    
    for i in 1:8760
        weather_week = copy(weather_forecast_list[i])
        # Calculate PV output
        weather_week = Make_PV_output(weather_week, noct_installed, module_height, wind_height, module_emissivity, module_absorption, module_surface_tilt, module_width, module_length, Berg_tilt, horizon)
        # Add the forecasted week to the list
        push!(weather_forecast_with_PV, weather_week)

    end

    # At the end, convert temperatures from [°C] to [°F]
    for i in 1:8760
        weather_forecast_with_PV[i].Temperature .= (weather_forecast_with_PV[i].Temperature .* (9/5)) .+ 32 # [°F]
    end

    # println(first(weather_forecast_with_PV))
    println("Forecasts for the year have been made")
    return weather_forecast_with_PV
    
end

# 3. Test functions to try out the error functions

function Test_smooth_uncertainty()
    """ Test the uncertainty function by plotting the original and simulated data.
    """
    file_name = "NREL_NSRD_HalfMoonBay25_60Min/137344_37.49_-122.42_1998.csv" # hourly data
    weather = DataFrame(CSV.File(file_name, header=3))

    # Make variability schedule
    std_schedule = make_std_schedule()
    # Extract indexes from 1 to 24*7
    std_schedule = std_schedule[1:24*7, :]

    # Apply error function to temperature data
    weekly_temp = weather.Temperature[1:24*7]
    new_temperature = Add_smooth_uncertainty_to_variable(weekly_temp, std_schedule.Temperature, false, false)
    plot([PlotlyJS.scatter(x=1:24*7, y=weekly_temp, mode="lines", name="Original"),PlotlyJS.scatter(x=1:24*7, y=new_temperature, mode="lines", name="Simulated")])

    # Apply error function to RH
    weekly_RH = weather[!,"Relative Humidity"][1:24*7]
    new_RH = Add_smooth_uncertainty_to_variable(weekly_RH, std_schedule[!,"Relative Humidity"], true, false)
    plot([PlotlyJS.scatter(x=1:24*7, y=weekly_RH, mode="lines", name="Original"),PlotlyJS.scatter(x=1:24*7, y=new_RH, mode="lines", name="Simulated")])

    # Apply error function to wind speed
    weekly_ws = weather[!,"Wind Speed"][1:24*7]
    new_ws = Add_smooth_uncertainty_to_variable(weekly_ws, std_schedule[!,"Wind Speed"], false, true)
    plot([PlotlyJS.scatter(x=1:24*7, y=weekly_ws, mode="lines", name="Original"),PlotlyJS.scatter(x=1:24*7, y=new_ws, mode="lines", name="Simulated")])

end

function Test_random_uncertainty()
    """ Test the uncertainty function by plotting the original and simulated data.
    """
    file_name = "NREL_NSRD_HalfMoonBay25_60Min/137344_37.49_-122.42_1998.csv" # hourly data
    weather = DataFrame(CSV.File(file_name, header=3))

    # Make variability schedule
    std_schedule = make_std_schedule()
    # Extract indexes from 1 to 24*7
    std_schedule = std_schedule[1:24*7, :]

    # Apply error function to temperature data
    weekly_temp = weather.Temperature[1:24*7]
    new_temperature = Add_random_uncertain_points_to_variable(weekly_temp, std_schedule.Temperature, false, false)
    plot([PlotlyJS.scatter(x=1:24*7, y=weekly_temp, mode="lines", name="Original"),PlotlyJS.scatter(x=1:24*7, y=new_temperature, mode="lines", name="Simulated")])

    # Apply error function to RH
    weekly_RH = weather[!,"Relative Humidity"][1:24*7]
    new_RH = Add_random_uncertain_points_to_variable(weekly_RH, std_schedule[!,"Relative Humidity"], true, false)
    plot([PlotlyJS.scatter(x=1:24*7, y=weekly_RH, mode="lines", name="Original"),PlotlyJS.scatter(x=1:24*7, y=new_RH, mode="lines", name="Simulated")])

    # Apply error function to wind speed
    weekly_ws = weather[!,"Wind Speed"][1:24*7]
    new_ws = Add_random_uncertain_points_to_variable(weekly_ws, std_schedule[!,"Wind Speed"], false, true)
    plot([PlotlyJS.scatter(x=1:24*7, y=weekly_ws, mode="lines", name="Original"),PlotlyJS.scatter(x=1:24*7, y=new_ws, mode="lines", name="Simulated")])
end

function Test_noisy_uncertainty()
    """ Test the uncertainty function by plotting the original and simulated data.
    """
    file_name = "NREL_NSRD_HalfMoonBay25_60Min/137344_37.49_-122.42_1998.csv" # hourly data
    weather = DataFrame(CSV.File(file_name, header=3))

    # Make variability schedule
    std_schedule = make_std_schedule()
    # extract indexes from 1 to 24*7
    std_schedule = std_schedule[1:24*7, :]

    # Apply error function to temperature data
    weekly_temp = weather.Temperature[1:24*7]
    new_temperature = Add_noise_to_variable(weekly_temp, std_schedule.Temperature, false, false)
    plot([PlotlyJS.scatter(x=1:24*7, y=weekly_temp, mode="lines", name="Original"),PlotlyJS.scatter(x=1:24*7, y=new_temperature, mode="lines", name="Simulated")])

    # Apply error function to RH
    weekly_RH = weather[!,"Relative Humidity"][1:24*7]
    new_RH = Add_noise_to_variable(weekly_RH, std_schedule[!,"Relative Humidity"], true, false)
    plot([PlotlyJS.scatter(x=1:24*7, y=weekly_RH, mode="lines", name="Original"),PlotlyJS.scatter(x=1:24*7, y=new_RH, mode="lines", name="Simulated")])

    # Apply error function to wind speed
    weekly_ws = weather[!,"Wind Speed"][1:24*7]
    new_ws = Add_noise_to_variable(weekly_ws, std_schedule[!,"Wind Speed"], false, true)
    plot([PlotlyJS.scatter(x=1:24*7, y=weekly_ws, mode="lines", name="Original"),PlotlyJS.scatter(x=1:24*7, y=new_ws, mode="lines", name="Simulated")])

end

function Test_make_forecast_full()
    """ Test the Make_weather_forecast_full_year function by printing the first week of the forecast.
    """
    # Define the errorfunction_dict
    
    # errorfunction_dict = Dict(
    #     "Temperature" => (Add_smooth_uncertainty_to_variable, (0.9, 0.5, 5, false)),
    #     "Wind Speed" => (Add_random_uncertain_points_to_variable, (5, 0.2, false)),
    #     "Relative Humidity" => (Add_noise_to_variable, (0.5, 0.3, true))
    # )
    
    errorfunction_dict = Dict(
        "Temperature" => (Add_smooth_uncertainty_to_variable, (0.9, 0.5, 5, false, false)),
        #"Wind Speed" => (Add_smooth_uncertainty_to_variable, (0.9, 0.1, 5, false, true)),
        "Relative Humidity" => (Add_noise_to_variable, (0.5, 0.3, true, false))
    )
    weather = DataFrame(CSV.File("NREL_NSRD_HalfMoonBay25_60Min/137344_37.49_-122.42_1998.csv", header=3))
    weather_next = DataFrame(CSV.File("NREL_NSRD_HalfMoonBay25_60Min/137344_37.49_-122.42_1999.csv", header=3))

    # Make variability schedule
    std_schedule = make_std_schedule()

    # Make forecast for the whole year
    weather_forecast_list = Make_weather_forecast_full_year(weather, weather_next, std_schedule, errorfunction_dict, 24*7)

    # Parameters for PV output calculation
    begin
        # Time Horizon Parameters
        begin
            TimeStart = 1;
            TimeEnd = 8760;
            f_run = 1 # run frequency [hour/time]
            Opt_Horizon = 168 # [hr]
            stepsize = 60*(1/f_run) # [min]
            δt = stepsize/60 # [hr] Declare stepzize for the optimization program
            NumRun = (TimeEnd-TimeStart+1) - Opt_Horizon + 1; # This is the max number of runs allowed
            # NumRun = 1
        end     
        # Berg Envelope Parameters
        begin
            L_wall = 20 # [ft] Length of the Wall
            H_wall = 8 # [ft] Height of the Wall
            R_wall = 15 # [ft^2·°F·hr/BTU] Thermal Resistance of the Wall
            R_floor = 17 # [ft^2·°F·hr/BTU] Thermal Resistance of the Floor
            H_Ceiling = 80/12 # [ft] Interior Height of the Ceiling (80 inch) 
            Area = L_wall * L_wall # [ft^2] Floor/Ceiling Area 
            Berg_tilt = 15 # [DEGREES] The Berg structure is tilted 15 degrees east of true south
            UA = L_wall*H_wall*4*(1/R_wall) + Area*2*(1/R_floor) # [BTU/(hr⋅°F)]
            TC = 581.7 # [BTU/°F] Total Thermal Capacitance of the Structure - calibrated using passive model 4.6
            Volume = 2330 # [ft^3]
            #=
            UA = 89, which is consistent with given documents. Note that this is excluding infiltration, since infiltration is considered 
            as a heat transfer q_Infiltration instead of extra UA.
            =#
        end
        # Radiation Parameters
        begin
            e_Berg = 0.75 # currently not used
            SHGC = 0.0215 # calibrated using passive model 4.6
            RCC = 0.0408 # calibrated using passive model 4.6
        end
        # Environmental Parameters
        begin
            wind_height = 9.144 # [m] The height above ground at which wind speed is measured. The PVWatts default is 9.144 m.
            albedo = 0.18 # [1] albedo of grass
            n_air = 1 # [1] Index of reflaction of air
            n_AR = 1.3 # [1] Index of reflaction of AR coating
            n_glass = 1.526 # [1] Index of reflaction of glass as in PVWatt model for standard modules
            standard_meridian_longitude = -120; # [DEGREES] Longitude of standard meridian
            longitude = -122.4286 # [DEGREES] Longitude of site
            latitude = 37.4636 # [DEGREES] Latitude of site
        end
        # Infiltration Parameters
        begin
            # Stack coefficient for building(house height = 1 story)
            Cs = 0.000145 # [(L/s)^2/(cm^4*K)]
            # Wind coefficient for building(house height =1 story, shelter class=1(no obstructions or local shielding))
            Cw = 0.000319 # [(L/s)^2/(cm^4*(m/s)^2)] (between 0.000319 and 0.000246 (house height =1 story, shelter class=2(typical shelter for an insolated rural house)))
            
            # Effective leakage area measured during Blower Door Test 
            ELA =  38.3 # [in^2] (between 38.3(Dan's interpolation) and 47.1(Jessie's interpolation))
            Al = ELA * 6.4516 # [cm^2]
    
            T_d = 273.15 - 1.67 # [°K] 99% Design Temperature for Half Moon Bay
            T_indoor_constant = 22 # [°C] Constant Indoor Temperature (for the simplicity of a linear model)
        end
        # Daylight Saving
        begin
            dls_start = DateTime(2024, 3, 10, 2)
            dls_end = DateTime(2024, 11, 3, 0)
            timezone = "America/Los_Angeles"  # Specify local timezone
        end
        # Interior Standard
        begin
            # Temperature standard 
            SetPointT_Low = 68 # [°F] Lower Bound of Set Point Temperature
            SetPointT_High = 75 # [°F] Upper Bound of Set Point Temperature
    
            # Relative Humidity standard (currently not used due to lack of humidity data)
            # SetPointW_Low = 0.4 # LEED Regulation
            # SetPointW_High = 0.6 # LEED Regulation
    
            # Ventilation
            Ventilation = 15 # [CFM/PPL] Building Standard
            Ra = 0.06 # [CFM/ft^2] Area Outdoor Air Rate; Source:CEE226E Week 5 - Energy Modeling Questions - Slide 9
            Rp = 5 # [CFM/PPL] People Outdoor Air Rate; Source:CEE226E Week 5 - Energy Modeling Questions - Slide 9   
        end
        # Lighting, Plugs, and Occupancy
        begin
            PeakLighting = 0.1 # [KW] Dan's suggestion
            PeakPlugLoad = 0.1 # [KW] Computer = 40W, Phone = 10W, 2*(40 + 10) = 100 [W] = 0.1 [kW]
            MaxOccupancy = 4 # [PPL] 4-MAN Office Room Plan
            PersonLatentHeat = 200; # [BTU/hr/PPL]  CEE226E Slide
            PersonSensibleHeat = 300; # [BTU/hr/PPL]  CEE226E Slide
            TotalPersonHeat = PersonSensibleHeat + PersonLatentHeat # [BTU/hr/PPL]
            # TotalPersonHeat = 300 # [BTU/hr/PPL] Dan's suggestion
        end
        # Equipments Parameters
        begin
            # Device Capacity Parameters
            PVSize = 3; # [kW] PV DC Power Capacity
            BatterySize = 5; # [kWh] Battery Energy Capacity
            InverterSize = 15 # [kW] Max Continuous AC Output Power
            PCM_H_Size = 10 * 3412.14; # [BTU] PCM Heating Storage Energy Capacity
            PCM_C_Size = 1 * 3412.14; # [BTU] PCM Cooling Storage Energy Capacity
    
            # Solar PV Parameters
            noct_installed = 45 # [°C] The “installed” nominal operating cell temperature. PVWatts assumes this value to be 45 C for rack-mounted arrays and 49 C for roof mount systems with restricted air flow around the module.
            module_height = 5 # [m] The height above ground of the center of the module. The PVWatts default is 5.0.
            module_width = 0.31579 # [m] Module width. The default value of 0.31579 meters in combination with the default module_length gives a hydraulic diameter of 0.5.
            module_length = 1.2 # [m] Module length. The default value of 1.2 meters in combination with the default module_width gives a hydraulic diameter of 0.5.
            module_emissivity = 0.84 # [1] The effectiveness of the module at radiating thermal energy.
            module_absorption = 0.83 # [1] The fraction of incident irradiance that is converted to thermal energy in the module. 
            module_surface_tilt = 27 # [DEGREES] Module tilt from horizontal. If not provided, the default value of 30 degrees is used.
            Γ_t = -0.47/100 # [1/°C] Temperature coefficient for standard module
            T_ref = 25 # [°C] Reference cell temperature
            η_PV = 0.86 # [1] PV DC efficiency after system loss
            P_dc0 = 1 # [kW/kW] Rated capacity at standard conditions
            η_PVIV = 0.94 # [1] PV(DC) to Home(AC) inverter efficiency
            
            # Battery parameters
            BatteryLoss = 0.00001 # [/hr] Battery Leakage Rate
            MaxDischarge = 0.8 # [1] Max Depth of Discharge in Battery (80%)
            η = 0.98 # [1] Battery Inverter Efficiency
    
            # LBNL System Parameters
            # Heat Pump heating and cooling capacity
            Cap_HP_C = 3.85 # [kW] Max input electrical power for heat pump (cooling)
            Cap_HP_H = 4 # [kW] Max input electrical power for heat pump (heating)
            
            # Heat Pump heating and cooling coefficients of performance
            HP_a = 6.08
            HP_b = -0.0941
            HP_c = 0.000464
    
            C_HP_OP = 0.02 * Cap_HP_H # [$/(kW*YR)] Operational Cost of Heat Pump
    
            # Thermal Energy Storage - Phase Changing Materials
            Cap_PCM_H_Charging = Cap_HP_H * 4 * 3412.14 # [BTU/hr] Max Charging Power Capacity of PCM Heating Storage
            Cap_PCM_H_Discharging = Cap_PCM_H_Charging # [BTU/hr] Max Discharging Power Capacity of PCM Heating Storage
            Cap_PCM_C_Charging = Cap_HP_C * 4 * 3412.14 # [BTU/hr] Max Charging Power Capacity of PCM Cooling Storage
            Cap_PCM_C_Discharging = Cap_PCM_C_Charging # [BTU/hr] Max Discharging Power Capacity of PCM Cooling Storage
    
            η_PCM_H = 0.99 # [1] PCM Heating Storage Efficiency
            η_PCM_C = 0.99 # [1] PCM Cooling Storage Efficiency
    
            C_PCM_H_OP = 0.02 * PCM_H_Size # [$/(kWh*YR)] Operational Cost of PCM Heating Storage
            C_PCM_C_OP = 0.02 * PCM_C_Size # [$/(kWh*YR)] Operational Cost of PCM Cooling Storage
    
            # Standard Operating Power of Heat Pump and PCM Thermal Storages 
            P_H2HP = Cap_HP_H # [kW] default constant electrical power consumption for heat pump (heating)
            P_H2C = Cap_HP_C # [kW] default constant electrical power consumption for heat pump (cooling)
            P_HP2PCM_H = Cap_PCM_H_Charging # [BTU/hr] default constant heat charging rate of PCM Heating Storage
            P_C2PCM_C = Cap_PCM_C_Charging # [BTU/hr] default constant heat charging rate of PCM Cooling Storage
            P_PCM_H2H = Cap_PCM_H_Discharging # [BTU/hr] default constant heat discharging rate of PCM Heating Storage
            P_PCM_C2H = Cap_PCM_C_Discharging # [BTU/hr] default constant heat discharging rate of PCM Cooling Storage
            
            # Big M Method
            M = 10000
        end
    end

    weather_forecast_with_PV = Make_PV_forecast_full_year(weather_forecast_list,
                                    noct_installed, module_height, wind_height, 
                                    module_emissivity, module_absorption, 
                                    module_surface_tilt, module_width, module_length, Berg_tilt, 
                                    24*7)

    println(first(weather_forecast_with_PV))

end

# 4. Plotting functions to visualize the error functions (not called in the main program)

function Plot_example_Temperature()
    """ Plot examples of temperature forecasts for different error function for the first week of the year.
    """

    file_name = "NREL_NSRD_HalfMoonBay25_60Min/137344_37.49_-122.42_1998.csv" # hourly data
    weather = DataFrame(CSV.File(file_name, header=3))

    # Make variability schedule
    std_schedule = make_std_schedule()
    # Extract indexes from 1 to 24*7
    std_schedule = std_schedule[1:24*7, :]

    # Apply error function to temperature data
    weekly_temp = weather.Temperature[1:24*7]
    new_temperature_smooth = Add_smooth_uncertainty_to_variable(weekly_temp, std_schedule.Temperature, 0.95, 2.0, 7, false, false)
    new_temperature_random = Add_random_uncertain_points_to_variable(weekly_temp, std_schedule.Temperature, 6, 0.5, false, false)
    new_temperature_noisy = Add_noise_to_variable(weekly_temp, std_schedule.Temperature, 0.5, 0.5, false, false)

    plot([PlotlyJS.scatter(x=1:24*7, y=new_temperature_noisy, mode="lines", name="Noisy", line=attr(width=1.7, color="rgb(255, 170, 180)")),
        PlotlyJS.scatter(x=1:24*7, y=new_temperature_random, mode="lines", name="Random", line=attr(width=1.7, color="rgb(120,220,160)")),
        PlotlyJS.scatter(x=1:24*7, y=new_temperature_smooth, mode="lines", name="Smooth", line=attr(width=1.7, color="rgb(0,180,240)")),
        PlotlyJS.scatter(x=1:24*7, y=weekly_temp, mode="lines", name="Actual", line=attr(width=1.5, color="black")),
        
        ],
        Layout(
        plot_bgcolor="white",   # Set the plot background color to white
        paper_bgcolor="white",  # Set the paper background color to white
        xaxis=attr(
            title="Time (hours)",  # Set the x-axis title
            showline=true,         # Show axis line
            linecolor="black",     # Axis line color
            showgrid=false,        # Hide grid lines
            tickcolor="black"      # Tick color
        ),
        yaxis=attr(
            title="Temperature (°C)",  # Set the y-axis title
            showline=true,             # Show axis line
            linecolor="black",         # Axis line color
            showgrid=false,            # Hide grid lines
            tickcolor="black"          # Tick color
        ),
        font=attr(
            family="Arial",   # Set the font family for the entire plot
            size=12,          # Set the font size for the entire plot
            color="black"     # Set the font color for the entire plot
        ),
        width=500,
        height=400  
    )
    )

end

function Plot_example_WindSpeed()
    """ Plot examples of temperature forecasts for different error function for the first week of the year.
    """

    file_name = "NREL_NSRD_HalfMoonBay25_60Min/137344_37.49_-122.42_1998.csv" # hourly data
    weather = DataFrame(CSV.File(file_name, header=3))

    # Make variability schedule
    std_schedule = make_std_schedule()
    # Extract indexes from 1 to 24*7
    std_schedule = std_schedule[1:24*7, :]

    # Apply error function to temperature data
    weekly_ws = weather[!,"Wind Speed"][1:24*7]
    new_ws_smooth = Add_smooth_uncertainty_to_variable(weekly_ws, std_schedule[!,"Wind Speed"], 0.95, 2.0, 5, false, true)
    new_ws_random = Add_random_uncertain_points_to_variable(weekly_ws, std_schedule[!,"Wind Speed"], 6, 0.5, false, true)
    new_ws_noisy = Add_noise_to_variable(weekly_ws, std_schedule[!,"Wind Speed"], 0.5, 0.5, false, true)

    plot([PlotlyJS.scatter(x=1:24*7, y=new_ws_noisy, mode="lines", name="Noisy", line=attr(width=1.7, color="rgb(255, 170, 180)")),
        PlotlyJS.scatter(x=1:24*7, y=new_ws_random, mode="lines", name="Random", line=attr(width=1.7, color="rgb(120,220,160)")),
        PlotlyJS.scatter(x=1:24*7, y=new_ws_smooth, mode="lines", name="Smooth", line=attr(width=1.7, color="rgb(0,180,240)")),
        PlotlyJS.scatter(x=1:24*7, y=weekly_ws, mode="lines", name="Actual", line=attr(width=1.5, color="black")),
        
        ],
        Layout(
        plot_bgcolor="white",   # Set the plot background color to white
        paper_bgcolor="white",  # Set the paper background color to white
        xaxis=attr(
            title="Time (hours)",  # Set the x-axis title
            showline=true,         # Show axis line
            linecolor="black",     # Axis line color
            showgrid=false,        # Hide grid lines
            tickcolor="black"      # Tick color
        ),
        yaxis=attr(
            title="Wind speed (m/s)",  # Set the y-axis title
            showline=true,             # Show axis line
            linecolor="black",         # Axis line color
            showgrid=false,            # Hide grid lines
            tickcolor="black"          # Tick color
        ),
        font=attr(
            family="Arial",   # Set the font family for the entire plot
            size=13,          # Set the font size for the entire plot
            color="black"     # Set the font color for the entire plot
        ),
        width=500,
        height=400
    )
    )

end

# 5. Functions to quantify the error in the forecast

function Calculate_RMSE(weather_forecast_with_PV, weather, weather_next, variable)
    """ Calculate the RMSE for a given variable in the forecast.
    Inputs: 
    - weather_forecast_with_PV: the list of weather forecasts for each week (or time horizon) for the whole year
    - weather: the actual weather data for the year
    """
    begin
        # Copy data
        wf = copy(weather_forecast_with_PV)
        w = copy(weather)
        wn = copy(weather_next)
        # Extract the variable from the weather DataFrame
        horizon = size(wf[1])[1] # typically 1 week
        forecasted_variable_list = [wf[i][!, variable] for i in 1:8760]
        
        # Repeat actual_variable at each hour for the upcoming horizon time
        weather_all = vcat(w, wn[1:horizon, :])
        actual_variable = weather_all[!, variable]
        actual_variable_list = [actual_variable[i:(i-1+horizon),:] for i in 1:8760]

        # Calculate the RMSE as the average squared error over a sample of length equal to the time horizon
        RMSE = 0
        for i in 1:8760
            RMSE += mean((forecasted_variable_list[i] .- actual_variable_list[i]).^2)
        end
        RMSE = sqrt(RMSE/8760)

    end
    return RMSE
end

function Quantify_error_in_forecast(weather_forecast_with_PV, weather, weather_next)
    """ Quantify the error in the forecast by calculating the RMSE for the temperature and PV output.
    Inputs: 
    - weather_forecast_with_PV: the list of weather forecasts for each week (or time horizon) for the whole year.
    Note that the temperature is in [°F] in the weather_forecast_with_PV list.
    - weather: the actual weather data for the year
    """
    begin
        # Calculate RMSE for temperature
        wf = copy(weather_forecast_with_PV)
        w = copy(weather)
        wn = copy(weather_next)
        # Change temperatures from weather_forecast_with_PV back to [°C] for the calculation
        for i in 1:8760
            wf[i].Temperature .= (wf[i].Temperature .- 32) .* (5/9) # [°C]
        end
        RMSE_Temp = Calculate_RMSE(wf, w, wn, "Temperature")

        # If very low, change to zero (this is due to computational errors when converting back and forth between °C and °F)
        if RMSE_Temp < 1e-6
            RMSE_Temp = 0.0
        end
        println("RMSE for temperature: ", RMSE_Temp)

        # Calculate RMSE for wind speed
        RMSE_WS = Calculate_RMSE(wf, w, wn, "Wind Speed")
        if RMSE_WS < 1e-6
            RMSE_WS = 0.0
        end
        println("RMSE for wind speed: ", RMSE_WS)

        # Calculate RMSE for relative humidity
        RMSE_RH = Calculate_RMSE(wf, w, wn, "Relative Humidity")
        if RMSE_RH < 1e-6
            RMSE_RH = 0.0
        end
        println("RMSE for relative humidity: ", RMSE_RH)

    end
    return RMSE_Temp, RMSE_WS, RMSE_RH
end
