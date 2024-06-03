######## Utility functions for the project ########

function create_datetime_from_row(row)
    """ Create a datetime variable for each time instance
    """
    # Corrected access to DataFrame row fields
    return DateTime(row.Year, row.Month, row.Day, row.Hour, row.Minute)
end

function Moving_average(data, step)
    """ Moving average function. 
        NB: this function does not return an array of the same size as the input
    """
    return [mean(data[i-step+1:i]) for i in step:length(data)]
end

function Smoothing(data, step)
    """ Performs smoothing of data using a moving average method of size step.
        Returns a new array of the same length as data.
    """
        smoothed_data = [data[1]]
        for i in step:length(data)
            push!(smoothed_data, mean(data[i-step+1:i]))
        end
        if length(smoothed_data) < length(data)
            append!(smoothed_data, data[length(smoothed_data)+1:end])
        end
        return smoothed_data
end

# Function to save a plot as a PNG file in the specified folder
function save_plot(plot, path, filename, format="png")
    # Create the full file path with the specified filename and format
    full_path = joinpath(path, string(filename, ".", format))
    
    # Save the plot as an image in the desired format
    savefig(plot, full_path)
end

# Function to store and retrieve results in/from a CSV file
function log_run(file, mpc_name, scenario_name, year, run_time, hours_forecast, loss_of_load)
    # Create a DataFrame for the current run
    df = DataFrame(MPC=[mpc_name], Scenario=[scenario_name], Year=[year], RunTime=[run_time], ForecastHorizon=[hours_forecast], LossOfLoad=[loss_of_load])
    
    # Write to CSV file, appending to it
    open(file, "a") do file
        CSV.write(file, df, append=true)
    end
end

function clean_and_parse(problem_string)
    # Remove any extraneous characters
    clean_str = replace(problem_string, "[" => "", "]" => "")
    # Split and parse the string into an array of Float64
    parse.(Float64, split(clean_str, ","))
end

function retrieve_loss_of_load(file)
    df = DataFrame(CSV.File(file))

    # Convert LossOfLoad column back to array
    loss_of_load_dict = Dict()

    for row in eachrow(df)
        mpc_name = row.MPC
        scenario_name = row.Scenario
        year = row.Year
        forecast_horizon = row.ForecastHorizon
        loss_of_load = clean_and_parse(row.LossOfLoad)
        loss_of_load_dict[(mpc_name, scenario_name, year, forecast_horizon)] = (loss_of_load)
    end

    return loss_of_load_dict
end

######## Data preparation function ########

function DataPreparation(Weather, Schedules, J_States)
    """ Prepares data once PV has been calculated.
    Inputs: Weather that already contains PV, 
    Schedules is either SimpleSchedule or ComplexSchedule, 
    J_States is the initial states of the system.
    Returns all the quantities needed for optimization.
    """
    begin
        # Set timesteps 
        NumTime = size(Weather)[1] # 168 hrs for now
        
        Weather.CFM = zeros(NumTime)
        Weather.RadHeatGain = zeros(NumTime)
        Weather.RadCool = zeros(NumTime)

        for i = 1:NumTime
            Weather.CFM[i] = Al*sqrt(Cs*(T_indoor_constant + 273.15 - T_d) + Cw*Weather[i, "Wind Speed"]^2) * 2.11888 # [ft^3/min] 
            Weather.RadHeatGain[i] = calculate_solarheatgain(Weather[i,:datetime], Weather[i,:DNI], Weather[i,:DHI], Weather[i,:GHI]) # [BTU/hr]
            day_status = IsDay(Weather.SolarTime[i], Weather.datetime[i]) 
            Weather.RadCool[i] = Q_radiativeCooling(T_indoor_constant + 273.15, Weather.Temperature[i], Weather[i,:"Dew Point"], Weather.SolarTime[i], day_status)
        end

        # Declare time series variables
        TemperatureAmbient = Weather.Temperature # [°F];
        TemperatureAmbientC = (TemperatureAmbient.-32).*(5/9) # [°C]
        PercentLighting = Schedules[:,1];
        PercentPlug = Schedules[:,2];
        PercentOccupied = Schedules[:,3];
        PVGeneration = Weather.PV;
        RadHeatGain = Weather.RadHeatGain;
        RadCooling = Weather.RadCool;
        CFM = Weather.CFM;
        InStorageBattery_1 = J_States[1] # [kWh]
        InStoragePCM_H_1 = J_States[2]
        InStoragePCM_C_1 = J_States[3]
        TemperatureIndoor_1 = J_States[4] # [°F]
    end
    return NumTime, TemperatureAmbient, TemperatureAmbientC, PercentLighting, PercentPlug, PercentOccupied, PVGeneration, RadHeatGain, RadCooling, CFM, InStorageBattery_1, InStoragePCM_H_1, InStoragePCM_C_1, TemperatureIndoor_1

end

