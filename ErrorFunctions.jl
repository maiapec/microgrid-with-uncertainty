begin 
    using Random
    using PlotlyJS
    using DataFrames
    using CSV
    using Dates
    using Statistics
    using PyCall
    using Base
    include("Util.jl")
end

# 1. Function to quantify the variability for each weather variable (except PV variables)

function make_std_schedule(directory_weather_data="NREL_NSRD_HalfMoonBay25_60Min/")
    """
    Computes the standard deviation across all years of the variables (except PV variables), for each hour.
    Takes as input the directory of the folder containing all the data files.
    Returns a DataFrame with the smoothed standard deviation of the variables for each hour of one year+one week 
    (to be able to make forecast fot the whole year).
    """
    
    df = []
    pd = pyimport("pandas")

    # Load the data from the weather files of all years
    for i in 1998:2022
        file_path = joinpath(directory_weather_data, "137344_37.49_-122.42_$(i).csv")
        weather = DataFrame(CSV.File(file_path, header=3))

        # Create the Datetime column
        weather.datetime = [create_datetime_from_row(row) for row in eachrow(weather)]

        # Convert Julia DataFrame to Pandas DataFrame
        columns = names(weather)
        data = [getproperty(weather, col) for col in columns]
        pandas_df = pd.DataFrame(Dict(zip(columns, data)))

         # Ensure the datetime column is in the correct format and set as the index
        pandas_df["datetime"] = pd.to_datetime(pandas_df["datetime"])
        pandas_df = pandas_df.set_index("datetime")
        push!(df, pandas_df)
    end

    # Combine the data frames
    combined_df = pd.concat(df)
    
    # Calculate the standard deviation of the variables for each hour of the year
    variables = ["Temperature", "Relative Humidity", "Wind Speed", "Dew Point", "Cloud Type"]

    # Repeat hours from 1 to 8760, and then again, in total 25 times
    hour_of_year = repeat(1:8760, 25)
    # One value per hour of the year
    stdeviations = Dict(var => combined_df[var].groupby(hour_of_year).std() for var in variables)
    # Smoothe values
    smooth_stdeviations = Dict(var => stdeviations[var].rolling(window=2000, min_periods=1).mean() for var in variables)
    # Smoothe to have one value per week of the year
    week = 24 * 7
    week_of_year = repeat(1:8760) .÷ week .+ 1
    stdeviations_per_week = Dict(var => smooth_stdeviations[var].groupby(week_of_year).mean() for var in variables)

    # stdeviations_per_week has the values for each week.
    # Repeat values for all hours of one year
    stdeviations_allyear = Dict(var => repeat(collect(stdeviations_per_week[var]), inner=[week])[1:8760] for var in variables)

    # Add again the first week in the end to have all year + 1 week.
    stdeviations_allyear = Dict(var => vcat(stdeviations_allyear[var], stdeviations_allyear[var][1:week]) for var in variables)

    # Convert to DataFrame
    stdeviations_allyear = DataFrame(stdeviations_allyear)

    return stdeviations_allyear
end

# 2. Methods to simulate forecast data. They can be combined.

function Add_smooth_uncertainty_to_variable(var_actual, var_std_schedule, alpha = 0.9, beta=0.1, step_smoothing=5, isRH=false, isWS=false)
    """ Simulates a smooth forecast of var_actual based on the variability schedule var_std_schedule for this variable.
        NB: var_actual and var_std should have the same indexes.
        alpha = how much of the previous value simulated is kept in the forecast. Use a high value of alpha.
        beta = how much of the variability schedule is added to the forecast
        step_smoothing = the size of the moving average window used to smooth the forecast
        isRH = if true, the forecast is limited to 0-100%.
        isWS = if true, the forecast is limited to positive values.
    """
    l = length(var_actual)
    var_sim = similar(var_actual, Float64)
    var_sim[1] = var_actual[1]
    delta = [i / l for i in 1:l]
    for i in 2:l
        incre = (var_actual[i] - var_actual[i-1]) * (1 + beta*var_std_schedule[i]*rand([-1, 1]) * delta[i])
        var_update = alpha * var_sim[i-1] + (1 - alpha) * var_actual[i-1]
        if isRH
            incre = min(incre, 100 - var_update)
            incre = max(incre, 0 - var_update)
        end
        if isWS
            incre = max(incre, 0 - var_update)
        end
        var_sim[i] = var_update + incre
    end
    var_sim = Smoothing(var_sim, step_smoothing)
    # Keep only two decimals
    var_sim = round.(var_sim, digits=2)
    return var_sim
end

function Add_smooth_uncertainty_to_variable_only_plus(var_actual, var_std_schedule, alpha = 0.9, beta=0.1, step_smoothing=5, isRH=false, isWS=false)
    """ Simulates a smooth forecast of var_actual based on the variability schedule var_std_schedule for this variable.
        Systematically overforecasts. Note that for a set of parameters, the output is fully determined (no more randomness).

        NB: var_actual and var_std should have the same indexes.
        alpha = how much of the previous value simulated is kept in the forecast. Use a high value of alpha.
        beta = how much of the variability schedule is added to the forecast
        step_smoothing = the size of the moving average window used to smooth the forecast
        isRH = if true, the forecast is limited to 0-100%.
        isWS = if true, the forecast is limited to positive values.
    """
    l = length(var_actual)
    var_sim = similar(var_actual, Float64)
    var_sim[1] = var_actual[1]
    delta = [i / l for i in 1:l]
    for i in 2:l
        incre = (var_actual[i] - var_actual[i-1]) * (1 + beta*var_std_schedule[i]* 1 * delta[i])
        var_update = alpha * var_sim[i-1] + (1 - alpha) * var_actual[i-1]
        if isRH
            incre = min(incre, 100 - var_update)
            incre = max(incre, 0 - var_update)
        end
        if isWS
            incre = max(incre, 0 - var_update)
        end
        var_sim[i] = var_update + incre
    end
    var_sim = Smoothing(var_sim, step_smoothing)
    # Keep only two decimals
    var_sim = round.(var_sim, digits=2)
    return var_sim
end

function Add_smooth_uncertainty_to_variable_only_minus(var_actual, var_std_schedule, alpha = 0.9, beta=0.1, step_smoothing=5, isRH=false, isWS=false)
    """ Simulates a smooth forecast of var_actual based on the variability schedule var_std_schedule for this variable.
        Systematically underforecasts. Note that for a set of parameters, the output is fully determined (no more randomness).

        NB: var_actual and var_std should have the same indexes.
        alpha = how much of the previous value simulated is kept in the forecast. Use a high value of alpha.
        beta = how much of the variability schedule is added to the forecast
        step_smoothing = the size of the moving average window used to smooth the forecast
        isRH = if true, the forecast is limited to 0-100%.
    """
    l = length(var_actual)
    var_sim = similar(var_actual, Float64)
    var_sim[1] = var_actual[1]
    delta = [i / l for i in 1:l]
    for i in 2:l
        incre = (var_actual[i] - var_actual[i-1]) * (1 + beta*var_std_schedule[i]* (-1) * delta[i])
        var_update = alpha * var_sim[i-1] + (1 - alpha) * var_actual[i-1]
        if isRH
            incre = min(incre, 100 - var_update)
            incre = max(incre, 0 - var_update)
        end
        if isWS
            incre = max(incre, 0 - var_update)
        end
        var_sim[i] = var_update + incre
    end
    var_sim = Smoothing(var_sim, step_smoothing)
    # Keep only two decimals
    var_sim = round.(var_sim, digits=2)
    return var_sim
end

function Add_random_uncertain_points_to_variable(var_actual, var_std_schedule, freq=5, beta=0.3, isRH=false, isWS=false)
    """ Simulates a forecast of var_actual which simply adds random events to var_actual at a frequency approximately equal to freq.
        NB: var_actual and var_std should have the same indexes.
        freq = approximately delta time between two random events added artificially.
        beta = how much of variability is added to the forecast.
        isRH = if true, the forecast is limited to 0-100%.
        Random events do not increase with time.
    """
    l = length(var_actual)
    var_sim = similar(var_actual, Float64)
    var_sim[1] = var_actual[1]
    for i in 2:l
        var_sim[i] = var_actual[i]
        if rand(1:freq) == 1
            if isRH
                var_sim[i] = min(100, max(0, var_sim[i] + beta*var_std_schedule[i]*rand([-1, 1])))
            else
                var_sim[i] += beta*var_std_schedule[i]*rand([-1, 1])
            end
            if isWS
                var_sim[i] = max(0, var_sim[i])
            end
        end
    end
    # Keep only two decimals
    var_sim = round.(var_sim, digits=2)
    return var_sim
end

function Add_noise_to_variable(var_actual, var_std_schedule, alpha=0.5, beta=0.3, isRH=false, isWS=false)
    """ Simulates a forecast of var_actual which adds noise to the actual values.
        NB: var_actual and var_std should have the same indexes.
        alpha = how much of the previous value simulated is kept in the forecast. 
        beta = how much of the variability schedule is added to the forecast. The value should be taken low.
        isRH = if true, the forecast is limited to 0-100%.
        NB: This function is equivalent to Add_random_uncertain_points_to_variable with freq=1.
        
        !!! Choose a relatively low value of alpha to ensure we're still following the actual data.
        !!! Choose a relatively low value of beta to ensure we're not adding too much noise.
    """
    l = length(var_actual)
    var_sim = similar(var_actual, Float64)
    var_sim[1] = var_actual[1]
    for i in 2:l
        incre = beta*var_std_schedule[i]*rand([-1, 1])
        var_update = alpha * var_sim[i-1] + (1 - alpha) * var_actual[i-1]
        if isRH
            incre = min(incre, 100 - var_update)
            incre = max(incre, 0 - var_update)
        end
        if isWS
            incre = max(incre, 0 - var_update)
        end
        var_sim[i] = var_update + incre
    end
    # Keep only two decimals
    var_sim = round.(var_sim, digits=2)
    return var_sim
end
