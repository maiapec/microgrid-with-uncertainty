# AA222 Project Model Predictive Control with uncertainty for Microgrid operation
# Complete Date: 06/01/2024

# **************
# Parameters to be chosen by the user:
begin
    simulated_years = 2014:2020 # Ex: 11 years are necessary to forecast for 10 years.
    Opt_Horizon = 168 # [hr]
    NumRun = 8760 # Can be adjusted or testing purposes, maximum 8760.
    # List of scenarios to be tested:
    namelist = ["Ground Truth","All Smooth Medium"]
    # Full name list would be: namelist = ["Ground Truth", "Only Smooth Temp","Only Smooth Wind","Only Smooth Dew Point","All Smooth Low","All Smooth Medium","All Smooth High","All Random Medium","All Noisy Medium"]
end
# **************


############ Initialize Tools ############
begin
    using JuMP
    using CSV
    using DataFrames
    # using PlotlyJS
    using Dates
    using XLSX
    using FileIO
    using Base
    using Random
    using Statistics
    using Gurobi
    using PyCall

    include("Util.jl")
    include("PassiveModel.jl")
    include("MPCAlgorithms.jl")
    include("ErrorFunctions.jl")
    include("MakeForecast.jl")
    include("MPCValidation.jl")
    include("MakeScenarios.jl")

    # Import necessary Python libraries
    pd = pyimport("pandas")

    # ENV["GRB_LICENSE_FILE"] = "C:\\Users\\Fred\\.julia\\environments\\v1.7\\gurobi.lic" # Fred's license
    # ENV["GRB_LICENSE_FILE"] = "/Library/gurobi1100/gurobi.lic" # Andreas' license
    ENV["GRB_LICENSE_FILE"] =  "/Library/gurobi1102/gurobi.lic" # Maia' license
end

############ Program Preparations ############
begin
    # Update automatically the date when this program is run.
    today_date = today()

    # Please update information of this program to automatically update the code name.
    version = 7.0

    code_name = "AA222_V$version._$today_date"

    # Create folder to later save data and plots
    begin
        # Define the path for the new folder
        folder_path = "Results/$code_name"

        # Use mkpath() to create the folder if it doesn't exist
        mkpath(folder_path)
    end

    # Set up logging
    log_file = "mpc_run_log.csv"
    if !isfile(log_file)
        open(log_file, "w") do file
            # Write header manually to avoid the error
            header = ["MPC", "Scenario", "Year", "RunTime", "ForecastHorizon", "LossOfLoad"] 
            # header = vcat(header, ["LossOfLoad_$i" for i in 1:1]...)
            write(file, join(header, ",") * "\n")
        end
    end
end

############ Load Data ############
begin
    # Load the weather data for Half Moon Bay from 2011 to 2020 (10 years of data)
    weather_directory = "NREL_NSRD_HalfMoonBay25_60Min" # hourly data
    weather_years = [] # Empty array of weather dataframes

    # Load the data from the weather files of all years
    for i in simulated_years
        file_path = joinpath(weather_directory, "137344_37.49_-122.42_$(i).csv")
        weather_file = DataFrame(CSV.File(file_path, header=3))
        push!(weather_years, weather_file) # Append the weather data to the array
    end
end

############ Declare Parameters ############ 

# Time Horizon Parameters
begin
    TimeStart = 1;
    TimeEnd = 8760;
    f_run = 1 # run frequency [hour/time]
    stepsize = 60*(1/f_run) # [min]
    δt = stepsize/60 # [hr] Declare stepzize for the optimization program
    # NumRun = (TimeEnd-TimeStart+1) # Max number of runs allowed
end   

# Further constant parameters
begin   
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

############ Generate Schedules ############
# Define load schedules for lighting, plugs, and occupancy
begin
    # Instructions
    begin
        # The are two schedules, format: Schedules[Lighting, Plugs, Occupancy]
        # All variables are unitless(0.5 = 50%). They are then used to multiply peak lighting/occupancy/plugs heat gain in the code.

        # SimpleSchedule:
        # About lighting schedule: full lighting intensity from 6 am to 10 pm daily
        # About plugs schedule: full plugs intensity from 6 am to 10 pm daily
        # About occupancy schedule: full occupancy from 6 am to 10 pm daily
    
        # ComplexSchedule:
        # About lighting schedule: full lighting intensity from 6 am to 10 pm daily
        # About plugs schedule: equipments (PC and phone) are charged from 7 pm to 12 am daily (5 hours)
        # About occupancy schedule: 50% occupany from 8 am to 8 pm, 100% occupancy during night time.
    end
    # Make SimpleSchedule
    # !!! Schedules have been modified to account for one more week at the end of the year
    begin
        SimpleSchedule = zeros(8760+24*7, 3);
        for i = 0:371
            for j = 10:17
                SimpleSchedule[i*24 + j, 1] = 1
                SimpleSchedule[i*24 + j, 2] = 1
                SimpleSchedule[i*24 + j, 3] = 1
            end
        end
    end
    # Make ComplexSchedule
    begin
        ComplexSchedule = zeros(8760+24*7, 3);
        for i = 0:371
            for j = 1:6
                ComplexSchedule[i*24 + j, 1] = 0
                ComplexSchedule[i*24 + j, 2] = 0
                ComplexSchedule[i*24 + j, 3] = 1
            end    
            for j = 7:8
                ComplexSchedule[i*24 + j, 1] = 1
                ComplexSchedule[i*24 + j, 2] = 0
                ComplexSchedule[i*24 + j, 3] = 1
            end  
            for j = 9:19
                ComplexSchedule[i*24 + j, 1] = 1
                ComplexSchedule[i*24 + j, 2] = 0
                ComplexSchedule[i*24 + j, 3] = 0.5
            end  
            for j = 20:20
                ComplexSchedule[i*24 + j, 1] = 1
                ComplexSchedule[i*24 + j, 2] = 1
                ComplexSchedule[i*24 + j, 3] = 0.5
            end  
            for j = 21:22
                ComplexSchedule[i*24 + j, 1] = 1
                ComplexSchedule[i*24 + j, 2] = 1
                ComplexSchedule[i*24 + j, 3] = 1
            end  
            for j = 23:24
                ComplexSchedule[i*24 + j, 1] = 0
                ComplexSchedule[i*24 + j, 2] = 1
                ComplexSchedule[i*24 + j, 3] = 1
            end  
        end
    end
end

############ Functions and Scenarios Choice ############
begin
    # Define MPC functions
    mpc_functions = Dict(
        "MPCDirect" => (Optimize_MPCDirect_full, ()), # MPC Direct which does not require more parameters
        "MPCSmooth" => (Optimize_MPCSmooth_full, (1,0.1)), # Parameters for hours penalized and penalty set on deviation of actions
    )
    weather_scenarios = get_specific_scenarios(namelist)
    # To get the complete list of scenarios, replace the previous two lines with the following:
    # weather_scenarios = get_all_scenarios()

    Schedule = SimpleSchedule
end

############## Run cross-validation ##############
println("------------- New Run ----------------\n")
cross_validate_mpc(mpc_functions, weather_scenarios, simulated_years, weather_years, TimeStart, NumRun, Opt_Horizon, Schedule, log_file)
println("\n------------- End of Run -------------")