######## MPC Algorithms ########

# Baseline/Control Algorithm

function Optimize_Baseline(Weather, Schedules, J_States)

    ########## Data Preparation  ##########  
    begin
        # Add PV output first
        horizon0 = 8760
        Weather = Make_PV_output(Weather, noct_installed, module_height, wind_height, module_emissivity, module_absorption, module_surface_tilt, module_width, module_length, Berg_tilt, horizon0)
        # Change to °F
        Weather.Temperature .= (Weather.Temperature .* (9/5)) .+ 32 # [°F]
        # Then prepare the rest of the data
        NumTime, TemperatureAmbient, TemperatureAmbientC, PercentLighting, PercentPlug, PercentOccupied, PVGeneration, RadHeatGain, RadCooling, CFM, InStorageBattery_1, InStoragePCM_H_1, InStoragePCM_C_1, TemperatureIndoor_1 = DataPreparation(Weather, Schedules, J_States)
    end

    ########## Declare model  ##########
    begin
        # Define the model name and solver. In this case, model name is "m"
        # Set path to license (for those using Gurobi)
        m = Model(Gurobi.Optimizer)
    end

    ######## Decision variables ########
    begin
        @variable(m, PV2H[1:NumTime] >= 0); # [kW] electrical power transfer from PV to home (Berg)

        @variable(m, PV2G[1:NumTime] >= 0); # [kW] electrical power transfer from PV to ground (curtailment)

        @variable(m, PV2B[1:NumTime] >= 0); # [kW] electrical power transfer from PV to battery

        @variable(m, B2H[1:NumTime] >= 0); # [kW] electrical power transfer from battery to home (Berg)

        @variable(m, H2HP[1:NumTime] >= 0); # [kW] electrical power transfer from home (Berg) to heat pump heating mode

        @variable(m, HP2H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from heat pump heating mode to home (Berg)

        @variable(m, H2C[1:NumTime] >= 0); # [kW] electrical power transfer from home (Berg) to heat pump cooling unit

        @variable(m, C2H[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from heat pump cooling mode to home (Berg)
        
        @variable(m, HP_OPTION[1:NumTime], Bin) # Heating Mode = 1, Cooling Mode = 0

        @variable(m, TES_Charging_OPTION[1:NumTime], Bin) # Charge hot TES = 1, Charge cold TES = 0

        @variable(m, TES_Discharging_OPTION[1:NumTime], Bin) # Discharge hot TES = 1, Discharge cold TES = 0

        @variable(m, TES_CD_OPTION[1:NumTime], Bin) # Charge TES = 1, Discharge TES = 0

        @variable(m, B_OPTION[1:NumTime], Bin) # Charge battery = 1, Discharge battery = 0

        @variable(m, HP2PCM_H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from heat pump heating mode to PCM heating storage

        @variable(m, C2PCM_C[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from heat pump cooling mode to PCM cooling storage

        @variable(m, PCM_H2H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from PCM heating storage to home (Berg)

        @variable(m, PCM_C2H[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from PCM cooling storage to home (Berg)

        @variable(m, InStorageBattery[1:NumTime] >= 0); # [kWh] Battery Remaining Charge

        @variable(m, InStoragePCM_H[1:NumTime] >= 0); # [BTU] PCM Heating Remaining Charge

        @variable(m, InStoragePCM_C[1:NumTime] >= 0); # [BTU] PCM Cooling Remaining Charge

        @variable(m, TemperatureIndoor[1:NumTime] >= 0); # [°F] Indoor Air Temperature

        @variable(m, G2H[1:NumTime] >= 0); # [kW] Loss of Load
    end

    ############ Objective Functions #############
    begin
        # Penalty for loss of load [$]
        @expression(m, penalty_l, δt * sum(G2H[t] for t = 1:NumTime))

        # Total Cost over Optimization Horizon [$]
        @objective(m, Min, penalty_l);
    end

    ############# Expressions ############
    begin
        # DeltaTemp 
        @expression(m, TempDelta[t=1:NumTime], TemperatureAmbient[t] - TemperatureIndoor[t]); # [°F]

        # Electricity usage from lighting 
        @expression(m, E_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t]); # [kW]

        # Electricity usage from plugs
        @expression(m, E_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t]); # [kW]

        # Total electricity usage
        @expression(m, E_total[t=1:NumTime], E_Lighting[t] + E_Plugs[t]); # [kW]

        # Calculate Ventilation
        @expression(m, CFMVen[t=1:NumTime], min(Rp * PercentOccupied[t] * MaxOccupancy + Ra * Area, PercentOccupied[t] * MaxOccupancy * Ventilation)); # [ft^3/min]

        # Heat gain from lighting
        @expression(m, Q_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t] * 3412.14); # [BTU/hr]

        # Heat gain from plugs
        @expression(m, Q_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t] * 3412.14); # [BTU/hr]

        # Heat gain from occupancy
        @expression(m, Q_Occupancy[t=1:NumTime], MaxOccupancy * PercentOccupied[t] * TotalPersonHeat); # [BTU/hr]

        # Heat gain from infiltration
        # 1.08 = specific heat capacity of air at STP:0.24 [BTU/(lb*°F)] * air density at STP:0.075 [lb/ft^3] * 60 [min/hr] 
        @expression(m, Q_Infiltration[t=1:NumTime], 1.08 * TempDelta[t] * CFM[t]); # [BTU/hr]

        # Heat gain from ventilation
        @expression(m, Q_Ventilation[t=1:NumTime], 1.08 * TempDelta[t] * CFMVen[t]); # [BTU/hr]

        # Heat gain from lighting, plugs, occupancy, ventilation, infiltration
        @expression(m, Q_Others[t=1:NumTime], Q_Infiltration[t] + Q_Ventilation[t] + Q_Occupancy[t] + Q_Lighting[t] + Q_Plugs[t]); # [BTU/hr]

        # Heat gain through structural evenlope
        @expression(m, Q_Envelope[t=1:NumTime], UA * TempDelta[t]); # [BTU/hr]

        # Heat gain from solar radiation
        @expression(m, Q_Rad[t=1:NumTime], SHGC * RadHeatGain[t]); # [BTU/hr]

        # Radiative cooling
        @expression(m, Q_RadCool[t=1:NumTime], RCC * RadCooling[t]); # [BTU/hr]

        # Detailed COP of heating and cooling (linearized using TemperatureIndoor = 22 [°C])
        # COP (Heating Home)
        @expression(m, COP_Heating_H[t=1:NumTime], HP_a + HP_b * (T_indoor_constant - TemperatureAmbientC[t]) + HP_c * (T_indoor_constant - TemperatureAmbientC[t])^2);
        
        # COP (Heating PCM H)
        @expression(m, COP_Heating_PCM[t=1:NumTime], HP_a + HP_b * (48 - TemperatureAmbientC[t]) + HP_c * (48 -TemperatureAmbientC[t])^2); 

        # COP (Cooling Home)
        @expression(m, COP_Cooling_H[t=1:NumTime], HP_a + HP_b * (TemperatureAmbientC[t] - T_indoor_constant) + HP_c * (TemperatureAmbientC[t] - T_indoor_constant)^2);
        
        # COP (Cooling PCM C)
        @expression(m, COP_Cooling_PCM[t=1:NumTime], HP_a + HP_b * (TemperatureAmbientC[t] - 11) + HP_c * (TemperatureAmbientC[t] - 11)^2); 
        
    end

    ############# Constraints ############
    begin
        # Initialization Constraints with States from last timestep
        begin
            # Set point temperature range constraints
            @constraint(m, TemperatureIndoor[1] == TemperatureIndoor_1); # [°F]

            # PCM heating storage initialization constraint, node at PCM heating storage
            @constraint(m, InStoragePCM_H[1] == InStoragePCM_H_1); # [BTU]

            # PCM cooling storage initialization constraint, node at PCM cooling storage
            @constraint(m, InStoragePCM_C[1] == InStoragePCM_C_1); # [BTU]

            # Battery storage initialization constraint
            @constraint(m, InStorageBattery[1] == InStorageBattery_1); # [kWh] 
        end

        # Set point temperature range constraints 
        @constraint(m, [t=1:NumTime], TemperatureIndoor[t] <= SetPointT_High); # [°F]

        @constraint(m, [t=1:NumTime], TemperatureIndoor[t] >= SetPointT_Low); # [°F]

        # Internal temperature balance evolution constraint
        @constraint(m, [t=1:NumTime-1], TemperatureIndoor[t+1] == TemperatureIndoor[t] + 
        (δt/TC)*(Q_Others[t] + Q_Envelope[t] + Q_Rad[t] - Q_RadCool[t] + HP2H[t] - C2H[t] + PCM_H2H[t] - PCM_C2H[t])); # [°F]
        
        # PV energy balance constraint, node at PV
        @constraint(m, [t=1:NumTime], PVGeneration[t] * PVSize ==  PV2B[t] + PV2H[t] + PV2G[t]); # [kW]
    
        # House electricity load constraint, node at house, battery efficiency modeled
        @constraint(m, [t=1:NumTime], E_total[t] + H2HP[t] + H2C[t] == PV2H[t] * η_PVIV + B2H[t] * η + G2H[t]); # [kW]

        # Battery storage balance constraint, node at battery, battery leakage modeled, battery efficiency modeled
        @constraint(m, [t=1:NumTime-1], InStorageBattery[t+1] == InStorageBattery[t] * δt * (1 - BatteryLoss) +  δt * (PV2B[t] * η - B2H[t])); # [kWh]

        # Battery discharging constraint, node at battery 
        @constraint(m, [t=1:NumTime], δt * B2H[t] <= InStorageBattery[t]); # [kWh]
        
        # Battery power inverter constraint, node at battery (inverter power constraint) 
        @constraint(m, [t=1:NumTime], B2H[t] + PV2B[t] <= InverterSize); # [kW]
        
        # Battery storage size constraint, node at battery
        @constraint(m, [t=1:NumTime], InStorageBattery[t] <= BatterySize); # [kWh]
        
        # Battery storage max discharge constraint, node at battery, always at least 20% full 
        @constraint(m, [t=1:NumTime], InStorageBattery[t] >= BatterySize * (1-MaxDischarge)); # [kWh]
        
        # Heating and Cooling Constraints (Needs Remodeling)
        begin
            # Detailed COP of heating and cooling (requires nonlinear optimization solver)
            # Heating (from kW to BTU/hr), node at heat pump heating option
            @constraint(m, [t=1:NumTime], H2HP[t] * 3412.14 == HP2PCM_H[t]/COP_Heating_PCM[t] + HP2H[t]/COP_Heating_H[t]);
            
            # Cooling (from kWh to BTU), node at heat pump cooling option
            @constraint(m, [t=1:NumTime], H2C[t] * 3412.14 == C2PCM_C[t]/COP_Cooling_PCM[t] + C2H[t]/COP_Cooling_H[t]);

            #=
            # Heating (from kW to BTU/hr), node at heat pump heating option
            @constraint(m, [t=1:NumTime], H2HP[t] * 3412.14 == (HP2PCM_H[t] + HP2H[t])/COP_H); # [BTU/hr]

            # Cooling (from kW to BTU/hr), node at heat pump cooling option
            @constraint(m, [t=1:NumTime], H2C[t] * 3412.14 == (C2PCM_C[t] + C2H[t])/COP_C); # [BTU/hr]
            =#

            # Heat pump power capacity constraint, node at heat pump heating mode
            @constraint(m, [t=1:NumTime], H2HP[t] <= Cap_HP_H); # [kW]

            # Heat pump power capacity constraint, node at heat pump cooling mode
            @constraint(m, [t=1:NumTime], H2C[t] <= Cap_HP_C); # [kW]

            # PCM heating storage balance constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime-1], InStoragePCM_H[t+1] == InStoragePCM_H[t] + δt * (HP2PCM_H[t] - PCM_H2H[t])); # [BTU]
            
            # PCM cooling storage balance constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime-1], InStoragePCM_C[t+1] == InStoragePCM_C[t] + δt * (C2PCM_C[t] - PCM_C2H[t])); # [BTU]
            
            # PCM heating storage discharging constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime], δt * PCM_H2H[t] <= InStoragePCM_H[t]); # [BTU]
            
            # PCM cooling storage discharging constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime], δt * PCM_C2H[t] <= InStoragePCM_C[t]); # [BTU]
            
            # PCM heating storage size constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime], InStoragePCM_H[t] <= PCM_H_Size); # [BTU]

            # PCM cooling storage size constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime], InStoragePCM_C[t] <= PCM_C_Size); # [BTU]
        end
        # Single Operation Choice Constraints
        begin
            # If HP_OPTION = 1, heating mode on, heating only
            @constraint(m, [t=1:NumTime], H2HP[t] <= M * HP_OPTION[t]) # [kW]          

            # If HP_OPTION = 0, cooling mode on, cooling only
            @constraint(m, [t=1:NumTime], H2C[t] <= M * (1 - HP_OPTION[t])) # [kW]    

            # If TES_Charging_OPTION = 1, charge hot TES only
            @constraint(m, [t=1:NumTime], HP2PCM_H[t] <= M * TES_Charging_OPTION[t]) # [BTU/hr]     
    
            # If TES_Charging_OPTION = 0, charge cold TES only
            @constraint(m, [t=1:NumTime], C2PCM_C[t] <= M * (1 - TES_Charging_OPTION[t])) # [BTU/hr]

            # If TES_Discharging_OPTION = 1, discharge hot TES only
            @constraint(m, [t=1:NumTime], PCM_H2H[t] <= M * TES_Discharging_OPTION[t]) # [BTU/hr]          
    
            # If TES_Discharging_OPTION = 0, discharge cold TES only
            @constraint(m, [t=1:NumTime], PCM_C2H[t] <= M * (1 - TES_Discharging_OPTION[t])) # [BTU/hr]

            # Total TES Charging 
            @expression(m, TES_Charge[t=1:NumTime], HP2PCM_H[t] + C2PCM_C[t]); # [BTU/hr]

            # Total TES Discharging 
            @expression(m, TES_Discharge[t=1:NumTime], PCM_H2H[t] + PCM_C2H[t]); # [BTU/hr]

            # If TES_CD_OPTION = 1, charge only
            @constraint(m, [t=1:NumTime], TES_Charge[t] <= M * TES_CD_OPTION[t]) # [BTU/hr]          

            # If TES_CD_OPTION = 0, discharge only
            @constraint(m, [t=1:NumTime], TES_Discharge[t] <= M * (1 - TES_CD_OPTION[t])) # [BTU/hr] 

            # If B_OPTION = 1, charge only
            @constraint(m, [t=1:NumTime], PV2B[t] <= M * B_OPTION[t]) # [kW]          

            # If B_OPTION = 0, discharge only
            @constraint(m, [t=1:NumTime], B2H[t] <= M * (1 - B_OPTION[t])) # [kW] 
        end
    end

    # display(latex_formulation(m)) # Display the model in latex formulation

    ########### Solve  ##########
    optimize!(m);

    ########### Model Results  ##########
    LossOfLoad = δt * value.(G2H)

    # MAYBE RETURN SOMETHING ELSE
    return LossOfLoad
end

# Direct MPC Algorithm

function Optimize_MPCDirect_onestep(Weather, Schedules, J_States)
    """ One step of the Direct MPC algorithm.
        Inputs: 
        - Weather for one week (or the chosen time horizon), which already contains PV and where temperatures have already been converted to °F 
        (typically one element of the output of of Make_weather_forecast_full_year).
        - Schedule, either SimpleSchedule or ComplexSchedule, over the time horizon.
        - J_States, the states of the system before performing this step.
    """
    ########## Data preparation  ##########
    begin
        # Directly use the function DataPreparation
        NumTime, TemperatureAmbient, TemperatureAmbientC, PercentLighting, PercentPlug, PercentOccupied, PVGeneration, RadHeatGain, RadCooling, CFM, InStorageBattery_1, InStoragePCM_H_1, InStoragePCM_C_1, TemperatureIndoor_1 = DataPreparation(Weather, Schedules, J_States)
    end

    ########## Declare model  ##########
    begin
        # Define the model name and solver. In this case, model name is "m".
        # Path to license must have been defined.
        m = Model(Gurobi.Optimizer)
    end
    ######## Decision variables ########
    begin
        @variable(m, PV2H[1:NumTime] >= 0); # [kW] electrical power transfer from PV to home (Berg)

        @variable(m, PV2G[1:NumTime] >= 0); # [kW] electrical power transfer from PV to ground (curtailment)

        @variable(m, PV2B[1:NumTime] >= 0); # [kW] electrical power transfer from PV to battery

        @variable(m, B2H[1:NumTime] >= 0); # [kW] electrical power transfer from battery to home (Berg)

        @variable(m, H2HP[1:NumTime] >= 0); # [kW] electrical power transfer from home (Berg) to heat pump heating mode

        @variable(m, HP2H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from heat pump heating mode to home (Berg)

        @variable(m, H2C[1:NumTime] >= 0); # [kW] electrical power transfer from home (Berg) to heat pump cooling unit

        @variable(m, C2H[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from heat pump cooling mode to home (Berg)
        
        @variable(m, HP_OPTION[1:NumTime], Bin) # Heating Mode = 1, Cooling Mode = 0

        @variable(m, TES_Charging_OPTION[1:NumTime], Bin) # Charge hot TES = 1, Charge cold TES = 0

        @variable(m, TES_Discharging_OPTION[1:NumTime], Bin) # Discharge hot TES = 1, Discharge cold TES = 0

        @variable(m, TES_CD_OPTION[1:NumTime], Bin) # Charge TES = 1, Discharge TES = 0

        @variable(m, B_OPTION[1:NumTime], Bin) # Charge battery = 1, Discharge battery = 0

        @variable(m, HP2PCM_H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from heat pump heating mode to PCM heating storage

        @variable(m, C2PCM_C[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from heat pump cooling mode to PCM cooling storage

        @variable(m, PCM_H2H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from PCM heating storage to home (Berg)

        @variable(m, PCM_C2H[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from PCM cooling storage to home (Berg)

        @variable(m, InStorageBattery[1:NumTime] >= 0); # [kWh] Battery Remaining Charge

        @variable(m, InStoragePCM_H[1:NumTime] >= 0); # [BTU] PCM Heating Remaining Charge

        @variable(m, InStoragePCM_C[1:NumTime] >= 0); # [BTU] PCM Cooling Remaining Charge

        @variable(m, TemperatureIndoor[1:NumTime] >= 0); # [°F] Indoor Air Temperature

        @variable(m, G2H[1:NumTime] >= 0); # [kW] Loss of Load
    end
    ############ Objective Functions #############
    begin
        # Penalty for loss of load [$]
        @expression(m, penalty_l, δt * sum(G2H[t] for t = 1:NumTime))

        # Total Cost over Optimization Horizon [$]
        @objective(m, Min, penalty_l);
    end
    ############# Expressions ############
    begin
        # DeltaTemp 
        @expression(m, TempDelta[t=1:NumTime], TemperatureAmbient[t] - TemperatureIndoor[t]); # [°F]

        # Electricity usage from lighting 
        @expression(m, E_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t]); # [kW]

        # Electricity usage from plugs
        @expression(m, E_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t]); # [kW]

        # Total electricity usage
        @expression(m, E_total[t=1:NumTime], E_Lighting[t] + E_Plugs[t]); # [kW]

        # Calculate Ventilation
        @expression(m, CFMVen[t=1:NumTime], min(Rp * PercentOccupied[t] * MaxOccupancy + Ra * Area, PercentOccupied[t] * MaxOccupancy * Ventilation)); # [ft^3/min]

        # Heat gain from lighting
        @expression(m, Q_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t] * 3412.14); # [BTU/hr]

        # Heat gain from plugs
        @expression(m, Q_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t] * 3412.14); # [BTU/hr]

        # Heat gain from occupancy
        @expression(m, Q_Occupancy[t=1:NumTime], MaxOccupancy * PercentOccupied[t] * TotalPersonHeat); # [BTU/hr]

        # Heat gain from infiltration
        # 1.08 = specific heat capacity of air at STP:0.24 [BTU/(lb*°F)] * air density at STP:0.075 [lb/ft^3] * 60 [min/hr] 
        @expression(m, Q_Infiltration[t=1:NumTime], 1.08 * TempDelta[t] * CFM[t]); # [BTU/hr]

        # Heat gain from ventilation
        @expression(m, Q_Ventilation[t=1:NumTime], 1.08 * TempDelta[t] * CFMVen[t]); # [BTU/hr]

        # Heat gain from lighting, plugs, occupancy, ventilation, infiltration
        @expression(m, Q_Others[t=1:NumTime], Q_Infiltration[t] + Q_Ventilation[t] + Q_Occupancy[t] + Q_Lighting[t] + Q_Plugs[t]); # [BTU/hr]

        # Heat gain through structural evenlope
        @expression(m, Q_Envelope[t=1:NumTime], UA * TempDelta[t]); # [BTU/hr]

        # Heat gain from solar radiation
        @expression(m, Q_Rad[t=1:NumTime], SHGC * RadHeatGain[t]); # [BTU/hr]

        # Radiative cooling
        @expression(m, Q_RadCool[t=1:NumTime], RCC * RadCooling[t]); # [BTU/hr]

        # Detailed COP of heating and cooling (linearized using TemperatureIndoor = 22 [°C])
        # COP (Heating Home)
        @expression(m, COP_Heating_H[t=1:NumTime], HP_a + HP_b * (T_indoor_constant - TemperatureAmbientC[t]) + HP_c * (T_indoor_constant - TemperatureAmbientC[t])^2);
        
        # COP (Heating PCM H)
        @expression(m, COP_Heating_PCM[t=1:NumTime], HP_a + HP_b * (48 - TemperatureAmbientC[t]) + HP_c * (48 -TemperatureAmbientC[t])^2); 

        # COP (Cooling Home)
        @expression(m, COP_Cooling_H[t=1:NumTime], HP_a + HP_b * (TemperatureAmbientC[t] - T_indoor_constant) + HP_c * (TemperatureAmbientC[t] - T_indoor_constant)^2);
        
        # COP (Cooling PCM C)
        @expression(m, COP_Cooling_PCM[t=1:NumTime], HP_a + HP_b * (TemperatureAmbientC[t] - 11) + HP_c * (TemperatureAmbientC[t] - 11)^2); 
        
    end
    ############# Constraints ############
    begin
        # Initialization Constraints with States from last timestep
        begin
            # Set point temperature range constraints
            @constraint(m, TemperatureIndoor[1] == TemperatureIndoor_1); # [°F]

            # PCM heating storage initialization constraint, node at PCM heating storage
            @constraint(m, InStoragePCM_H[1] == InStoragePCM_H_1); # [BTU]

            # PCM cooling storage initialization constraint, node at PCM cooling storage
            @constraint(m, InStoragePCM_C[1] == InStoragePCM_C_1); # [BTU]

            # Battery storage initialization constraint
            @constraint(m, InStorageBattery[1] == InStorageBattery_1); # [kWh] 
        end

        # Set point temperature range constraints 
        @constraint(m, [t=1:NumTime], TemperatureIndoor[t] <= SetPointT_High); # [°F]

        @constraint(m, [t=1:NumTime], TemperatureIndoor[t] >= SetPointT_Low); # [°F]

        # Internal temperature balance evolution constraint
        @constraint(m, [t=1:NumTime-1], TemperatureIndoor[t+1] == TemperatureIndoor[t] + 
        (δt/TC)*(Q_Others[t] + Q_Envelope[t] + Q_Rad[t] - Q_RadCool[t] + HP2H[t] - C2H[t] + PCM_H2H[t] - PCM_C2H[t])); # [°F]
        
        # PV energy balance constraint, node at PV
        @constraint(m, [t=1:NumTime], PVGeneration[t] * PVSize ==  PV2B[t] + PV2H[t] + PV2G[t]); # [kW]
    
        # House electricity load constraint, node at house, battery efficiency modeled
        @constraint(m, [t=1:NumTime], E_total[t] + H2HP[t] + H2C[t] == PV2H[t] * η_PVIV + B2H[t] * η + G2H[t]); # [kW]

        # Battery storage balance constraint, node at battery, battery leakage modeled, battery efficiency modeled
        @constraint(m, [t=1:NumTime-1], InStorageBattery[t+1] == InStorageBattery[t] * δt * (1 - BatteryLoss) +  δt * (PV2B[t] * η - B2H[t])); # [kWh]

        # Battery discharging constraint, node at battery 
        @constraint(m, [t=1:NumTime], δt * B2H[t] <= InStorageBattery[t]); # [kWh]
        
        # Battery power inverter constraint, node at battery (inverter power constraint) 
        @constraint(m, [t=1:NumTime], B2H[t] + PV2B[t] <= InverterSize); # [kW]
        
        # Battery storage size constraint, node at battery
        @constraint(m, [t=1:NumTime], InStorageBattery[t] <= BatterySize); # [kWh]
        
        # Battery storage max discharge constraint, node at battery, always at least 20% full 
        @constraint(m, [t=1:NumTime], InStorageBattery[t] >= BatterySize * (1-MaxDischarge)); # [kWh]
        
        # Heating and Cooling Constraints (Needs Remodeling)
        begin
            # Detailed COP of heating and cooling (requires nonlinear optimization solver)
            # Heating (from kW to BTU/hr), node at heat pump heating option
            @constraint(m, [t=1:NumTime], H2HP[t] * 3412.14 == HP2PCM_H[t]/COP_Heating_PCM[t] + HP2H[t]/COP_Heating_H[t]);
            
            # Cooling (from kWh to BTU), node at heat pump cooling option
            @constraint(m, [t=1:NumTime], H2C[t] * 3412.14 == C2PCM_C[t]/COP_Cooling_PCM[t] + C2H[t]/COP_Cooling_H[t]);

            #=
            # Heating (from kW to BTU/hr), node at heat pump heating option
            @constraint(m, [t=1:NumTime], H2HP[t] * 3412.14 == (HP2PCM_H[t] + HP2H[t])/COP_H); # [BTU/hr]

            # Cooling (from kW to BTU/hr), node at heat pump cooling option
            @constraint(m, [t=1:NumTime], H2C[t] * 3412.14 == (C2PCM_C[t] + C2H[t])/COP_C); # [BTU/hr]
            =#

            # Heat pump power capacity constraint, node at heat pump heating mode
            @constraint(m, [t=1:NumTime], H2HP[t] <= Cap_HP_H); # [kW]

            # Heat pump power capacity constraint, node at heat pump cooling mode
            @constraint(m, [t=1:NumTime], H2C[t] <= Cap_HP_C); # [kW]

            # PCM heating storage balance constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime-1], InStoragePCM_H[t+1] == InStoragePCM_H[t] + δt * (HP2PCM_H[t] - PCM_H2H[t])); # [BTU]
            
            # PCM cooling storage balance constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime-1], InStoragePCM_C[t+1] == InStoragePCM_C[t] + δt * (C2PCM_C[t] - PCM_C2H[t])); # [BTU]
            
            # PCM heating storage discharging constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime], δt * PCM_H2H[t] <= InStoragePCM_H[t]); # [BTU]
            
            # PCM cooling storage discharging constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime], δt * PCM_C2H[t] <= InStoragePCM_C[t]); # [BTU]
            
            # PCM heating storage size constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime], InStoragePCM_H[t] <= PCM_H_Size); # [BTU]

            # PCM cooling storage size constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime], InStoragePCM_C[t] <= PCM_C_Size); # [BTU]
        end
        # Single Operation Choice Constraints
        begin
            # If HP_OPTION = 1, heating mode on, heating only
            @constraint(m, [t=1:NumTime], H2HP[t] <= M * HP_OPTION[t]) # [kW]          

            # If HP_OPTION = 0, cooling mode on, cooling only
            @constraint(m, [t=1:NumTime], H2C[t] <= M * (1 - HP_OPTION[t])) # [kW]    

            # If TES_Charging_OPTION = 1, charge hot TES only
            @constraint(m, [t=1:NumTime], HP2PCM_H[t] <= M * TES_Charging_OPTION[t]) # [BTU/hr]     
    
            # If TES_Charging_OPTION = 0, charge cold TES only
            @constraint(m, [t=1:NumTime], C2PCM_C[t] <= M * (1 - TES_Charging_OPTION[t])) # [BTU/hr]

            # If TES_Discharging_OPTION = 1, discharge hot TES only
            @constraint(m, [t=1:NumTime], PCM_H2H[t] <= M * TES_Discharging_OPTION[t]) # [BTU/hr]          
    
            # If TES_Discharging_OPTION = 0, discharge cold TES only
            @constraint(m, [t=1:NumTime], PCM_C2H[t] <= M * (1 - TES_Discharging_OPTION[t])) # [BTU/hr]

            # Total TES Charging 
            @expression(m, TES_Charge[t=1:NumTime], HP2PCM_H[t] + C2PCM_C[t]); # [BTU/hr]

            # Total TES Discharging 
            @expression(m, TES_Discharge[t=1:NumTime], PCM_H2H[t] + PCM_C2H[t]); # [BTU/hr]

            # If TES_CD_OPTION = 1, charge only
            @constraint(m, [t=1:NumTime], TES_Charge[t] <= M * TES_CD_OPTION[t]) # [BTU/hr]          

            # If TES_CD_OPTION = 0, discharge only
            @constraint(m, [t=1:NumTime], TES_Discharge[t] <= M * (1 - TES_CD_OPTION[t])) # [BTU/hr] 

            # If B_OPTION = 1, charge only
            @constraint(m, [t=1:NumTime], PV2B[t] <= M * B_OPTION[t]) # [kW]          

            # If B_OPTION = 0, discharge only
            @constraint(m, [t=1:NumTime], B2H[t] <= M * (1 - B_OPTION[t])) # [kW] 
        end
    end

    # display(latex_formulation(m)) # Display the model in latex formulation

    ########### Solve  ##########
    optimize!(m); 
    ########### Model Results  ##########
    # Optimal Actions
    begin
        BestActions = zeros(13);
        BestActions[1] = value.(PV2H[1]);
        BestActions[2] = value.(PV2G[1]);
        BestActions[3] = value.(PV2B[1]);
        BestActions[4] = value.(B2H[1]);
        BestActions[5] = value.(G2H[1]);
        BestActions[6] = value.(H2HP[1]);
        BestActions[7] = value.(HP2H[1]);
        BestActions[8] = value.(H2C[1]);
        BestActions[9] = value.(C2H[1]);
        BestActions[10] = value.(HP2PCM_H[1]);
        BestActions[11] = value.(C2PCM_C[1]);
        BestActions[12] = value.(PCM_H2H[1]);
        BestActions[13] = value.(PCM_C2H[1]);
    end
    # States
    begin
        New_J_States = zeros(4)
        New_J_States[1] = value.(InStorageBattery[2])
        New_J_States[2] = value.(InStoragePCM_H[2])
        New_J_States[3] = value.(InStoragePCM_C[2])
        New_J_States[4] = value.(TemperatureIndoor[2])
    end
    # Cost
    CurrentCost = δt * value.(G2H[1])

    return BestActions, New_J_States, CurrentCost
end

function Optimize_MPCDirect_full(Weather_forecast_with_PV, Schedule, Opt_Horizon, TimeStart, NumRun, current_states)
    """ Performs the MPC algorithm for the full run.
        Inputs: 
        - Weather_forecast_with_PV: list of weather forecasts of length Opt_horizon (typically one week) for the full year. Should be the output of Make_weather_forecast_full_year.
        - Schedule: either SimpleSchedule or ComplexSchedule, over the full year.
        - Opt_Horizon: length of the forecast horizon (typically one week).
        - TimeStart, NumRun: from the larger model.
        - current_states: Initial states of the system.
        Outputs: 
        - Costs (time series for the loss of load), AccumulatedCosts (accumulated loss of load), 
        - AllStates (time series for the consecutive states of the system), ActionMatrix (matrix of operation schedules)
        - TotalCost, runtime
    """
    # Initialization
    begin
        Costs = zeros(NumRun)
        AccumulatedCosts = zeros(NumRun)
        AllStates = zeros(NumRun, 4)
        ActionMatrix = zeros(NumRun, 13)
        TotalCost = 0.0
    end

    # Run the MPC algorithm
    begin
        for i = TimeStart:f_run:NumRun
    
            println("Current time: $i")
            week_weather_forecast = copy(Weather_forecast_with_PV[i])
            println("Extracted forecast")
            
            operation_schedule, new_states, currentcost = Optimize_MPCDirect_onestep(week_weather_forecast, Schedule[i : i+Opt_Horizon-1, :], current_states);
            println("Optimized: $i")
            print(new_states)

            # Update the current states of heat, chill, and battery storage as well as indoor temperature.
            for x = 1:4
                global current_states[x] = new_states[x];
            end
            
            # Add the cost
            Costs[i] = currentcost;
            AccumulatedCosts[i] = sum(Costs[n] for n = 1:i)
            AllStates[i, :] = new_states
            ActionMatrix[i, :] = operation_schedule
            TotalCost += currentcost
            
        end
    end

    return Costs, AccumulatedCosts, AllStates, ActionMatrix, TotalCost

end
    
# Smooth MPC Algorithm

function Optimize_MPCSmooth_onestep(Weather, Schedules, J_States, PreviousActions, HoursPenalized=1, ActionPenalty=1)

    ########## Data Preparations  ##########  
    begin
        NumTime, TemperatureAmbient, TemperatureAmbientC, PercentLighting, PercentPlug, PercentOccupied, PVGeneration, RadHeatGain, RadCooling, CFM, InStorageBattery_1, InStoragePCM_H_1, InStoragePCM_C_1, TemperatureIndoor_1 = DataPreparation(Weather, Schedules, J_States)
    end  
    ########## Declare model  ##########
    begin
        # Define the model name and solver. In this case, model name is "m"
        # m = Model(Clp.Optimizer)
        # m = Model(Ipopt.Optimizer)
        begin
            # Set path to license (for those using Gurobi)
            
            
            m = Model(Gurobi.Optimizer)
            # set_optimizer_attribute(m, "OptimalityTol", 1e-3)   # Optimality tolerance
            set_optimizer_attribute(m, "OutputFlag", 0)         # Suppress solver output
        end
    end
    ######## Decision variables ########
    begin
        @variable(m, PV2H[1:NumTime] >= 0); # [kW] electrical power transfer from PV to home (Berg)
        PV2H_prev = PreviousActions[1, :]

        @variable(m, PV2G[1:NumTime] >= 0); # [kW] electrical power transfer from PV to ground (curtailment)
        PV2G_prev = PreviousActions[2, :]

        @variable(m, PV2B[1:NumTime] >= 0); # [kW] electrical power transfer from PV to battery
        PV2B_prev = PreviousActions[3, :]

        @variable(m, B2H[1:NumTime] >= 0); # [kW] electrical power transfer from battery to home (Berg)
        B2H_prev = PreviousActions[4, :]

        @variable(m, G2H[1:NumTime] >= 0); # [kW] Loss of Load
        G2H_prev = PreviousActions[5, :]

        @variable(m, H2HP[1:NumTime] >= 0); # [kW] electrical power transfer from home (Berg) to heat pump heating mode
        H2HP_prev = PreviousActions[6, :]

        @variable(m, HP2H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from heat pump heating mode to home (Berg)
        HP2H_prev = PreviousActions[7, :]

        @variable(m, H2C[1:NumTime] >= 0); # [kW] electrical power transfer from home (Berg) to heat pump cooling unit
        H2C_prev = PreviousActions[8, :]

        @variable(m, C2H[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from heat pump cooling mode to home (Berg)
        C2H_prev = PreviousActions[9, :]

        @variable(m, HP_OPTION[1:NumTime], Bin) # Heating Mode = 1, Cooling Mode = 0

        @variable(m, TES_Charging_OPTION[1:NumTime], Bin) # Charge hot TES = 1, Charge cold TES = 0

        @variable(m, TES_Discharging_OPTION[1:NumTime], Bin) # Discharge hot TES = 1, Discharge cold TES = 0

        @variable(m, TES_CD_OPTION[1:NumTime], Bin) # Charge TES = 1, Discharge TES = 0

        @variable(m, B_OPTION[1:NumTime], Bin) # Charge battery = 1, Discharge battery = 0

        @variable(m, HP2PCM_H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from heat pump heating mode to PCM heating storage
        HP2PCM_H_prev = PreviousActions[10, :]

        @variable(m, C2PCM_C[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from heat pump cooling mode to PCM cooling storage
        C2PCM_C_prev = PreviousActions[11, :]

        @variable(m, PCM_H2H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from PCM heating storage to home (Berg)
        PCM_H2H_prev = PreviousActions[12, :]

        @variable(m, PCM_C2H[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from PCM cooling storage to home (Berg)
        PCM_C2H_prev = PreviousActions[13, :]

        @variable(m, InStorageBattery[1:NumTime] >= 0); # [kWh] Battery Remaining Charge

        @variable(m, InStoragePCM_H[1:NumTime] >= 0); # [BTU] PCM Heating Remaining Charge

        @variable(m, InStoragePCM_C[1:NumTime] >= 0); # [BTU] PCM Cooling Remaining Charge

        @variable(m, TemperatureIndoor[1:NumTime] >= 0); # [°F] Indoor Air Temperature

    end
    ############ Objective Functions #############
    begin
        # Penalty for loss of load [$]
        @expression(m, penalty_l, δt * sum(G2H[t] for t = 1:NumTime))

        # # Total Cost over Optimization Horizon [$]
        # @objective(m, Min, penalty_l);

        @expression(m, penalty_a, δt * ActionPenalty * sum((PV2H[t-1] - PV2H_prev[t])^2 + (G2H[t-1] - G2H_prev[t])^2
        + (PV2G[t-1] - PV2G_prev[t])^2 + (PV2B[t-1] - PV2B_prev[t])^2 + (B2H[t-1] - B2H_prev[t])^2 
        + (H2HP[t-1] - H2HP_prev[t])^2 + (HP2H[t-1] - HP2H_prev[t])^2 + (H2C[t-1] - H2C_prev[t])^2 
        + (C2H[t-1] - C2H_prev[t])^2 + (HP2PCM_H[t-1] - HP2PCM_H_prev[t])^2 + (C2PCM_C[t-1] - C2PCM_C_prev[t])^2 
        + (PCM_H2H[t-1] - PCM_H2H_prev[t])^2 + (PCM_C2H[t-1] - PCM_C2H_prev[t])^2 for t = 2:HoursPenalized+1)
        )

        # Total Cost over Optimization Horizon [$] with penalty on deviation from previous actions
        @objective(m, Min, penalty_l + penalty_a*ActionPenalty);
    end
    ############# Expressions ############
    begin
        # DeltaTemp 
        @expression(m, TempDelta[t=1:NumTime], TemperatureAmbient[t] - TemperatureIndoor[t]); # [°F]

        # Electricity usage from lighting 
        @expression(m, E_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t]); # [kW]

        # Electricity usage from plugs
        @expression(m, E_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t]); # [kW]

        # Total electricity usage
        @expression(m, E_total[t=1:NumTime], E_Lighting[t] + E_Plugs[t]); # [kW]

        # Calculate Ventilation
        @expression(m, CFMVen[t=1:NumTime], min(Rp * PercentOccupied[t] * MaxOccupancy + Ra * Area, PercentOccupied[t] * MaxOccupancy * Ventilation)); # [ft^3/min]

        # Heat gain from lighting
        @expression(m, Q_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t] * 3412.14); # [BTU/hr]

        # Heat gain from plugs
        @expression(m, Q_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t] * 3412.14); # [BTU/hr]

        # Heat gain from occupancy
        @expression(m, Q_Occupancy[t=1:NumTime], MaxOccupancy * PercentOccupied[t] * TotalPersonHeat); # [BTU/hr]

        # Heat gain from infiltration
        # 1.08 = specific heat capacity of air at STP:0.24 [BTU/(lb*°F)] * air density at STP:0.075 [lb/ft^3] * 60 [min/hr] 
        @expression(m, Q_Infiltration[t=1:NumTime], 1.08 * TempDelta[t] * CFM[t]); # [BTU/hr]

        # Heat gain from ventilation
        @expression(m, Q_Ventilation[t=1:NumTime], 1.08 * TempDelta[t] * CFMVen[t]); # [BTU/hr]

        # Heat gain from lighting, plugs, occupancy, ventilation, infiltration
        @expression(m, Q_Others[t=1:NumTime], Q_Infiltration[t] + Q_Ventilation[t] + Q_Occupancy[t] + Q_Lighting[t] + Q_Plugs[t]); # [BTU/hr]

        # Heat gain through structural evenlope
        @expression(m, Q_Envelope[t=1:NumTime], UA * TempDelta[t]); # [BTU/hr]

        # Heat gain from solar radiation
        @expression(m, Q_Rad[t=1:NumTime], SHGC * RadHeatGain[t]); # [BTU/hr]

        # Radiative cooling
        @expression(m, Q_RadCool[t=1:NumTime], RCC * RadCooling[t]); # [BTU/hr]

        # Detailed COP of heating and cooling (linearized using TemperatureIndoor = 22 [°C])
        # COP (Heating Home)
        @expression(m, COP_Heating_H[t=1:NumTime], HP_a + HP_b * (T_indoor_constant - TemperatureAmbientC[t]) + HP_c * (T_indoor_constant - TemperatureAmbientC[t])^2);
        
        # COP (Heating PCM H)
        @expression(m, COP_Heating_PCM[t=1:NumTime], HP_a + HP_b * (48 - TemperatureAmbientC[t]) + HP_c * (48 -TemperatureAmbientC[t])^2); 

        # COP (Cooling Home)
        @expression(m, COP_Cooling_H[t=1:NumTime], HP_a + HP_b * (TemperatureAmbientC[t] - T_indoor_constant) + HP_c * (TemperatureAmbientC[t] - T_indoor_constant)^2);
        
        # COP (Cooling PCM C)
        @expression(m, COP_Cooling_PCM[t=1:NumTime], HP_a + HP_b * (TemperatureAmbientC[t] - 11) + HP_c * (TemperatureAmbientC[t] - 11)^2); 
        
    end
    ############# Constraints ############
    begin
        # Initialization Constraints with States from last timestep
        begin
            # Set point temperature range constraints
            @constraint(m, TemperatureIndoor[1] == TemperatureIndoor_1); # [°F]

            # PCM heating storage initialization constraint, node at PCM heating storage
            @constraint(m, InStoragePCM_H[1] == InStoragePCM_H_1); # [BTU]

            # PCM cooling storage initialization constraint, node at PCM cooling storage
            @constraint(m, InStoragePCM_C[1] == InStoragePCM_C_1); # [BTU]

            # Battery storage initialization constraint
            @constraint(m, InStorageBattery[1] == InStorageBattery_1); # [kWh] 
        end

        # Set point temperature range constraints 
        @constraint(m, [t=1:NumTime], TemperatureIndoor[t] <= SetPointT_High); # [°F]

        @constraint(m, [t=1:NumTime], TemperatureIndoor[t] >= SetPointT_Low); # [°F]

        # Internal temperature balance evolution constraint
        @constraint(m, [t=1:NumTime-1], TemperatureIndoor[t+1] == TemperatureIndoor[t] + 
        (δt/TC)*(Q_Others[t] + Q_Envelope[t] + Q_Rad[t] - Q_RadCool[t] + HP2H[t] - C2H[t] + PCM_H2H[t] - PCM_C2H[t])); # [°F]
        
        # PV energy balance constraint, node at PV
        @constraint(m, [t=1:NumTime], PVGeneration[t] * PVSize ==  PV2B[t] + PV2H[t] + PV2G[t]); # [kW]
    
        # House electricity load constraint, node at house, battery efficiency modeled
        @constraint(m, [t=1:NumTime], E_total[t] + H2HP[t] + H2C[t] == PV2H[t] * η_PVIV + B2H[t] * η + G2H[t]); # [kW]

        # Battery storage balance constraint, node at battery, battery leakage modeled, battery efficiency modeled
        @constraint(m, [t=1:NumTime-1], InStorageBattery[t+1] == InStorageBattery[t] * δt * (1 - BatteryLoss) +  δt * (PV2B[t] * η - B2H[t])); # [kWh]

        # Battery discharging constraint, node at battery 
        @constraint(m, [t=1:NumTime], δt * B2H[t] <= InStorageBattery[t]); # [kWh]
        
        # Battery power inverter constraint, node at battery (inverter power constraint) 
        @constraint(m, [t=1:NumTime], B2H[t] + PV2B[t] <= InverterSize); # [kW]
        
        # Battery storage size constraint, node at battery
        @constraint(m, [t=1:NumTime], InStorageBattery[t] <= BatterySize); # [kWh]
        
        # Battery storage max discharge constraint, node at battery, always at least 20% full 
        @constraint(m, [t=1:NumTime], InStorageBattery[t] >= BatterySize * (1-MaxDischarge)); # [kWh]
        
        # Heating and Cooling Constraints (Needs Remodeling)
        begin
            # Detailed COP of heating and cooling (requires nonlinear optimization solver)
            # Heating (from kW to BTU/hr), node at heat pump heating option
            @constraint(m, [t=1:NumTime], H2HP[t] * 3412.14 == HP2PCM_H[t]/COP_Heating_PCM[t] + HP2H[t]/COP_Heating_H[t]);
            
            # Cooling (from kWh to BTU), node at heat pump cooling option
            @constraint(m, [t=1:NumTime], H2C[t] * 3412.14 == C2PCM_C[t]/COP_Cooling_PCM[t] + C2H[t]/COP_Cooling_H[t]);

            #=
            # Heating (from kW to BTU/hr), node at heat pump heating option
            @constraint(m, [t=1:NumTime], H2HP[t] * 3412.14 == (HP2PCM_H[t] + HP2H[t])/COP_H); # [BTU/hr]

            # Cooling (from kW to BTU/hr), node at heat pump cooling option
            @constraint(m, [t=1:NumTime], H2C[t] * 3412.14 == (C2PCM_C[t] + C2H[t])/COP_C); # [BTU/hr]
            =#

            # Heat pump power capacity constraint, node at heat pump heating mode
            @constraint(m, [t=1:NumTime], H2HP[t] <= Cap_HP_H); # [kW]

            # Heat pump power capacity constraint, node at heat pump cooling mode
            @constraint(m, [t=1:NumTime], H2C[t] <= Cap_HP_C); # [kW]

            # PCM heating storage balance constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime-1], InStoragePCM_H[t+1] == InStoragePCM_H[t] + δt * (HP2PCM_H[t] - PCM_H2H[t])); # [BTU]
            
            # PCM cooling storage balance constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime-1], InStoragePCM_C[t+1] == InStoragePCM_C[t] + δt * (C2PCM_C[t] - PCM_C2H[t])); # [BTU]
            
            # PCM heating storage discharging constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime], δt * PCM_H2H[t] <= InStoragePCM_H[t]); # [BTU]
            
            # PCM cooling storage discharging constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime], δt * PCM_C2H[t] <= InStoragePCM_C[t]); # [BTU]
            
            # PCM heating storage size constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime], InStoragePCM_H[t] <= PCM_H_Size); # [BTU]

            # PCM cooling storage size constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime], InStoragePCM_C[t] <= PCM_C_Size); # [BTU]
        end
        # Single Operation Choice Constraints
        begin
            # If HP_OPTION = 1, heating mode on, heating only
            @constraint(m, [t=1:NumTime], H2HP[t] <= M * HP_OPTION[t]) # [kW]          

            # If HP_OPTION = 0, cooling mode on, cooling only
            @constraint(m, [t=1:NumTime], H2C[t] <= M * (1 - HP_OPTION[t])) # [kW]    

            # If TES_Charging_OPTION = 1, charge hot TES only
            @constraint(m, [t=1:NumTime], HP2PCM_H[t] <= M * TES_Charging_OPTION[t]) # [BTU/hr]     
    
            # If TES_Charging_OPTION = 0, charge cold TES only
            @constraint(m, [t=1:NumTime], C2PCM_C[t] <= M * (1 - TES_Charging_OPTION[t])) # [BTU/hr]

            # If TES_Discharging_OPTION = 1, discharge hot TES only
            @constraint(m, [t=1:NumTime], PCM_H2H[t] <= M * TES_Discharging_OPTION[t]) # [BTU/hr]          
    
            # If TES_Discharging_OPTION = 0, discharge cold TES only
            @constraint(m, [t=1:NumTime], PCM_C2H[t] <= M * (1 - TES_Discharging_OPTION[t])) # [BTU/hr]

            # Total TES Charging 
            @expression(m, TES_Charge[t=1:NumTime], HP2PCM_H[t] + C2PCM_C[t]); # [BTU/hr]

            # Total TES Discharging 
            @expression(m, TES_Discharge[t=1:NumTime], PCM_H2H[t] + PCM_C2H[t]); # [BTU/hr]

            # If TES_CD_OPTION = 1, charge only
            @constraint(m, [t=1:NumTime], TES_Charge[t] <= M * TES_CD_OPTION[t]) # [BTU/hr]          

            # If TES_CD_OPTION = 0, discharge only
            @constraint(m, [t=1:NumTime], TES_Discharge[t] <= M * (1 - TES_CD_OPTION[t])) # [BTU/hr] 

            # If B_OPTION = 1, charge only
            @constraint(m, [t=1:NumTime], PV2B[t] <= M * B_OPTION[t]) # [kW]          

            # If B_OPTION = 0, discharge only
            @constraint(m, [t=1:NumTime], B2H[t] <= M * (1 - B_OPTION[t])) # [kW] 
        end
    end

    # display(latex_formulation(m)) # Display the model in latex formulation

    ########### Solve  ##########
    optimize!(m); 
    ########### Model Results  ##########
    # Optimal Actions
    begin
        BestActions = zeros(13, HoursPenalized+1);
        BestActions[1, :] = value.(PV2H[1:HoursPenalized+1]);
        BestActions[2, :] = value.(PV2G[1:HoursPenalized+1]);
        BestActions[3, :] = value.(PV2B[1:HoursPenalized+1]);
        BestActions[4, :] = value.(B2H[1:HoursPenalized+1]);
        BestActions[5, :] = value.(G2H[1:HoursPenalized+1]);
        BestActions[6, :] = value.(H2HP[1:HoursPenalized+1]);
        BestActions[7, :] = value.(HP2H[1:HoursPenalized+1]);
        BestActions[8, :] = value.(H2C[1:HoursPenalized+1]);
        BestActions[9, :] = value.(C2H[1:HoursPenalized+1]);
        BestActions[10,:] = value.(HP2PCM_H[1:HoursPenalized+1]);
        BestActions[11,:] = value.(C2PCM_C[1:HoursPenalized+1]);
        BestActions[12,:] = value.(PCM_H2H[1:HoursPenalized+1]);
        BestActions[13,:] = value.(PCM_C2H[1:HoursPenalized+1]);
    end
    # States
    begin
        New_J_States = zeros(4)
        New_J_States[1] = value.(InStorageBattery[2])
        New_J_States[2] = value.(InStoragePCM_H[2])
        New_J_States[3] = value.(InStoragePCM_C[2])
        New_J_States[4] = value.(TemperatureIndoor[2])
    end
    # Cost
    CurrentCost = δt * value.(G2H[1])

    return BestActions, New_J_States, CurrentCost
end

function Optimize_MPCSmooth_full(Weather_forecast_with_PV, Schedule, Opt_Horizon, TimeStart, NumRun, current_states, hours_penalized=1, deviation_penalty=1)
    """ Performs the MPC algorithm for the full run.
        Inputs: 
        - Weather_forecast_with_PV: list of weather forecasts of length Opt_horizon (typically one week) for the full year. Should be the output of Make_weather_forecast_full_year.
        - Schedule: either SimpleSchedule or ComplexSchedule, over the full year.
        - Opt_Horizon: length of the forecast horizon (typically one week).
        - TimeStart, NumRun: from the larger model.
        - current_states: Initial states of the system.
        Outputs: 
        - Costs (time series for the loss of load), AccumulatedCosts (accumulated loss of load), 
        - AllStates (time series for the consecutive states of the system), ActionMatrix (matrix of operation schedules)
        - TotalCost, runtime
    """

    # Initialization
    begin
        Costs = zeros(NumRun)
        AccumulatedCosts = zeros(NumRun)
        AllStates = zeros(NumRun, 4)
        ActionMatrix = zeros(13, NumRun)
        TotalCost = 0.0
    end

    # Run the MPC algorithm
    begin
        operational_schedule = zeros(13, 2) 
        for i = TimeStart:f_run:NumRun
            local previous_actions = operational_schedule

            println("Current time: $i")
            week_weather_forecast = copy(Weather_forecast_with_PV[i])
            println("Extracted forecast")
            
            operation_schedule, new_states, currentcost = Optimize_MPCSmooth_onestep(week_weather_forecast, Schedule[i : i+Opt_Horizon-1, :], current_states, previous_actions, hours_penalized, deviation_penalty);
            println("Optimized: $i")
            print(new_states)

            # Update the current states of heat, chill, and battery storage as well as indoor temperature.
            for x = 1:4
                global current_states[x] = new_states[x];
            end
            
            # Add the cost
            Costs[i] = currentcost;
            AccumulatedCosts[i] = sum(Costs[n] for n = 1:i)
            AllStates[i, :] = new_states
            ActionMatrix[:, i] = operation_schedule[:, 1]
            TotalCost += currentcost
            previous_actions = operation_schedule
            
        end
    end

    return Costs, AccumulatedCosts, AllStates, ActionMatrix, TotalCost
end
