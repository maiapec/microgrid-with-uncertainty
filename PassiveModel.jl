# This code is the property of Fred Fan, Civil and Environmental Engineering Department, Stanford University.

# Create a datetime variable for each time instance
function create_datetime_from_row(row)
    # Corrected access to DataFrame row fields
    return DateTime(row.Year, row.Month, row.Day, row.Hour, row.Minute)
end

# Define a moving average function
function Moving_average(x, w)
    return [mean(x[i-w+1:i]) for i in w:length(x)]
end

# Define a function to calculate solar time for each time instance
function add_solar_time(df)
    solartimes = zeros(size(df)[1])
    for i = 1:size(df)[1]
        local_time = df[i, :datetime]
        # Convert local time to Local Solar Time
        LT = hour(local_time) + minute(local_time)/60 + second(local_time)/3600

        # Day of the year
        N = dayofyear(local_time)

        # Calculate the Equation of Time (Eqt)
        if 1 ≤ N ≤ 106
            Eqt = -14.2 * sin(π * (N+7)/111)
        elseif 107 ≤ N ≤ 166
            Eqt = 4 * sin(π * (N-106)/59)
        elseif 167 ≤ N ≤ 246
            Eqt = 6.5 * sin(π * (N-166)/80)
        elseif 247 ≤ N ≤ 366
            Eqt = 16.4 * sin(π * (N-247)/113)  
        end

        # Calculate the Solar Time (T_solar)
        T_solar = LT + Eqt/60 + (4*(standard_meridian_longitude - longitude))/60 # [DEGREES]
        solartimes[i] = T_solar
    end
    df[!, :SolarTime] = solartimes
    return df 
end 

# Define a function to calculate solar irradiance on a certain surface for each time instance
function calculate_irradiance(local_time, DNI, DHI, GHI, surface_azimuth, tilt_angle)
    # Convert local time to Local Solar Time
    LT = hour(local_time) + minute(local_time)/60 + second(local_time)/3600

    # Day of the year
    N = dayofyear(local_time)

    # Calculate the Equation of Time (Eqt)
    if 1 ≤ N ≤ 106
        Eqt = -14.2 * sin(π * (N+7)/111)
    elseif 107 ≤ N ≤ 166
        Eqt = 4 * sin(π * (N-106)/59)
    elseif 167 ≤ N ≤ 246
        Eqt = 6.5 * sin(π * (N-166)/80)
    elseif 247 ≤ N ≤ 366
        Eqt = 16.4 * sin(π * (N-247)/113)  
    end

    # Calculate the Solar Time (T_solar)
    T_solar = LT + Eqt/60 + (4*(standard_meridian_longitude - longitude))/60 # [DEGREES]

    # Solar Declination (δ)
    δ = 23.45 * sin(deg2rad((360/365) * (284 + N))) # [DEGREES]

    # Solar Time Angle (ω)
    ω = 15 * (T_solar - 12) # [DEGREES]

    # Latitude (φ)
    φ = deg2rad(latitude) # [RADIANS]

    # Solar Elevation Angle (α)
    α = asin(sin(φ) * sin(deg2rad(δ)) + cos(φ) * cos(deg2rad(δ)) * cos(deg2rad(ω))) # [RADIANS]
    
    # Solar Azimuth Angle (Ψ)
    if ω < 0
        Ψ = acos((sin(deg2rad(δ)) - sin(α) * sin(φ)) / (cos(α) * cos(φ))) # [RADIANS]
    else
        Ψ = 2 * π - acos((sin(deg2rad(δ)) - sin(α) * sin(φ)) / (cos(α) * cos(φ))) # [RADIANS]
    end

    # Surface tilt β
    β = deg2rad(tilt_angle) # [RADIANS]
    # Surface azimuth γ
    γ = deg2rad(surface_azimuth) # [RADIANS]
    # Angle of Incidence (θ)
    θ = acos(sin(α) * cos(β) + cos(α) * sin(β) * cos(γ - Ψ)) # [RADIANS]

    # Beam Irradiance I_b
    # no direct beam irradiance if the angle of incidence is negative
    if cos(θ) < 0
        I_b = 0
    else
        I_b = DNI * cos(θ) # [W/m^2]
    end
    
    # Diffuse irradiance I_d
    I_d = max(DHI * (1+cos(β))/2 + 0.5*GHI*(0.012*(π/2 - α)-0.04)*(1-cos(β)),0) # [W/m^2]
    
    # Reflected irradiance I_r
    I_r = max(GHI * albedo * (1 - cos(β)) / 2, 0) # [W/m^2]
    
    # Total irradiance I_total
    I_total = I_d + I_r + I_b # [W/m^2]
    return I_total
end

# Define a function to calculate total solar irradiance received for each time instance
function calculate_solarheatgain(datetime_obj, DNI, DHI, GHI)
    
    I_t = calculate_irradiance(datetime_obj, DNI, DHI, GHI, 0, 0) # [W/m^2] top side

    I_e = calculate_irradiance(datetime_obj, DNI, DHI, GHI, 90-Berg_tilt, 90) # [W/m^2] east side

    I_s = calculate_irradiance(datetime_obj, DNI, DHI, GHI, 180-Berg_tilt, 90) # [W/m^2] south side

    I_w = calculate_irradiance(datetime_obj, DNI, DHI, GHI, 270-Berg_tilt, 90) # [W/m^2] west side

    I_n = calculate_irradiance(datetime_obj, DNI, DHI, GHI, 0-Berg_tilt, 90) # [W/m^2] north side

    # Calculate total solar heat gain from all 5 sides (3.412[BTU/hr/W], 10.764[ft^2/m^2])
    TotalSolarHeatGain = 3.412*(I_t * L_wall * L_wall + (I_e + I_s + I_w + I_n) * L_wall * H_wall)/10.764 # [BTU/hr]

    return TotalSolarHeatGain
end

# Define a function to calculate radiative cooling for each time instance
function Q_radiativeCooling(Ti, Ta, DewPoint, solartime, day_status)
    
    # Calculate dew point temperature in Celcius
    temp_dewpoint = DewPoint # [°C]

    # Emissivity of sky based on dew point temperature in Celcius
    e_sky = 0.741 + 0.0062 * temp_dewpoint # [1] 
    
    # Total exposed area subject to radiative cooling
    exposed_area = (L_wall*L_wall + 4*L_wall*H_wall)/10.764 # [m^2] 
    
    # Radiative Cooling using Stefan–Boltzmann Law
    sigma = 5.67*10^-8 # [W/(m^2*K^4)] Stefan–Boltzmann constant

    # Total exposed area subject to radiative cooling by side
    exposed_area = [L_wall*H_wall, L_wall*H_wall, L_wall*H_wall, L_wall*H_wall, L_wall*L_wall]/10.764 # [m^2] 
    total_exposed_area = sum(exposed_area) # [m^2] 
    
    # Surface temperature by side 
    Ts_list = estimate_Ts(Ti, Ta, solartime, day_status) # [°F]
    Ts = (5/9)*(Ts_list .- 32).+273.15 # [°K]

    T_sky = (5/9)*(Ta - 32)
    # Calculate the net heat loss due to radiation which is the difference in radiation from Berg to ambient and radiation from ambient to Berg.
    Q_radcool = 3.412 * sigma * (e_Berg * sum((Ts.^4).*exposed_area) - e_sky * total_exposed_area *(T_sky+273.15)^4) #[BTU/hr]

    return Q_radcool
end   

# Define a function to estimate structural surface temperatures for each time instance
function estimate_Ts(Ti, Ta, solartime, day_status)
    # If it is at night, the surface temperature is the average between ambient and indoor temperature.
    if day_status == 0
        Ts = (Ta + Ti)/2
        return [Ts, Ts, Ts, Ts, Ts]
    end

    # East side, peaks at 10 am
    Ts_E = max(Ti, Ti * (-10/7) * cos(2*pi*(solartime+2)/24))
    
    # South side, peaks at 2 pm
    Ts_S = max(Ti, Ti * (-10/7) * cos(2*pi*(solartime-2)/24))
    
    # West side, peaks at 4:30 pm
    Ts_W = max(Ti, Ti * (-10/7) * cos(2*pi*(solartime-4.5)/24))
    
    # North side, relatively cool
    Ts_N = Ti * (-6/7) * cos(2*pi*(solartime-48)/72)

    # Top side, not as hot, peaks at 2 pm
    Ts_R = max(Ti, Ti * (-8/7) * cos(2*pi*(solartime-50)/72))

    return [Ts_E, Ts_S, Ts_W, Ts_N, Ts_R]
end

# Define a function to check if the current time instance is daytime (by checking if sun is up) for each time instance
function IsDay(solar_time, local_time)
    # Calculate the Solar Time (T_solar)
    T_solar = solar_time # [DEGREES]
    N = dayofyear(local_time)
    # Solar Declination (δ)
    δ = 23.45 * sin(deg2rad((360/365) * (284 + N))) # [DEGREES]

    # Solar Time Angle (ω)
    ω = 15 * (T_solar - 12) # [DEGREES]

    # Latitude (φ)
    φ = deg2rad(latitude) # [RADIANS]

    # Solar Elevation Angle (α)
    α = asin(sin(φ) * sin(deg2rad(δ)) + cos(φ) * cos(deg2rad(δ)) * cos(deg2rad(ω))) # [RADIANS]

    day_status = 0
    
    if α>0
        day_status = 1
    end
    
    return day_status
end    

# Define a function to calculate dew point temperature for each time instance (currently not used, since NSRDB has dew point data)
function calculate_dew_point(temp_amb, RH)
    # Convert ambient temperature from Fahrenheit to Celcius
    T = (5/9) * (temp_amb - 32) # [°C]

    # Calculate the saturation vapor pressure (es) at the temperature T
    es = 6.112 * exp((17.67 * T) / (T + 243.5))
    
    # Calculate the actual vapor pressure (e)
    e = (RH / 100) * es
    
    # Calculate the dew point temperature (Td)
    Td = (243.5 * log(e / 6.112)) / (17.67 - log(e / 6.112)) # [°C]
    
    return Td
end

# Define a function to calculate PV generation (capacity factor) for each time instance 
function calculate_PV(local_time, DNI, DHI, GHI, surface_azimuth, tilt_angle, T_cell)
    # Convert local time to Local Solar Time
    LT = hour(local_time) + minute(local_time)/60 + second(local_time)/3600

    # Day of the year (N)
    N = dayofyear(local_time)

    # Calculate the Equation of Time (Eqt)
    if 1 ≤ N ≤ 106
        Eqt = -14.2 * sin(π * (N+7)/111)
    elseif 107 ≤ N ≤ 166
        Eqt = 4 * sin(π * (N-106)/59)
    elseif 167 ≤ N ≤ 246
        Eqt = 6.5 * sin(π * (N-166)/80)
    elseif 247 ≤ N ≤ 366 #leap year
        Eqt = 16.4 * sin(π * (N-247)/113)  
    end

    # Calculate the Solar Time (T_solar)
    T_solar = LT + Eqt/60 + (4*(standard_meridian_longitude - longitude))/60 # [DEGREES]

    # Solar Declination (δ)
    δ = 23.45 * sin(deg2rad((360/365) * (284 + N))) # [DEGREES]

    # Solar Time Angle (ω)
    ω = 15 * (T_solar - 12) # [DEGREES]

    # Latitude (φ)
    φ = deg2rad(latitude) # [RADIANS]

    # Solar Elevation Angle (α)
    α = asin(sin(φ) * sin(deg2rad(δ)) + cos(φ) * cos(deg2rad(δ)) * cos(deg2rad(ω))) # [RADIANS]
    
    # Solar Azimuth Angle (Ψ)
    if ω < 0
        Ψ = acos((sin(deg2rad(δ)) - sin(α) * sin(φ)) / (cos(α) * cos(φ))) # [RADIANS]
    else
        Ψ = 2 * π - acos((sin(deg2rad(δ)) - sin(α) * sin(φ)) / (cos(α) * cos(φ))) # [RADIANS]
    end

    # Surface tilt (β)
    β = deg2rad(tilt_angle) # [RADIANS]
    # Surface azimuth (γ)
    γ = deg2rad(surface_azimuth) # [RADIANS]
    # Angle of Incidence (θ)
    θ = acos(sin(α) * cos(β) + cos(α) * sin(β) * cos(γ - Ψ)) # [RADIANS]

    # Beam Irradiance (I_b)
    # no direct beam irradiance if the angle of incidence is negative
    if cos(θ) < 0
        I_b = 0
    else
        I_b = DNI * cos(θ) # [W/m^2]
    end
    
    # Diffuse irradiance (I_d)
    I_d = max(DHI * (1+cos(β))/2 + 0.5*GHI*(0.012*(π/2 - α)-0.04)*(1-cos(β)),0) # [W/m^2]
    
    # Reflected irradiance (I_r)
    I_r = max(GHI * albedo * (1 - cos(β)) / 2, 0) # [W/m^2]
    
    # Total irradiance (I_poa)
    I_poa = I_d + I_r + I_b # [W/m^2]

    # Calculate transmittence 
    begin
        # Calculate angle of refraction into AR coating (θ_2) 
        θ_2 = asin((n_air / n_AR) * sin(θ)) # [RADIANS]
        
        # Calculate transmittance through AR coating (τ_AR) 
        τ_AR = 1 - 0.5 * ((sin(θ_2 - θ)^2) / (sin(θ_2 + θ)^2)) + ((tan(θ_2 - θ)^2) / ((tan(θ_2 + θ)^2))) # [1]

        # Calculate angle of refraction into glass cover (θ_3) 
        θ_3 = asin((n_AR / n_glass) * sin(θ_2)) # [RADIANS]

        # Calculate transmittance through glass cover (τ_glass) 
        τ_glass = 1 - 0.5 * ((sin(θ_3 - θ_2)^2) / (sin(θ_3 + θ_2)^2)) + ((tan(θ_3 - θ_2)^2) / ((tan(θ_3 + θ_2)^2))) # [1]
        
        # Calculate effective transmittance through AR coated module (τ_cover)
        τ_cover = τ_AR * τ_glass # [1]
    end
    
    # Calculate transmitted POA irradiance [I_tr]
    I_tr = I_poa * τ_cover # [W/m^2]

    # Calculate PV DC power output [P_dc]
    P_dc = (I_tr/1000) * P_dc0 * (1 + Γ_t * (T_cell - T_ref)) # [kW/kW]

    # Calculate PV DC power output after system loss [P_PV]
    P_PV = P_dc * η_PV # [kW/kW]

    return P_PV 
end 