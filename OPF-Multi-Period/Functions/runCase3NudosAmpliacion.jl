# --------------------------------
# Functions/runCase3NudosAmpliacion.jl
# Case: 3 Nudos Ampliación (solar, gas, bio, wind) with UC
# --------------------------------
using CSV
using DataFrames
using Dates
using JuMP: value          

include(joinpath(@__DIR__, "..", "OPF", "DC_OPF", "DC_OPF_UC.jl"))
using .DC_OPF_UC_Module: DC_OPF_UC


function run_case_3nudos_ampliacion(; solver_name::String)

    bMVA  = 100
    hours = 24
    nN    = 3
    nL    = 3
    s     = solver_name

    case_path = joinpath(@__DIR__, "..", "Cases", "3 Nudos Ampliación")

    dGen   = CSV.read(joinpath(case_path, "generatorData.csv"), DataFrame)

    println("\n------ dGen read from CSV ------")
    println(dGen)
    println()

    println("Cm       = ", dGen.Cm_EUR_MWh)
    println("CnL      = ", dGen.CnL_EUR_h)
    println("C_start  = ", dGen.C_start_EUR)
    println("Pmin_MW  = ", dGen.Pmin_MW)
    println("Pmax_MW  = ", dGen.Pmax_MW)
    println("RR_MW/h  = ", dGen.RR_MW_per_h)
    println("Tmut_h   = ", dGen.Tmut_h)
    println("Tmdt_h   = ", dGen.Tmdt_h)
    println("status   = ", dGen.status)
    println("bus      = ", dGen.bus)
    println("---------------------------------\n")

    dDem   = CSV.read(joinpath(case_path, "nodeData.csv"),      DataFrame)
    dSolar = CSV.read(joinpath(case_path, "solarData.csv"),     DataFrame)
    dWind  = CSV.read(joinpath(case_path, "windData.csv"),      DataFrame)
    dLine  = CSV.read(joinpath(case_path, "LineData.csv"),      DataFrame)

    # Build Node :: Vector{DataFrame}
    Node = Vector{DataFrame}(undef, hours)
    for h in 1:hours
        df_h = filter(:hour => ==(h), dDem)
        Node[h] = select(df_h, :bus_i, :type, :Pd)
    end

    # DC-OPF with UC
    m, solGen, solFlows, solVoltage, solCosts, solCurt =
        DC_OPF_UC(dLine, dGen, Node, nN, nL, bMVA, s, hours, dSolar, dWind)

    # 
    y_st = m[:y_st]   # órdenes de arranque
    y_sg = m[:y_sg]   # start generating

    # Base
    base = DataFrame(
        hour = repeat(1:hours, inner = nN),
        bus  = repeat(1:nN, hours),)

    base.P_gen_MW        = zeros(Float64, nrow(base))   
    base.Pd_MW           = zeros(Float64, nrow(base))
    base.P_curt_RES_MW   = zeros(Float64, nrow(base))
    base.P_solar_av_MW   = zeros(Float64, nrow(base))
    base.P_wind_av_MW    = zeros(Float64, nrow(base))
    base.P_solar_use_MW  = zeros(Float64, nrow(base))
    base.P_wind_use_MW   = zeros(Float64, nrow(base))

    # Solar long
    solar_long = DataFrame(hour = Int[], bus = Int[], P_solar_av_MW = Float64[])
    for h in 1:hours
        col = Symbol("h$h")
        for r in 1:nrow(dSolar)
            push!(solar_long, (
                hour = h,
                bus  = dSolar.bus[r],
                P_solar_av_MW = dSolar[r, col],
            ))
        end
    end

    # Wind long
    wind_long = DataFrame(hour = Int[], bus = Int[], P_wind_av_MW = Float64[])
    for h in 1:hours
        col = Symbol("h$h")
        for r in 1:nrow(dWind)
            push!(wind_long, (
                hour = h,
                bus  = dWind.bus[r],
                P_wind_av_MW = dWind[r, col],
            ))
        end
    end

    # base
    for row in 1:nrow(base)
        h = base.hour[row]
        b = base.bus[row]

        # Total thermal generation in the bus
        idxs_g = findall((solGen.hour .== h) .& (solGen.bus .== b))
        if !isempty(idxs_g)
            base.P_gen_MW[row] = sum(solGen.powerGen[idxs_g])
        end

        # Demand
        idx_d = findfirst((dDem.hour .== h) .& (dDem.bus_i .== b))
        if idx_d !== nothing
            base.Pd_MW[row] = dDem.Pd[idx_d]
        end

        # Curtailment (solar in bus 1 and wind in bus 3)
        idx_c = findfirst((solCurt.hour .== h) .& (solCurt.bus .== b))
        if idx_c !== nothing
            base.P_curt_RES_MW[row] = solCurt.P_curtailment[idx_c]
        end

        # Solar available
        idx_s = findfirst((solar_long.hour .== h) .& (solar_long.bus .== b))
        if idx_s !== nothing
            base.P_solar_av_MW[row] = solar_long.P_solar_av_MW[idx_s]
        end

        # Viento available
        idx_w = findfirst((wind_long.hour .== h) .& (wind_long.bus .== b))
        if idx_w !== nothing
            base.P_wind_av_MW[row] = wind_long.P_wind_av_MW[idx_w]
        end

        # Solar used (bus 1)
        if b == 1
            base.P_solar_use_MW[row] =
                base.P_solar_av_MW[row] - base.P_curt_RES_MW[row]
        end

        # Wind used (bus 3)
        if b == 3
            base.P_wind_use_MW[row] =
                base.P_wind_av_MW[row] - base.P_curt_RES_MW[row]
        end
    end

    sort!(base, [:hour, :bus])
    total_cost = solCosts[solCosts.hour .== -1, :operation_cost][1]

    # ---------- Results tables per BUS ----------

    println("\n============ RESULTS 3 NUDOS AMPLIACIÓN (UC) ============\n")

    # BUS 1: Solar + Demand
    println(">>> BUS 1 (Solar + Demand)")
    bus1_df = filter(:bus => x -> x == 1, base)
    bus1_df = select(bus1_df,
        :hour,
        :P_solar_av_MW,
        :P_solar_use_MW,
        :Pd_MW,
        :P_curt_RES_MW,)
    show(bus1_df, allrows = true, allcols = true)
    println("\n------------------------------------------------------------\n")

    # BUS 2: Gas + Biomass + Demand
    println(">>> BUS 2 (Generación térmica Gas + Biomasa + Demanda)")
    
    bus2_df = DataFrame(
        hour          = 1:hours,
        P_gen_gas_MW  = zeros(Float64, hours),
        P_gen_bio_MW  = zeros(Float64, hours),
        P_gen_total_MW = zeros(Float64, hours),
        Pd_MW         = zeros(Float64, hours),
        y_st_gas      = zeros(Int, hours),
        y_sg_gas      = zeros(Int, hours),
        y_st_bio      = zeros(Int, hours),
        y_sg_bio      = zeros(Int, hours),)

    for h in 1:hours

        # --- Demand BUS 2 ---
        idx_d = findfirst((dDem.hour .== h) .& (dDem.bus_i .== 2))
        if idx_d !== nothing
            bus2_df.Pd_MW[h] = dDem.Pd[idx_d]
        end

        # --- GAS (gen 1) ---
        idx_gas = findfirst((solGen.hour .== h) .& (solGen.gen .== 1))
        if idx_gas !== nothing
            bus2_df.P_gen_gas_MW[h] = solGen.powerGen[idx_gas]
        end

        # --- BIOMASS (gen 2) ---
        idx_bio = findfirst((solGen.hour .== h) .& (solGen.gen .== 2))
        if idx_bio !== nothing
            bus2_df.P_gen_bio_MW[h] = solGen.powerGen[idx_bio]
        end

        # Total
        bus2_df.P_gen_total_MW[h] =
            bus2_df.P_gen_gas_MW[h] + bus2_df.P_gen_bio_MW[h]

        # Gas = generator 1
        bus2_df.y_st_gas[h] = Int(round(value(y_st[1, h])))
        bus2_df.y_sg_gas[h] = Int(round(value(y_sg[1, h])))

        # Biomasa = generator 2
        bus2_df.y_st_bio[h] = Int(round(value(y_st[2, h])))
        bus2_df.y_sg_bio[h] = Int(round(value(y_sg[2, h])))
    end

    show(bus2_df, allrows = true, allcols = true)
    println("\n------------------------------------------------------------\n")

    # BUS 3: Wind + Demand
    println(">>> BUS 3 (Wind + Demand)")
    bus3_df = filter(:bus => x -> x == 3, base)
    bus3_df = select(bus3_df,
        :hour,
        :P_wind_av_MW,
        :P_wind_use_MW,
        :Pd_MW,
        :P_curt_RES_MW,
    )
    show(bus3_df, allrows = true, allcols = true)
    println("\n------------------------------------------------------------\n")

    println("\nTotal system cost: $(round(total_cost, digits = 3)) €")
    println("============================================================\n")

    # ---------- Save CSV in Results ----------
    results_folder = joinpath(@__DIR__, "..", "Results")
    if !isdir(results_folder)
        mkdir(results_folder)
    end

    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")

    CSV.write(joinpath(results_folder, "3NudosAmpliacion_bus1_$timestamp.csv"), bus1_df)
    CSV.write(joinpath(results_folder, "3NudosAmpliacion_bus2_$timestamp.csv"), bus2_df)
    CSV.write(joinpath(results_folder, "3NudosAmpliacion_bus3_$timestamp.csv"), bus3_df)

    println("Files saved in $results_folder")
    println(" - 3NudosAmpliacion_bus1_$timestamp.csv")
    println(" - 3NudosAmpliacion_bus2_$timestamp.csv")
    println(" - 3NudosAmpliacion_bus3_$timestamp.csv\n")

    return m, base, solGen, solFlows, solVoltage, solCosts, solCurt
end
