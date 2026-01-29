# --------------------------------
# Functions/runCase3NudosUC.jl
# --------------------------------
using CSV
using DataFrames
using Dates

include(joinpath(@__DIR__, "..", "OPF", "DC_OPF", "DC_OPF_UC.jl"))
using .DC_OPF_UC_Module: DC_OPF_UC

function run_case_3nudos_uc(; solver_name::String)

    bMVA  = 100
    hours = 24
    nN    = 3
    nL    = 2
    s     = solver_name

    case_path = joinpath(@__DIR__, "..", "Cases", "3 Nudos UC")

    dGen   = CSV.read(joinpath(case_path, "generatorData.csv"), DataFrame)
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

    # Run UC
    m, solGen, solFlows, solVoltage, solCosts, solCurt =
        DC_OPF_UC(dLine, dGen, Node, nN, nL, bMVA, s, hours, dSolar, dWind)

    # ==========================
    # Build base table
    # ==========================

    base = DataFrame(
        hour = repeat(1:hours, inner = nN),
        bus  = repeat(1:nN, hours),)

    base.P_gen_MW        = zeros(Float64, nrow(base))
    base.Pd_MW           = zeros(Float64, nrow(base))
    base.P_curt_solar_MW = zeros(Float64, nrow(base))
    base.P_solar_av_MW   = zeros(Float64, nrow(base))
    base.P_solar_use_MW  = zeros(Float64, nrow(base))

    # solar long
    solar_long = DataFrame(hour = Int[], bus = Int[], P_solar_av_MW = Float64[])
    for h in 1:hours
        col = Symbol("h$h")
        for r in 1:nrow(dSolar)
            push!(solar_long, (
                hour = h,
                bus  = dSolar.bus[r],
                P_solar_av_MW = dSolar[r, col],))
        end
    end

    # base
    for row in 1:nrow(base)
        h = base.hour[row]
        b = base.bus[row]

        # Gas generation:
        idx_g = findfirst((solGen.hour .== h) .& (solGen.bus .== b))
        if idx_g !== nothing
            base.P_gen_MW[row] = solGen.powerGen[idx_g]
        end

        # Demand
        idx_d = findfirst((dDem.hour .== h) .& (dDem.bus_i .== b))
        if idx_d !== nothing
            base.Pd_MW[row] = dDem.Pd[idx_d]
        end

        # Curtailment
        idx_c = findfirst((solCurt.hour .== h) .& (solCurt.bus .== b))
        if idx_c !== nothing
            base.P_curt_solar_MW[row] = solCurt.P_curtailment[idx_c]
        end

        # Solar available
        idx_s = findfirst((solar_long.hour .== h) .& (solar_long.bus .== b))
        if idx_s !== nothing
            base.P_solar_av_MW[row] = solar_long.P_solar_av_MW[idx_s]
        end

        # Solar used
        base.P_solar_use_MW[row] = base.P_solar_av_MW[row] - base.P_curt_solar_MW[row]
    end

    sort!(base, [:hour, :bus])

    total_cost = solCosts[solCosts.hour .== -1, :operation_cost][1]

    # ==========================
    # Print BUS Tables
    # ==========================

    println("\n==================== RESULTS 3 NUDOS UC ====================\n")

    # BUS 1 (Solar)
    println(">>> BUS 1 (Solar)")
    bus1_df = filter(:bus => x -> x == 1, base)
    bus1_df = select(bus1_df, :hour, :P_solar_av_MW, :P_solar_use_MW, :P_curt_solar_MW)
    show(bus1_df, allrows = true, allcols = true)
    println("\n------------------------------------------------------------\n")

    # BUS 2 (Gas)
    println(">>> BUS 2 (Gas)")
    bus2_df = filter(:bus => x -> x == 2, base)
    bus2_df = select(bus2_df, :hour, :P_gen_MW)
    show(bus2_df, allrows = true, allcols = true)
    println("\n------------------------------------------------------------\n")

    # BUS 3 (Demand)
    println(">>> BUS 3 (Demand)")
    bus3_df = filter(:bus => x -> x == 3, base)
    bus3_df = select(bus3_df, :hour, :Pd_MW)
    show(bus3_df, allrows = true, allcols = true)
    println("\n------------------------------------------------------------\n")

    println("\nTotal system cost: $(round(total_cost, digits = 3)) €")
    println("============================================================\n")

    # ==========================
    # Save in Results/
    # ==========================

    results_folder = joinpath(@__DIR__, "..", "Results")
    if !isdir(results_folder)
        mkdir(results_folder)
    end

    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")

    CSV.write(joinpath(results_folder, "3NudosUC_bus1_$timestamp.csv"), bus1_df)
    CSV.write(joinpath(results_folder, "3NudosUC_bus2_$timestamp.csv"), bus2_df)
    CSV.write(joinpath(results_folder, "3NudosUC_bus3_$timestamp.csv"), bus3_df)

    println("Files saved in $results_folder")
    println(" - 3NudosUC_bus1_$timestamp.csv")
    println(" - 3NudosUC_bus2_$timestamp.csv")
    println(" - 3NudosUC_bus3_$timestamp.csv\n")

    return m, base, solGen, solFlows, solVoltage, solCosts, solCurt
end
