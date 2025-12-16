# ---------------------------------------------------------
# runCase19Nudos.jl
# ---------------------------------------------------------

using CSV
using DataFrames


include(joinpath(@__DIR__, "..", "OPF", "DC_OPF", "DC_OPF_UC.jl"))
using .DC_OPF_UC_Module: DC_OPF_UC

function run_case_19nudos(; solver_name::String = "Gurobi")

    case_dir = joinpath(@__DIR__, "..", "Cases", "19 Nudos")

    # Data
    dLine  = CSV.read(joinpath(case_dir, "lineData.csv"),      DataFrame; delim=';')
    dGen   = CSV.read(joinpath(case_dir, "generatorData.csv"), DataFrame; delim=';')
    dSolar = CSV.read(joinpath(case_dir, "solarData.csv"),     DataFrame; delim=';')
    dWind  = CSV.read(joinpath(case_dir, "windData.csv"),      DataFrame; delim=';')

    # Demand data: 1 file per hour (24 files)
    hours = 24
    dNodes = Vector{DataFrame}(undef, hours)
    demand_dir = joinpath(case_dir, "demandData")

    for h in 1:hours
        fname = joinpath(demand_dir, "nodeData_$(h).csv")
        dNodes[h] = CSV.read(fname, DataFrame; delim=';')
    end

    # Nodes and lines
    nN = maximum(dNodes[1].bus_i)
    nL = nrow(dLine)
    bMVA = 100

    # Run DC-OPF-UC
    m, solGen, solFlows, solVoltage, solCosts, solCurt =
        DC_OPF_UC(dLine, dGen, dNodes, nN, nL, bMVA, solver_name, hours, dSolar, dWind)

        # ================= TABLE HOUR–BUS =================

    # 1. Matrix: DEMAND / SOLAR / WIND in MW 
    nBuses = nN
    P_Demand_MW = zeros(nBuses, hours)
    for t in 1:hours
        for row in eachrow(dNodes[t])
            P_Demand_MW[row.bus_i, t] = row.Pd
        end
    end

    G_Solar_MW = zeros(nBuses, hours)
    for h in 1:hours
        col = Symbol("h$h")
        for r in 1:nrow(dSolar)
            bus = dSolar.bus[r]
            G_Solar_MW[bus, h] = dSolar[r, col]
        end
    end

    G_Wind_MW = zeros(nBuses, hours)
    for h in 1:hours
        col = Symbol("h$h")
        for r in 1:nrow(dWind)
            bus = dWind.bus[r]
            G_Wind_MW[bus, h] = dWind[r, col]
        end
    end

    # 2. DataFrame 
    tabla_19nudos = DataFrame(
        hour              = Int[],
        bus               = Int[],
        P_thermal_MW      = Float64[],
        P_solar_MW        = Float64[],
        P_wind_MW         = Float64[],
        P_demand_MW       = Float64[],
        P_curtailment_MW  = Float64[],)

    # 3. Table filling
    for t in 1:hours
        for i in 1:nBuses
            # Thermal generation bus-hour
            mask = (solGen.hour .== t) .& (solGen.bus .== i)
            P_thermal = sum(solGen.powerGen[mask])  # ya viene en MW porque multiplicaste por bMVA

            # Curtailment bus-hour
            mask_c = (solCurt.hour .== t) .& (solCurt.bus .== i)
            P_curt = isempty(solCurt.P_curtailment[mask_c]) ? 0.0 :
                     solCurt.P_curtailment[mask_c][1]  

            push!(tabla_19nudos, (
                hour             = t,
                bus              = i,
                P_thermal_MW     = P_thermal,
                P_solar_MW       = G_Solar_MW[i, t],
                P_wind_MW        = G_Wind_MW[i, t],
                P_demand_MW      = P_Demand_MW[i, t],
                P_curtailment_MW = P_curt,))
        end
    end

    # 4. Total cost (solCosts, hour = -1)
    total_cost = solCosts.operation_cost[findfirst(solCosts.hour .== -1)]

    println("\nTotal system cost = $(round(total_cost, digits=2)) €")

    # 5. Save CSV in Results
    results_dir = joinpath(@__DIR__, "..", "Results")

    CSV.write(joinpath(results_dir, "tabla_19nudos_generation.csv"), tabla_19nudos)
    CSV.write(joinpath(results_dir, "costs_19nudos.csv"), solCosts)


    return m, solGen, solFlows, solVoltage, solCosts, solCurt, tabla_19nudos

end
