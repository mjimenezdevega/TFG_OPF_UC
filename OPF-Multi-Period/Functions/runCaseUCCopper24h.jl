# -----------------------------------------
# Functions/runCaseUCCopper24h.jl
# -----------------------------------------

include(joinpath(@__DIR__, "..", "Cases", "UC Copper 24h", "uc_copper_24h_csv.jl"))

using JuMP
const MOI = JuMP.MOI
using Dates
using CSV
using DataFrames

# Solver
map_solver(s::AbstractString) = begin
    s2 = lowercase(strip(s))
    if s2 == "gurobi"
        :Gurobi
    else
        :HiGHS
    end
end

"""
    run_case_24h(; solver_name="HiGHS") -> model

Run case UC+ED 24h "Copper Plate".
Lee:
- demandData.csv  (col: PD_MW)
- generatorData.csv
- solarData.csv   (col: PSolar_MW)
- windData.csv    (col: PWind_MW)
"""
function run_case_24h(; solver_name::AbstractString = "HiGHS")
    CASE   = joinpath(@__DIR__, "..", "Cases", "UC Copper 24h")
    demand = joinpath(CASE, "demandData.csv")
    gens   = joinpath(CASE, "generatorData.csv")
    solar  = joinpath(CASE, "solarData.csv")
    wind   = joinpath(CASE, "windData.csv")

    model, df_detail = run_example(
        demand_csv = demand,
        gens_csv   = gens,
        solar_csv  = solar,
        wind_csv   = wind,
        optimizer  = map_solver(solver_name),)

    
    save_results_detail!(df_detail)

    return model
end

"""
    save_results_detail!(df_detail)

Save in Results/CSV with:
hour, tech, Ng_on, Nst_cmd, Nsg_start, Nsd_stop, Pg_MW, PD_MW,
PRES_avail_MW, Pcurt_MW, cost_tech_EUR, cost_hour_total_EUR
"""
function save_results_detail!(df_detail::DataFrame)
    outdir = joinpath(@__DIR__, "..", "Results")
    isdir(outdir) || mkpath(outdir)

    fname = joinpath(outdir, "uc_copper_24h_DETAIL_" * Dates.format(now(), "yyyymmdd_HHMMSS") * ".csv")
    CSV.write(fname, df_detail)
    @info "Detailed results saved in $(fname)"
end