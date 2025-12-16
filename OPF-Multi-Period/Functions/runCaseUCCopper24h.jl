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
    run_case_UC_Copper_24h(; solver_name="HiGHS") -> model

Run case UC+ED 24h "plato de cobre".
`solver_name` admits "HiGHS" or "Gurobi". 
"""
function run_case_24h(; solver_name::AbstractString = "HiGHS")
    CASE  = joinpath(@__DIR__, "..", "Cases", "UC Copper 24h")
    demand = joinpath(CASE, "demandData.csv")
    gens   = joinpath(CASE, "generatorData.csv")
    cf     = joinpath(CASE, "res_cf_24h.csv")
    meta   = joinpath(CASE, "renewables_meta.csv")

    # run_example (model, df)
    model, df = run_example(
        demand_csv   = demand,
        gens_csv     = gens,
        res_cf_csv   = cf,
        res_meta_csv = meta,
        optimizer    = map_solver(solver_name),
    )

    # Save results
    data = load_data(demand_csv=demand, gens_csv=gens, res_cf_csv=cf, res_meta_csv=meta)
    save_results!(model, data)

    return model
end

"""
    save_results!(model, data)

Write a CSV in `Results/` with (hour, tech, Pg, Ng).
"""
function save_results!(model, data)
    (; G, T) = data

    rows = DataFrame(
        hour = Int[],
        tech = String[],
        Pg   = Float64[],
        Ng   = Union{Missing,Int}[],
    )

    for t in 1:T, g in G
        pg = JuMP.value(model[:Pg][g, t])
        ng = JuMP.value(model[:Ng][g, t])

        pg = isnan(pg) ? 0.0 : pg
        ng_out = isnan(ng) ? missing : Int(round(ng))

        push!(rows, (t, String(g), pg, ng_out))
    end

    outdir = joinpath(@__DIR__, "..", "Results")
    isdir(outdir) || mkpath(outdir)
    fname = joinpath(outdir, "uc_copper_24h_" * Dates.format(now(), "yyyymmdd_HHMMSS") * ".csv")
    CSV.write(fname, rows)
    @info "Results saved in $(fname)"
end
