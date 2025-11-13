# -----------------------------------------
# Functions/runCase24h.jl
# Conecta el caso 24h (UC-Copper) con el menú
# -----------------------------------------

# Trae run_example y load_data del caso 24h
include(joinpath(@__DIR__, "..", "Cases", "Caso 24h", "uc_copper_24h_csv.jl"))

using JuMP
const MOI = JuMP.MOI
using Dates
using CSV
using DataFrames


# Convierte el texto del menú a símbolo de solver
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

Lanza el caso UC+ED 24h "plato de cobre" leyendo los CSV de `Cases/Caso 24h`.
`solver_name` admite "HiGHS" o "Gurobi". Devuelve el modelo JuMP ya optimizado.
"""
function run_case_24h(; solver_name::AbstractString = "HiGHS")
    CASE  = joinpath(@__DIR__, "..", "Cases", "Caso 24h")
    demand = joinpath(CASE, "demand_24h.csv")
    gens   = joinpath(CASE, "generators_tech.csv")
    cf     = joinpath(CASE, "res_cf_24h.csv")
    meta   = joinpath(CASE, "renewables_meta.csv")

    m = run_example(
        demand_csv   = demand,
        gens_csv     = gens,
        res_cf_csv   = cf,
        res_meta_csv = meta,
        optimizer    = map_solver(solver_name),
    )

    # Guardado opcional y robusto (solo si hay solución)
    try
        data = load_data(demand_csv=demand, gens_csv=gens, res_cf_csv=cf, res_meta_csv=meta)
        save_results!(m, data)
    catch err
        @warn "No se pudieron guardar resultados: $err"
    end

    return m
end

# ----------------- helpers -----------------

"""
    save_results!(model, data)

Escribe un CSV en `Results/` con (hour, tech, Pg, Ng).
Solo guarda si el estado del solver es factible/óptimo.
"""
function save_results!(model, data)
    ts = JuMP.termination_status(model)
    if ts ∉ (MOI.OPTIMAL, MOI.FEASIBLE_POINT, MOI.LOCALLY_SOLVED)
        @warn "Estado del modelo: $ts. No hay solución primal; no se guardan resultados."
        return
    end

    (; G, T) = data
    rows = DataFrame(hour=Int[], tech=String[], Pg=Float64[], Ng=Union{Missing,Int}[])
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
    @info "Resultados guardados en $(fname)"
end
