# ============================================================
# UC + ED "Copper plate" (24 hours) 
# ============================================================

using JuMP
using CSV
using DataFrames
using Printf

# Solvers
using HiGHS
const _GUROBI_OK = try
    @eval using Gurobi
    true
catch
    false
end

# ------------------------------------------------------------
# Data from CSV
# ------------------------------------------------------------
"""
    load_data(; demand_csv="demandData.csv",
               gens_csv="generatorData.csv",
               solar_csv="solarData.csv",
               wind_csv="windData.csv")

NamedTuple:
(; T, PD, G, Ntotal, Pub, Pmin, Pmax, Cm, CnL, Cstart, RR, Tst, Tmut, Tmdt, Ng0, Pg0,
   PRES_avail)
"""
function load_data(; demand_csv::AbstractString="demandData.csv",
                     gens_csv::AbstractString="generatorData.csv",
                     solar_csv::AbstractString="solarData.csv",
                     wind_csv::AbstractString="windData.csv")

    # --- Demand 24h ---
    dfPD = CSV.read(demand_csv, DataFrame)
    if !(:PD_MW in propertynames(dfPD))
        error("ERROR in demandData.csv:\nColumn PD_MW must exist.\nColumns: $(String.(propertynames(dfPD)))")
    end
    PD = Vector{Float64}(dfPD.PD_MW)
    T  = length(PD)
    @assert T == 24 "24 lines were expected (hours) '$(demand_csv)'"

    # --- THERMAL GENERATORS (per group) ---
    dfG = CSV.read(gens_csv, DataFrame)
    required_cols = [:tech, :Ng_total, :Pub_MW, :Pmin_MW, :Pmax_MW,
                     :Cm_EUR_MWh, :CnL_EUR_h, :C_start_EUR,
                     :RR_MW_per_h, :Tst_h, :Tmut_h, :Tmdt_h,
                     :Ng0_units, :Pg0_total_MW]
    for c in required_cols
        @assert c ∈ propertynames(dfG) "'$(gens_csv)' must have column $(c)"
    end

    G       = Symbol.(dfG.tech)
    Ntotal  = Dict(Symbol(r.tech) => Int(r.Ng_total)         for r in eachrow(dfG))
    Pub     = Dict(Symbol(r.tech) => Float64(r.Pub_MW)       for r in eachrow(dfG))
    Pmin    = Dict(Symbol(r.tech) => Float64(r.Pmin_MW)      for r in eachrow(dfG))
    Pmax    = Dict(Symbol(r.tech) => Float64(r.Pmax_MW)      for r in eachrow(dfG))
    Cm      = Dict(Symbol(r.tech) => Float64(r.Cm_EUR_MWh)   for r in eachrow(dfG))
    CnL     = Dict(Symbol(r.tech) => Float64(r.CnL_EUR_h)    for r in eachrow(dfG))
    Cstart  = Dict(Symbol(r.tech) => Float64(r.C_start_EUR)  for r in eachrow(dfG))
    RR      = Dict(Symbol(r.tech) => Float64(r.RR_MW_per_h)  for r in eachrow(dfG))
    Tst     = Dict(Symbol(r.tech) => Int(r.Tst_h)            for r in eachrow(dfG))
    Tmut    = Dict(Symbol(r.tech) => Int(r.Tmut_h)           for r in eachrow(dfG))
    Tmdt    = Dict(Symbol(r.tech) => Int(r.Tmdt_h)           for r in eachrow(dfG))
    Ng0     = Dict(Symbol(r.tech) => Int(r.Ng0_units)        for r in eachrow(dfG))
    Pg0     = Dict(Symbol(r.tech) => Float64(r.Pg0_total_MW) for r in eachrow(dfG))

    # --- Renewables per hour (MW) ---
    @assert ispath(solar_csv) "'$(solar_csv)' not found"
    @assert ispath(wind_csv)  "'$(wind_csv)' not found"

    dfS = CSV.read(solar_csv, DataFrame)
    @assert :PSolar_MW ∈ propertynames(dfS) "'$(solar_csv)' must have PSolar_MW column"
    PSolar = Vector{Float64}(dfS.PSolar_MW)
    @assert length(PSolar) == T "'$(solar_csv)' must have $(T) (compatible with demand)"

    dfW = CSV.read(wind_csv, DataFrame)
    @assert :PWind_MW ∈ propertynames(dfW) "'$(wind_csv)' must have PWind_MW column"
    PWind = Vector{Float64}(dfW.PWind_MW)
    @assert length(PWind) == T "'$(wind_csv)' must have $(T) (compatible with demand)"

    PRES_avail = PSolar .+ PWind

    return (; T, PD, G, Ntotal, Pub, Pmin, Pmax, Cm, CnL, Cstart, RR, Tst, Tmut, Tmdt, Ng0, Pg0, PRES_avail)
end

# ------------------------------------------------------------
# UC+ED Model
# ------------------------------------------------------------
function build_uc_copper_24h(data; optimizer::Symbol = :HiGHS)
    (; T, PD, G, Ntotal, Pmin, Pmax, Cm, CnL, Cstart, RR,
       Tst, Tmut, Tmdt, Ng0,
       PRES_avail) = data

    model = Model()
    if optimizer == :Gurobi
        @assert _GUROBI_OK "Gurobi is not available in this project."
        set_optimizer(model, Gurobi.Optimizer)
    else
        set_optimizer(model, HiGHS.Optimizer)
    end

    # Slack costs
    VoLL   = 10_000.0
    εspill = 1e-3

    @variables(model, begin
        Pg[g in G, t in 1:T] >= 0.0

        0 <= Ng[g in G, t in 1:T]  <= Ntotal[g], Int
        0 <= Nsg[g in G, t in 1:T] <= Ntotal[g], Int
        0 <= Nsd[g in G, t in 1:T] <= Ntotal[g], Int
        0 <= Nst[g in G, t in 1:T] <= Ntotal[g], Int

        Pcurt[t in 1:T]  >= 0.0

        # Slacks (to avoid infeasibility)
        Ddef[t in 1:T]   >= 0.0  # load shedding
        Dspill[t in 1:T] >= 0.0  # thermal spill
    end)

    # Objective (same structure as DC-OPF-UC + slacks)
    @objective(model, Min,
        sum(
            sum(Cm[g]*Pg[g,t] + CnL[g]*Ng[g,t] + (t == 1 ? 0.0 : Cstart[g]*Nsg[g,t]) for g in G)
            + VoLL*Ddef[t] + εspill*Dspill[t]
            for t in 1:T))

    # Balance (with slacks)
    @constraint(model, [t in 1:T],
        sum(Pg[g,t] for g in G) + PRES_avail[t] - Pcurt[t]
        ==
        PD[t] - Ddef[t] + Dspill[t])

    # Pmin/Pmax
    @constraint(model, [g in G, t in 1:T], Pmin[g]*Ng[g,t] <= Pg[g,t])
    @constraint(model, [g in G, t in 1:T], Pg[g,t] <= Pmax[g]*Ng[g,t])

    # Curtailment (only renewables)
    @constraint(model, [t in 1:T], Pcurt[t] <= PRES_avail[t])

    # ------------------------------
    # Intertemporal constraints
    # ------------------------------

    # Ramps apply from hour 2
    @constraint(model, [g in G, t in 2:T],
        Pg[g,t] - Pg[g,t-1] <= RR[g]*(Ng[g,t] - Nsg[g,t]) + Pmin[g]*Nsg[g,t])
    @constraint(model, [g in G, t in 2:T],
        Pg[g,t] - Pg[g,t-1] >= -Pmin[g]*Nsd[g,t] - RR[g]*(Ng[g,t] - Nsg[g,t]))

    # Unit count dynamics
    @constraint(model, [g in G],
        Ng[g,1] == Ng0[g] + Nsg[g,1] - Nsd[g,1])
    @constraint(model, [g in G, t in 2:T],
        Ng[g,t] == Ng[g,t-1] + Nsg[g,t] - Nsd[g,t])

    # --- Start generating mapping ---
    # Hour 1 is "free": allow Nsg[g,1] = Nst[g,1] even if Tst[g] > 0
    @constraint(model, [g in G], Nsg[g,1] == Nst[g,1])

    # For t >= 2: apply Tst delay (if τ<1 -> no start possible)
    for g in G
        if Tst[g] > 0
            for t in 2:T
                τ = t - Tst[g]
                if τ >= 1
                    @constraint(model, Nsg[g,t] == Nst[g,τ])
                else
                    @constraint(model, Nsg[g,t] == 0)
                end
            end
        else
            # If Tst == 0, also map directly for all hours
            @constraint(model, [t in 2:T], Nsg[g,t] == Nst[g,t])
        end
    end

    # Minimum down time
    for g in G
        if Tmdt[g] > 0
            for t in 2:T
                i_start = max(1, t - Tmdt[g] + 1)
                @constraint(model,
                    Nst[g,t] <= Ntotal[g] - Ng[g,t-1] - sum(Nsd[g,i] for i in i_start:(t-1)))
            end
        end
    end

    # Minimum up time
    for g in G
        if Tmut[g] > 0
            for t in 2:T
                i_start = max(1, t - Tmut[g] + 1)
                @constraint(model,
                    Nsd[g,t] <= Ng[g,t-1] - sum(Nsg[g,i] for i in i_start:(t-1)))
            end
        end
    end

    return model
end

# ------------------------------------------------------------
# RUN + RESULTS
# ------------------------------------------------------------
function run_example(; demand_csv::AbstractString="demandData.csv",
                      gens_csv::AbstractString="generatorData.csv",
                      solar_csv::AbstractString="solarData.csv",
                      wind_csv::AbstractString="windData.csv",
                      optimizer::Symbol=:HiGHS)

    data = load_data(demand_csv=demand_csv, gens_csv=gens_csv,
                     solar_csv=solar_csv, wind_csv=wind_csv)

    model = build_uc_copper_24h(data; optimizer=optimizer)
    optimize!(model)

    println("Estatus: ", termination_status(model))
    if !has_values(model)
        println("The model has no solution.")
        return model, nothing
    end

    total_cost = objective_value(model)
    println("\n====================================")
    @printf("   Total cost 24h UC Copper= %.2f €\n", total_cost)
    println("====================================\n")

    T = data.T
    G = data.G

    df = DataFrame(
        hour                = Int[],
        tech                = String[],
        Ng_on               = Int[],
        Nst_cmd             = Int[],
        Nsg_start           = Int[],
        Nsd_stop            = Int[],
        Pg_MW               = Float64[],
        PD_MW               = Float64[],
        PRES_avail_MW       = Float64[],
        Pcurt_MW            = Float64[],
        Ddef_MW             = Float64[],
        Dspill_MW           = Float64[],
        cost_tech_EUR       = Float64[],
        cost_hour_total_EUR = Float64[],)

    # Slack costs
    VoLL   = 10_000.0
    εspill = 1e-3

    for t in 1:T
        PD_t    = data.PD[t]
        PRES_t  = data.PRES_avail[t]
        Pcurt_t = value(model[:Pcurt][t])
        Ddef_t  = value(model[:Ddef][t])
        Dsp_t   = value(model[:Dspill][t])

        # Total Cost hour (match objective: NO startup cost at t=1)
        cost_hour = 0.0
        for g in G
            cost_hour += data.Cm[g]*value(model[:Pg][g,t]) +
                         data.CnL[g]*value(model[:Ng][g,t]) +
                         (t == 1 ? 0.0 : data.Cstart[g]*value(model[:Nsg][g,t]))
        end
        cost_hour += VoLL*Ddef_t + εspill*Dsp_t

        for g in G
            Ng_on   = Int(round(value(model[:Ng][g,t])))
            Nst_cmd = Int(round(value(model[:Nst][g,t])))
            Nsg     = Int(round(value(model[:Nsg][g,t])))
            Nsd     = Int(round(value(model[:Nsd][g,t])))
            Pg_MW   = value(model[:Pg][g,t])

            # Per-tech cost (also NO startup cost at t=1)
            cost_tech = data.Cm[g]*Pg_MW + data.CnL[g]*Ng_on + (t == 1 ? 0.0 : data.Cstart[g]*Nsg)

            push!(df, (
                t, String(g),
                Ng_on, Nst_cmd, Nsg, Nsd,
                Pg_MW, PD_t, PRES_t, Pcurt_t, Ddef_t, Dsp_t,
                cost_tech, cost_hour))
        end
    end

    println("===== RESULTS (detailed) =====")
    show(df, allrows=true, allcols=true)
    println("\n==============================")

    return model, df
end