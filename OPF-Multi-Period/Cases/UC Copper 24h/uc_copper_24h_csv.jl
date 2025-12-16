# ============================================================
# UC + ED "plato de cobre" (24 horas)
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
    load_data(; demand_csv="demand_24h.csv",
               gens_csv="generators_tech.csv",
               res_cf_csv="res_cf_24h.csv",
               res_meta_csv="renewables_meta.csv")

Devuelve un NamedTuple con:
(; T, PD, G, Ntotal, Pub, Pmin, Pmax, Cm, CnL, RR, Tst, Tmut, Tmdt, Ng0, Pg0,
   PRES, cFRES, VoLL)
"""
function load_data(; demand_csv::AbstractString="demand_24h.csv",
                     gens_csv::AbstractString="generators_tech.csv",
                     res_cf_csv::AbstractString="res_cf_24h.csv",
                     res_meta_csv::AbstractString="renewables_meta.csv")

    # --- Demand 24h ---
    dfPD = CSV.read(demand_csv, DataFrame)
    @assert :PD_MW ∈ propertynames(dfPD) "'$(demand_csv)' debe tener columna PD_MW"
    PD = Vector{Float64}(dfPD.PD_MW)
    T  = length(PD)
    @assert T == 24 "Se esperaban 24 filas (horas) en '$(demand_csv)'"

    # --- THERMAL GENERATORS ---
    dfG = CSV.read(gens_csv, DataFrame)
    required_cols = [:tech, :Ng_total, :Pub_MW, :Pmin_MW, :Pmax_MW,
                     :Cm_EUR_MWh, :CnL_EUR_h, :RR_MW_per_h, :Tst_h, :Tmut_h, :Tmdt_h,
                     :Ng0_units, :Pg0_total_MW]
    for c in required_cols
        @assert c ∈ propertynames(dfG) "'$(gens_csv)' debe incluir la columna $(c)"
    end

    G       = Symbol.(dfG.tech)
    Ntotal  = Dict(Symbol(r.tech) => Int(r.Ng_total)        for r in eachrow(dfG))
    Pub     = Dict(Symbol(r.tech) => Float64(r.Pub_MW)      for r in eachrow(dfG))
    Pmin    = Dict(Symbol(r.tech) => Float64(r.Pmin_MW)     for r in eachrow(dfG))
    Pmax    = Dict(Symbol(r.tech) => Float64(r.Pmax_MW)     for r in eachrow(dfG))
    Cm      = Dict(Symbol(r.tech) => Float64(r.Cm_EUR_MWh)  for r in eachrow(dfG))
    CnL     = Dict(Symbol(r.tech) => Float64(r.CnL_EUR_h)   for r in eachrow(dfG))
    RR      = Dict(Symbol(r.tech) => Float64(r.RR_MW_per_h) for r in eachrow(dfG))
    Tst     = Dict(Symbol(r.tech) => Int(r.Tst_h)           for r in eachrow(dfG))
    Tmut    = Dict(Symbol(r.tech) => Int(r.Tmut_h)          for r in eachrow(dfG))
    Tmdt    = Dict(Symbol(r.tech) => Int(r.Tmdt_h)          for r in eachrow(dfG))
    Ng0     = Dict(Symbol(r.tech) => Int(r.Ng0_units)       for r in eachrow(dfG))
    Pg0     = Dict(Symbol(r.tech) => Float64(r.Pg0_total_MW) for r in eachrow(dfG))

    # --- Renewables VoLL ---
    PRES  = 0.0
    VoLL  = 10_000.0
    if ispath(res_meta_csv)
        dfM = CSV.read(res_meta_csv, DataFrame)
        if nrow(dfM) > 0
            names_lc = lowercase.(string.(names(dfM)))

            # helper
            find_col(cands) = begin
                for c in cands
                    i = findfirst(==(c), names_lc)
                    if i !== nothing
                        return i
                    end
                end
                return nothing
            end

            iP = find_col(["pres_mw", "pres"])
            iV = find_col(["voll"])

            if iP !== nothing
                PRES = Float64(dfM[!, names(dfM)[iP]][1])
            end
            if iV !== nothing
                VoLL = Float64(dfM[!, names(dfM)[iV]][1])
            end
        end
    end

    cFRES = zeros(Float64, T)
    if ispath(res_cf_csv)
        dfCF = CSV.read(res_cf_csv, DataFrame)
        @assert :cf ∈ propertynames(dfCF) "'$(res_cf_csv)' debe tener columna cf"
        cFRES = Vector{Float64}(dfCF.cf)
        @assert length(cFRES) == T "'$(res_cf_csv)' debe tener $(T) filas (coincidir con demanda)"
    end

    return (; T, PD, G, Ntotal, Pub, Pmin, Pmax, Cm, CnL, RR, Tst, Tmut, Tmdt, Ng0, Pg0,
             PRES, cFRES, VoLL)
end

# ------------------------------------------------------------
# UC+ED Model
# ------------------------------------------------------------

function build_uc_copper_24h(data; optimizer::Symbol = :HiGHS)
    (; T, PD, G, Ntotal, Pub, Pmin, Pmax, Cm, CnL, RR,
       Tst, Tmut, Tmdt, Ng0, Pg0,
       PRES, cFRES, VoLL) = data

    model = Model()
    if optimizer == :Gurobi
        @assert _GUROBI_OK "Gurobi no está disponible en este proyecto."
        set_optimizer(model, Gurobi.Optimizer)
    else
        set_optimizer(model, HiGHS.Optimizer)
    end

    # ---------- VARIABLES ----------
    @variables(model, begin
        # Power per technology (MW)
        Pg[g in G, t in 1:T] >= 0.0

        # Committed units and Start-up commands
        0 <= Ng[g in G, t in 1:T]  <= Ntotal[g], Int    # unidades encendidas N_g(t)
        0 <= Nsg[g in G, t in 1:T] <= Ntotal[g], Int    # arranques efectivos N_g^{sg}(t)
        0 <= Nsd[g in G, t in 1:T] <= Ntotal[g], Int    # paradas efectivas N_g^{sd}(t)
        0 <= Nst[g in G, t in 1:T] <= Ntotal[g], Int    # órdenes de arranque N_g^{st}(t)

        # Slacks demand
        Ddef[t in 1:T]   >= 0.0
        Dspill[t in 1:T] >= 0.0

        # Curtail 
        Pcurt[t in 1:T]  >= 0.0
    end)

    # Dshed = Dspill - Ddef 
    @expression(model, Dshed[t in 1:T], Dspill[t] - Ddef[t])

    # ---------- OBJETIVE ----------
    εspill = 1e-3   # symbolic cost for conventional surplus
    @objective(model, Min,
        sum(
            sum(CnL[g]*Ng[g,t] + Cm[g]*Pg[g,t] for g in G) +
            VoLL*Ddef[t] + εspill*Dspill[t]
            for t in 1:T))

    # ---------- RESTRICTIONS ----------

    # (3.2) Power balance
    @constraint(model, [t in 1:T],
        sum(Pg[g,t] for g in G) + PRES*cFRES[t] - Pcurt[t] ==
        PD[t] - Ddef[t] + Dspill[t])

    # (3.3) Pmin / Pmax in MW per unit 
    @constraint(model, [g in G, t in 1:T],
        Pmin[g]*Ng[g,t] <= Pg[g,t])
    @constraint(model, [g in G, t in 1:T],
        Pg[g,t] <= Pmax[g]*Ng[g,t])

    # Curtailment 
    @constraint(model, [t in 1:T],
        Pcurt[t] <= PRES*cFRES[t])

    # (3.4) Ramp limits (t >= 2)

    # Hour 1 is left 'free' with respect to pre-horizon ramping limits, 
    # allowing units to enter directly with any output between Pmin and Pmax, 
    # thus enabling them to meet the initial demand.
    
    @constraint(model, [g in G, t in 2:T],
        Pg[g,t] - Pg[g,t-1] <=
            RR[g]*(Ng[g,t] - Nsg[g,t]) + Pmin[g]*Nsg[g,t])

    @constraint(model, [g in G, t in 2:T],
        Pg[g,t] - Pg[g,t-1] >=
            -Pmin[g]*Nsd[g,t] - RR[g]*(Ng[g,t] - Nsg[g,t]))

    # (3.5) Number of units
    
    @constraint(model, [g in G],
        Ng[g,1] == Ng0[g] + Nsg[g,1] - Nsd[g,1])
    # t >= 2:
    @constraint(model, [g in G, t in 2:T],
        Ng[g,t] == Ng[g,t-1] + Nsg[g,t] - Nsd[g,t])

    # (3.6) Startup time:
    
    for g in G, t in 2:T
        τ = t - Tst[g]
        if τ >= 1
            @constraint(model, Nsg[g,t] == Nst[g,τ])
        else  
        end
    end
    
    # (3.7) Minimum down time t >= 2
    # Nst(t) ≤ Ntotal - N_g(t-1) - Σ_{i=t-Tmdt}^{t} Nsd(i)
    for g in G, t in 2:T
        i_start = max(1, t - Tmdt[g])
        @constraint(model,
            Nst[g,t] <= Ntotal[g] -
                        Ng[g,t-1] -
                        sum(Nsd[g,i] for i in i_start:t))
    end

    # (3.8) Minimum up time t >= 2
    for g in G, t in 2:T
        i_start = max(1, t - Tmut[g])
        rng = i_start:(t-1)
        @constraint(model,
            Nsd[g,t] <= Ng[g,t-1] -
                        (isempty(rng) ? 0 : sum(Nsg[g,i] for i in rng)))
    end

    
    @constraint(model, [g in G, t in 1:T], 0 <= Ng[g,t] <= Ntotal[g])

    return model
end



# ------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------
function run_example(; demand_csv::AbstractString="demand_24h.csv",
                      gens_csv::AbstractString="generators_tech.csv",
                      res_cf_csv::AbstractString="res_cf_24h.csv",
                      res_meta_csv::AbstractString="renewables_meta.csv",
                      optimizer::Symbol=:HiGHS)

    data = load_data(demand_csv=demand_csv, gens_csv=gens_csv,
                     res_cf_csv=res_cf_csv, res_meta_csv=res_meta_csv)

    model = build_uc_copper_24h(data; optimizer=optimizer)
    optimize!(model)

    println("Estatus: ", termination_status(model))

    if !has_values(model)
        println("The model has no solution.")
        return model, nothing
    end

    total_cost = objective_value(model)

    println("\n====================================")
    @printf("   Total cost 24h = %.2f €\n", total_cost)
    println("====================================\n")

    T = data.T
    G = data.G
    eps_spill = 1e-3

    # ------- Results Table -------
    df = DataFrame(
        hora                 = Int[],
        tech                 = String[],
        Ng_on                = Int[],
        Nsg_start            = Int[],
        Nsd_stop             = Int[],
        Pg_MW                = Float64[],
        PD_MW                = Float64[],
        PRES_avail_MW        = Float64[],
        Pcurt_MW             = Float64[],
        cost_EUR             = Float64[],
        cost_hora_total_EUR  = Float64[],)

    for t in 1:T
        PD_t    = data.PD[t]
        PRES_t  = data.PRES * data.cFRES[t]
        Pcurt_t = value(model[:Pcurt][t])

        # Total hour cost
        cost_hora = sum(
            data.CnL[g] * value(model[:Ng][g,t]) +
            data.Cm[g]  * value(model[:Pg][g,t])
            for g in G) + data.VoLL * value(model[:Ddef][t]) +
            eps_spill   * value(model[:Dspill][t])

        for g in G
            Ng_on     = Int(round(value(model[:Ng][g,t])))
            Nsg_start = Int(round(value(model[:Nsg][g,t])))
            Nsd_stop  = Int(round(value(model[:Nsd][g,t])))
            Pg_MW     = value(model[:Pg][g,t])

            cost_gen = data.CnL[g]*Ng_on + data.Cm[g]*Pg_MW

            push!(df, (
                t,
                String(g),
                Ng_on,
                Nsg_start,
                Nsd_stop,
                Pg_MW,
                PD_t,
                PRES_t,
                Pcurt_t,
                cost_gen,
                cost_hora))
        end
    end

    println("===== RESULTS =====")
    show(df, allrows=true, allcols=true)
    println("\n===============================")

    return model, df
end

if abspath(PROGRAM_FILE) == @__FILE__
    CASE = @__DIR__
    run_example(
        demand_csv   = joinpath(CASE, "demand_24h.csv"),
        gens_csv     = joinpath(CASE, "generators_tech.csv"),
        res_cf_csv   = joinpath(CASE, "res_cf_24h.csv"),
        res_meta_csv = joinpath(CASE, "renewables_meta.csv"),
        optimizer    = :HiGHS,)
end
