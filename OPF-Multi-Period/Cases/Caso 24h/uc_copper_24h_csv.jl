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
    Ntotal  = Dict(Symbol(r.tech) => Int(r.Ng_total)      for r in eachrow(dfG))
    Pub     = Dict(Symbol(r.tech) => Float64(r.Pub_MW)    for r in eachrow(dfG))
    Pmin    = Dict(Symbol(r.tech) => Float64(r.Pmin_MW)   for r in eachrow(dfG))
    Pmax    = Dict(Symbol(r.tech) => Float64(r.Pmax_MW)   for r in eachrow(dfG))
    Cm      = Dict(Symbol(r.tech) => Float64(r.Cm_EUR_MWh) for r in eachrow(dfG))
    CnL     = Dict(Symbol(r.tech) => Float64(r.CnL_EUR_h)  for r in eachrow(dfG))
    RR      = Dict(Symbol(r.tech) => Float64(r.RR_MW_per_h) for r in eachrow(dfG))
    Tst     = Dict(Symbol(r.tech) => Int(r.Tst_h)         for r in eachrow(dfG))
    Tmut    = Dict(Symbol(r.tech) => Int(r.Tmut_h)        for r in eachrow(dfG))
    Tmdt    = Dict(Symbol(r.tech) => Int(r.Tmdt_h)        for r in eachrow(dfG))
    Ng0     = Dict(Symbol(r.tech) => Int(r.Ng0_units)     for r in eachrow(dfG))
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
"""
    build_uc_copper_24h(data; optimizer::Symbol = :HiGHS)

Construye el modelo JuMP del UC+ED 24h (plato de cobre) con:
- Dshed simétrico: Ddef (déficit) y Dspill (excedente), alias Dshed = Dspill - Ddef
- Ng, Nsg, Nsd, Nst enteros y acotados
"""
function build_uc_copper_24h(data; optimizer::Symbol = :HiGHS)
    (; T, PD, G, Ntotal, Pub, Pmin, Pmax, Cm, CnL, RR, Tst, Tmut, Tmdt, Ng0, Pg0,
      PRES, cFRES, VoLL) = data

    model = Model()
    if optimizer == :Gurobi
        @assert _GUROBI_OK "Gurobi no está disponible en este proyecto (Pkg.add(\"Gurobi\"); configurar licencia)."
        set_optimizer(model, Gurobi.Optimizer)
    else
        set_optimizer(model, HiGHS.Optimizer)
    end

    # ---------- Variables ----------
    @variables(model, begin
        # Power per Technology
        Pg[g in G, t in 1:T] >= 0.0

        # Unidades encendidas, arranques, paradas y órdenes diferidas 
        0 <= Ng[g in G, t in 1:T]  <= Ntotal[g],   Int
        0 <= Nsg[g in G, t in 1:T] <= Ntotal[g],   Int
        0 <= Nsd[g in G, t in 1:T] <= Ntotal[g],   Int
        0 <= Nst[g in G, t in 1:T] <= Ntotal[g],   Int

        # Slack simétrico de demanda
        Ddef[t in 1:T]   >= 0.0     # déficit (load shedding)
        Dspill[t in 1:T] >= 0.0     # excedente (vertedero convencional)

        # Curtail renovable
        Pcurt[t in 1:T]  >= 0.0
    end)

    # Alias con signo (para logs/impresiones si quieres seguir “viendo” Dshed)
    @expression(model, Dshed[t in 1:T], Dspill[t] - Ddef[t])

    # ---------- Objetivo ----------
    εspill = 1e-3  # coste simbólico por “tirar” excedente convencional
    @objective(model, Min,
        sum(
            sum(CnL[g]*Ng[g,t] + Cm[g]*Pg[g,t] for g in G) +
            VoLL*Ddef[t] + εspill*Dspill[t]
            for t in 1:T
        )
    )

    # ---------- RESTRICCIONES ----------

    # (a) Balance de potencia con renovable agregada, curtail y slacks
    @constraint(model, [t in 1:T],
        sum(Pg[g,t] for g in G) + PRES*cFRES[t] - Pcurt[t] == PD[t] - Ddef[t] + Dspill[t]
    )

    # (b) Pmin / Pmax escalado por Ng
    @constraint(model, [g in G, t in 1:T], Pmin[g]*Pub[g]*Ng[g,t] <= Pg[g,t])
    @constraint(model, [g in G, t in 1:T], Pg[g,t] <= Pmax[g]*Pub[g]*Ng[g,t])

    # (c) Curtailment ≤ renovable disponible
    @constraint(model, [t in 1:T], Pcurt[t] <= PRES*cFRES[t])

    # (d) Rampas con start/stop (RR en MW/h por unidad; Pub[g] MW por unidad)
    # t >= 2
    @constraint(model, [g in G, t in 2:T],
        Pg[g,t] - Pg[g,t-1] <= RR[g]*Ng[g,t-1] + Pub[g]*Nsg[g,t]
    )
    @constraint(model, [g in G, t in 2:T],
        Pg[g,t-1] - Pg[g,t] <= RR[g]*Ng[g,t]   + Pub[g]*Nsd[g,t]
    )
    # t = 1 respecto a Pg0, Ng0
    @constraint(model, [g in G],
        Pg[g,1] - Pg0[g] <= RR[g]*Ng0[g] + Pub[g]*Nsg[g,1]
    )
    @constraint(model, [g in G],
        Pg0[g] - Pg[g,1] <= RR[g]*Ng[g,1] + Pub[g]*Nsd[g,1]
    )

    # (e) Evolución del número de unidades encendidas
    @constraint(model, [g in G, t in 1:T],
        Ng[g,t] == (t==1 ? Ng0[g] : Ng[g,t-1]) + Nsg[g,t] - Nsd[g,t]
    )

    # (g) Arranque diferido: Nsg(t) = Nst(t - Tst)
    for g in G, t in 1:T
        τ = t - Tst[g]
        if τ >= 1
            @constraint(model, Nsg[g,t] == Nst[g,τ])
        else
            @constraint(model, Nsg[g,t] == 0)
        end
    end

    # (h) Minimum down time
    for g in G, t in 1:T
        rng = max(1, t - Tmdt[g]) : t
        # no se puede arrancar si no se ha cumplido down-time desde el último apagado
        @constraint(model,
            Nst[g,t] <= Ntotal[g] - (t==1 ? Ng0[g] : Ng[g,t-1]) - sum(Nsd[g,i] for i in rng)
        )
    end

    # (i) Minimum up time
    for g in G, t in 1:T
        rng = max(1, t - Tmut[g]) : t-1
        @constraint(model,
            Nsd[g,t] <= (t==1 ? Ng0[g] : Ng[g,t-1]) - (isempty(rng) ? 0 : sum(Nsg[g,i] for i in rng))
        )
    end

    # Lógica de límites del número de unidades
    @constraint(model, [g in G, t in 1:T], 0 <= Ng[g,t] <= Ntotal[g])

    return model
end

# ------------------------------------------------------------
# Ejecutable “de ejemplo” 
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
    if has_values(model)
        @printf("Coste total (€): %.2f\n", objective_value(model))
        for t in 1:data.T
            println("\n=== Hora ", t, " ===")
            @printf("PD=%.1f MW | RES=%.1f MW\n", data.PD[t], data.PRES*data.cFRES[t])
            for g in data.G
                @printf("%s: Pg=%.1f MW, Ng=%d, Nsg=%d, Nsd=%d\n",
                        String(g),
                        value(model[:Pg][g,t]),
                        Int(round(value(model[:Ng][g,t]))),
                        Int(round(value(model[:Nsg][g,t]))),
                        Int(round(value(model[:Nsd][g,t]))))
            end
            @printf("Pcurt=%.1f MW, Dshed=%.1f MW\n",
                    value(model[:Pcurt][t]),
                    value(model[:Dshed][t]))
        end
    end
    return model
end

if abspath(PROGRAM_FILE) == @__FILE__
    CASE = @__DIR__            # carpeta Cases/Caso 24h
    run_example(
        demand_csv = joinpath(CASE, "demand_24h.csv"),
        gens_csv   = joinpath(CASE, "generators_tech.csv"),
        res_cf_csv = joinpath(CASE, "res_cf_24h.csv"),
        res_meta_csv = joinpath(CASE, "renewables_meta.csv"),
        optimizer  = :HiGHS,
    )
end



