# ---------------------------------------------------------
# OPF/DC_OPF/DC_OPF_UC.jl — DC-OPF with Unit Commitment
# ---------------------------------------------------------

module DC_OPF_UC_Module

using JuMP
using DataFrames
using LinearAlgebra
using SparseArrays
import MathOptInterface as MOI

import Gurobi
import HiGHS
import Ipopt

include(joinpath(@__DIR__, "Functions", "susceptanceMatrix.jl"))

function DC_OPF_UC(
    dLine::DataFrame,
    dGen::DataFrame,
    dNodes::Vector{DataFrame},
    nN::Int,
    nL::Int,
    bMVA::Int,
    solver::String,
    hours::Int,
    dSolar::DataFrame,
    dWind::DataFrame,)

    # === Thermal Generators Data ===
    nG      = nrow(dGen)
    gen_bus = Vector{Int}(dGen.bus)

    Pmin_pu = dGen.Pmin_MW      ./ bMVA
    Pmax_pu = dGen.Pmax_MW      ./ bMVA
    Cm      = dGen.Cm_EUR_MWh
    CnL     = dGen.CnL_EUR_h
    Cstart  = dGen.C_start_EUR

    RR_pu   = dGen.RR_MW_per_h  ./ bMVA          
    Pmsg_pu = dGen.Pmsg_MW      ./ bMVA          

    Tmut    = dGen.Tmut_h       # minimum up time
    Tmdt    = dGen.Tmdt_h       # minimum down time
    Tst     = dGen.Tst_h        # Start generator time

 
    # === Susceptance Matrix ===
    B = susceptanceMatrix(dLine, nN, nL)

    # === Demands ===
    P_Demand = zeros(nN, hours)
    for t in 1:hours
        for row in eachrow(dNodes[t])
            P_Demand[row.bus_i, t] = row.Pd / bMVA
        end
    end

    # === Solar ===
    G_Solar = zeros(nN, hours)
    for h in 1:hours
        col = Symbol("h$h")
        for r in 1:nrow(dSolar)
            bus = dSolar.bus[r]
            G_Solar[bus, h] = dSolar[r, col] / bMVA
        end
    end

    # === Wind ===
    G_Wind = zeros(nN, hours)
    for h in 1:hours
        col = Symbol("h$h")
        for r in 1:nrow(dWind)
            bus = dWind.bus[r]
            G_Wind[bus, h] = dWind[r, col] / bMVA
        end
    end

    # === Model ===
    m = if solver == "Gurobi"
        Model(Gurobi.Optimizer)
    elseif solver == "HiGHS"
        Model(HiGHS.Optimizer)
    elseif solver == "Ipopt"
        Model(Ipopt.Optimizer)
    else
        error("Solver no soportado en DC_OPF_UC")
    end
    set_silent(m)

    # --------------------------------------------------
    # VARIABLES
    # --------------------------------------------------
    @variable(m, P_g[1:nG, 1:hours] >= 0)
    @variable(m, θ[1:nN, 1:hours])
    @variable(m, P_Curt[1:nN, 1:hours] >= 0)

    @variable(m, y_g[1:nG, 1:hours],   Bin)  # ON/OFF
    @variable(m, y_sg[1:nG, 1:hours], Bin)   # y_g^{sg}(t) - start generating
    @variable(m, y_sd[1:nG, 1:hours], Bin)   # y_g^{sd}(t) - shutdown
    @variable(m, y_st[1:nG, 1:hours], Bin)   # y_g^{st}(t) - startup (order)

   
    @constraint(m, [g in 1:nG, t in 1:hours], y_sg[g,t] + y_sd[g,t] <= 1)

    # --------------------------------------------------
    # A.3 – Power limits
    # Pmin * y_g(t) ≤ P_g(t) ≤ Pmax * y_g(t)
    # --------------------------------------------------
    @constraint(m, [g in 1:nG, t in 1:hours],
        P_g[g,t] >= Pmin_pu[g] * y_g[g,t])
    @constraint(m, [g in 1:nG, t in 1:hours],
        P_g[g,t] <= Pmax_pu[g] * y_g[g,t])

    # --------------------------------------------------
    # A.5 – The commitment state of a unit is 'on' if it was generating in the previous time-step or has started 
    #       generating in the current time-step, unless it has been shut down in the current time-step:
    #
    # y_g(t) = y_g(t-1) + y_sg(t) - y_sd(t)
    # --------------------------------------------------
    @constraint(m, [g in 1:nG, t in 2:hours],
        y_g[g,t] == y_g[g,t-1] + y_sg[g,t] - y_sd[g,t])

    # --------------------------------------------------
    # A.4 – Ramps
    #
    # -Pmsg * y_sd(t) - Prr * y(t-1) ≤ P(t)-P(t-1)
    #                                ≤ Prr * y(t) + Pmsg * y_sg(t)
    # --------------------------------------------------
    @constraint(m, [g in 1:nG],
        P_g[g,1] <= RR_pu[g]*y_g[g,1] + Pmsg_pu[g]*y_sg[g,1])

    @constraint(m, [g in 1:nG, t in 2:hours],
        -(Pmsg_pu[g]*y_sd[g,t] + RR_pu[g]*y_g[g,t-1])
        <= P_g[g,t] - P_g[g,t-1])

    @constraint(m, [g in 1:nG, t in 2:hours],
        P_g[g,t] - P_g[g,t-1]
        <= RR_pu[g]*y_g[g,t] + Pmsg_pu[g]*y_sg[g,t])

    # --------------------------------------------------
    # A.6 – Start time Tst_g
    #
    # y_sg(t) = y_st(t - Tst_g)
    # --------------------------------------------------
    for g in 1:nG
        if Tst[g] > 0
            for t in 1:hours
                if t > Tst[g]
                    @constraint(m,
                        y_sg[g,t] == y_st[g, t - Tst[g]])
                end
            end
        end
    end

    # --------------------------------------------------
    # A.7 – Minimum down time
    #
    # y_st(t) ≤ [1 - y_g(t-1)] - Σ_{i} y_sd(i)
    # --------------------------------------------------
    for g in 1:nG
        if Tmdt[g] > 0
            for t in 2:hours
                first_t_down = max(1, t - Tmdt[g] + 1)
                @constraint(m,
                    y_st[g,t] <=
                        (1 - y_g[g,t-1]) -
                        sum(y_sd[g,τ] for τ in first_t_down:(t-1)))
            end
        end
    end

    # --------------------------------------------------
    # A.8 – Minimum up time
    #
    # y_sd(t) ≤ y_g(t-1) - Σ_{i} y_sg(i)
    # --------------------------------------------------
    for g in 1:nG
        if Tmut[g] > 0
            for t in 2:hours
                first_t_up = max(1, t - Tmut[g] + 1)
                @constraint(m,
                    y_sd[g,t] <=
                        y_g[g,t-1] -
                        sum(y_sg[g,τ] for τ in first_t_up:(t-1)))
            end
        end
    end

    # --------------------------------------------------
    # Nodal balance (DC)
    # --------------------------------------------------
    for t in 1:hours
        for i in 1:nN
            @constraint(m,
                sum(P_g[g,t] for g in 1:nG if gen_bus[g] == i) +
                (G_Solar[i,t] + G_Wind[i,t] - P_Curt[i,t]) -
                P_Demand[i,t]
                ==
                sum(B[i,j]*(θ[i,t] - θ[j,t]) for j in 1:nN))
        end
    end

    # Curtailment
    @constraint(m, [i in 1:nN, t in 1:hours],
        P_Curt[i,t] <= G_Solar[i,t] + G_Wind[i,t])

    # --------------------------------------------------
    # Line limits
    # --------------------------------------------------
    angmin = deg2rad.(dLine.angmin)
    angmax = deg2rad.(dLine.angmax)

    for k in 1:nL
        if dLine.status[k] != 0
            f  = dLine.fbus[k]
            to = dLine.tbus[k]
            @constraint(m, [t in 1:hours],
                angmin[k] <= θ[f,t] - θ[to,t] <= angmax[k])
            rateA = dLine.rateA[k] / bMVA
            Bij   = B[f,to]
            @constraint(m, [t in 1:hours],
                -rateA <= Bij*(θ[f,t] - θ[to,t]) <= rateA)
        end
    end

    # --------------------------------------------------
    # Slack (type == 3)
    # --------------------------------------------------
    for t in 1:hours
        for row in eachrow(dNodes[t])
            if row.type == 3
                @constraint(m, θ[row.bus_i, t] == 0)
            end
        end
    end

    # --------------------------------------------------
    # Costs:
    #   • Cm * P_g  → variable energy cost
    #   • CnL * y_g → “no load” cost
    #   • Cstart * y_sg → Start generating cost
    # --------------------------------------------------
    @objective(m, Min,
        sum(Cm[g]*P_g[g,t]*bMVA   for g in 1:nG, t in 1:hours) +
        sum(CnL[g]*y_g[g,t]       for g in 1:nG, t in 1:hours) +
        sum(Cstart[g]*y_sg[g,t]   for g in 1:nG, t in 1:hours))

    optimize!(m)

    # --------------------------------------------------
    # RESULTS
    # --------------------------------------------------
    solGen     = DataFrame(hour = Int[], gen = Int[], bus = Int[], powerGen = Float64[])
    solCurt    = DataFrame(hour = Int[], bus = Int[], P_curtailment = Float64[])
    solFlows   = DataFrame(hour = Int[], fbus = Int[], tbus = Int[], flow = Float64[])
    solVoltage = DataFrame(hour = Int[], bus = Int[], voltageNode = Float64[], angleDegrees = Float64[])
    solCosts   = DataFrame(hour = Int[], operation_cost = Float64[])

    if termination_status(m) in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ITERATION_LIMIT)
        for t in 1:hours
            for g in 1:nG
                push!(solGen, (
                    hour = t,
                    gen  = g,
                    bus  = gen_bus[g],
                    powerGen = value(P_g[g,t]) * bMVA))
            end
            for i in 1:nN
                push!(solCurt, (
                    hour = t,
                    bus  = i,
                    P_curtailment = value(P_Curt[i,t]) * bMVA))
            end
            for k in 1:nL
                f  = dLine.fbus[k]
                to = dLine.tbus[k]
                flow = B[f,to] * (value(θ[f,t]) - value(θ[to,t])) * bMVA
                push!(solFlows, (
                    hour = t,
                    fbus = f,
                    tbus = to,
                    flow = flow))
            end
            for i in 1:nN
                push!(solVoltage, (
                    hour = t,
                    bus  = i,
                    voltageNode = 1.0,
                    angleDegrees = rad2deg(value(θ[i,t]))))
            end
            cost_t = 0.0
            for g in 1:nG
                cost_t += Cm[g]*value(P_g[g,t])*bMVA +
                          CnL[g]*value(y_g[g,t]) +
                          Cstart[g]*value(y_sg[g,t])
            end
            push!(solCosts, (hour = t, operation_cost = cost_t))
        end
        push!(solCosts, (hour = -1, operation_cost = objective_value(m)))
    else
        println("No optimal solution found in DC_OPF_UC")
    end

    return m, solGen, solFlows, solVoltage, solCosts, solCurt
end

end 
