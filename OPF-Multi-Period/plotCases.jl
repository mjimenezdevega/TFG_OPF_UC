using Plots
using DataFrames

# ==========================================================
# Construir mezcla por hora
# ==========================================================
function build_hourly_mix(df::DataFrame)
    groups = groupby(df, :hora)

    mix = DataFrame(
        hora    = Int[],
        PD_MW   = Float64[],
        RES_MW  = Float64[],
        CC_MW   = Float64[],
        CR_MW   = Float64[],
        Biom_MW = Float64[],
        Nuc_MW  = Float64[],
    )

    for sub in groups
        h = first(sub.hora)
        PD = first(sub.PD_MW)
        PRES_av = first(sub.PRES_avail_MW)
        Pcurt   = first(sub.Pcurt_MW)
        RES_eff = PRES_av - Pcurt

        Pg_CC   = sum(sub.Pg_MW[sub.tech .== "CC"])
        Pg_CR   = sum(sub.Pg_MW[sub.tech .== "CR"])
        Pg_Biom = sum(sub.Pg_MW[sub.tech .== "Biomasa"])
        Pg_Nuc  = sum(sub.Pg_MW[sub.tech .== "Nuclear"])

        push!(mix, (h, PD, RES_eff, Pg_CC, Pg_CR, Pg_Biom, Pg_Nuc))
    end

    sort!(mix, :hora)
    return mix
end


# ==========================================================
# Plot principal del mix
# ==========================================================
function plot_mix(df::DataFrame)
    mix = build_hourly_mix(df)

    bar(
        mix.hora,
        [mix.RES_MW mix.CC_MW mix.Biom_MW mix.Nuc_MW mix.CR_MW],
        bar_position = :stack,
        label = ["Renovables" "Gas (CC)" "Biomasa" "Nuclear" "CR"],
        color = [:green :gray :brown :purple :black],
        xlabel = "Hora",
        ylabel = "Potencia [MW]",
        title  = "Curva de demanda y mix de generación",
        legend = :topright,
    )

    plot!(
        mix.hora,
        mix.PD_MW,
        label = "Demanda",
        linewidth = 3,
        color = :red,
    )
end
