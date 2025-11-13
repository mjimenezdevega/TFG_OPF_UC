# --------------------------------
# main.jl
# --------------------------------

# Cargas (si ya las usas en tu repo)
include("./Functions/loadLibraries.jl")
include("./Functions/loadFunctions.jl")

# Utilidades necesarias para el menú y UC-Copper
include(joinpath(@__DIR__, "Functions", "chooseOption.jl"))
include(joinpath(@__DIR__, "Functions", "clearTerminal.jl"))
include(joinpath(@__DIR__, "Functions", "runCase24h.jl"))
include(joinpath(@__DIR__, "Functions", "selectStudyCase.jl"))

# Si tu flujo DC/AC necesita includes específicos, añádelos aquí:
# include(joinpath(@__DIR__, "OPF", "DC_OPF.jl"))
# include(joinpath(@__DIR__, "OPF", "AC_OPF.jl"))
# include(joinpath(@__DIR__, "Functions", "resultManager.jl"))

function main()
    # boot()  # descomenta si lo usas

    endProgram = false
    while !endProgram
        clearTerminal()

        # Menú (3 tipos: UC-Copper / DC-OPF / AC-OPF)
        selectedCase, opfType, s = selectStudyCase()
        clearTerminal()

        println("\nGenerating OPF...")

        # Solo DC/AC necesitan extraer datos del caso
        if opfType == "DC-OPF" || opfType == "AC-OPF"
            println("\nExtracting data...")
            data = extractData(selectedCase)
            println("Data extracted.")
            # hours = data[8]  # si lo usas
        end

        if opfType == "DC-OPF"
            # Ajusta la firma a tu DC_OPF real
            m, solGen, solFlows, solVoltage, solCosts, solCurt, solStorage =DC_OPF(data[1], data[2], data[3], data[4], data[5], data[6], s)
            println("\nProblem solved (DC-OPF).")

            # Opcional:
            # resultManager(m, solGen, solFlows, solVoltage, solCosts, solCurt, solStorage, data[7], opfType, s)

        elseif opfType == "AC-OPF"
            # Ajusta la firma a tu AC_OPF real
            m, solGen, solFlows, solVoltage =
                AC_OPF(data[1], data[2], data[3], data[4], data[5], data[6], s)
            println("\nProblem solved (AC-OPF).")

            # Opcional:
            # resultManager(m, solGen, solFlows, solVoltage, nothing, nothing, nothing, data[7], opfType, s)

        elseif (selectedCase == "Caso 24h") && (opfType == "UC-Copper")
            model = run_case_24h(; solver_name = s)        # "HiGHS" o "Gurobi"
            println("\n--- Finished UC-Copper 24h ---")
            return model                                   # salimos (no usar resultManager aquí)

        else
            println("ERROR: Failed to load OPF type")
        end

        # ÚNICO bloque de continuar/salir (dentro del while)
        println("\nPress ENTER to continue or any other input to exit.")
        if readline() == ""
            endProgram = false
        else
            endProgram = true
        end
    end # while
end # function

# Ejecuta main si llamas este archivo como script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
