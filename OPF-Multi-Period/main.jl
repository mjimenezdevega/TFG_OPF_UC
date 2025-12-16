# --------------------------------
# main.jl
# --------------------------------

include("./Functions/loadLibraries.jl")
include("./Functions/loadFunctions.jl")

include(joinpath(@__DIR__, "Functions", "chooseOption.jl"))
include(joinpath(@__DIR__, "Functions", "clearTerminal.jl"))
include(joinpath(@__DIR__, "Functions", "runCaseUCCopper24h.jl"))
include(joinpath(@__DIR__, "Functions", "runCase3NudosUC.jl"))
include(joinpath(@__DIR__, "Functions", "runCase3NudosAmpliacion.jl"))
include(joinpath(@__DIR__, "Functions", "runCase19Nudos.jl"))
include(joinpath(@__DIR__, "Functions", "selectStudyCase.jl"))

function main()
    endProgram = false
    while !endProgram
        clearTerminal()

        selectedCase, opfType, s = selectStudyCase()
        clearTerminal()

        println("\nGenerating OPF...")

    

        # =======================
        # Solve based in type
        # =======================

        

        if (selectedCase == "3 Nudos UC") && (opfType == "DC-OPF")
            m, base, solGen, solFlows, solVoltage, solCosts, solCurt =
                run_case_3nudos_uc(; solver_name = s)
            println("\nProblem solved (DC-OPF - 3 Nudos UC).")

        elseif (selectedCase == "3 Nudos Ampliación") && (opfType == "DC-OPF")
            m, base, solGen, solFlows, solVoltage, solCosts, solCurt =
            run_case_3nudos_ampliacion(; solver_name = s)

            println("\nProblem solved (DC-OPF - 3 Nudos Ampliación).")
       
        elseif (selectedCase == "19 Nudos") && (opfType == "DC-OPF")
            m, solGen, solFlows, solVoltage, solCosts, solCurt =
                run_case_19nudos(; solver_name = s)
            println("\nProblem solved (DC-OPF-UC - 19 Nudos).")


        elseif opfType == "DC-OPF"
            m, solGen, solFlows, solVoltage, solCosts, solCurt, solStorage =
                DC_OPF(data[1], data[2], data[3], data[4], data[5], data[6], s)
            println("\nProblem solved (DC-OPF).")

        elseif opfType == "AC-OPF"
            m, solGen, solFlows, solVoltage =
                AC_OPF(data[1], data[2], data[3], data[4], data[5], data[6], s)
            println("\nProblem solved (AC-OPF).")

        elseif (selectedCase == "UC Copper 24h") && (opfType == "UC-Copper")
            model = run_case_24h(; solver_name = s)
            println("\n--- Finished UC-Copper 24h ---")
            return model

        else
            println("ERROR: Failed to load OPF type")
        end

        println("\nPress ENTER to continue or any other input to exit.")
        if readline() == ""
            endProgram = false
        else
            endProgram = true
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
