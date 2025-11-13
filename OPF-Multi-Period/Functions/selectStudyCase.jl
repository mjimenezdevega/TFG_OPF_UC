function selectStudyCase()
    println("*** selectStudyCase v3 (UC-Copper, DC-OPF, AC-OPF) ***")
    println(">>> cargado desde: ", @__FILE__)  

    while true
        caseList = readdir("Cases")
        selectedCase = chooseOption(caseList, "case")

        opfTypeList = ["UC-Copper", "DC-OPF", "AC-OPF"]
        opfType = chooseOption(opfTypeList, "type of OPF")

        if opfType == "UC-Copper"
            solversList = ["HiGHS", "Gurobi"]
        elseif opfType == "DC-OPF"
            solversList = ["HiGHS", "Gurobi", "Ipopt"]
        else
            solversList = ["Ipopt", "Couenne"]
        end
        solver = chooseOption(solversList, "solver")

        println("Summary:")
        println("Study case _______ : ", selectedCase)
        println("Type of OPF ______ : ", opfType)
        println("Optimizer ________ : ", solver)
        println("\nPress ENTER to continue or any other input to reselect.")
        response = readline()
        if response == ""
            return selectedCase, opfType, solver
        end
    end
end
