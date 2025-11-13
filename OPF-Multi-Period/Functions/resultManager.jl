# This function manages the model variable and the DataFrames of the Optimization solution

function resultManager(model, genSolution, flowSolution, voltageSolution, solCosts, solCurt, solStorage, mFilePath, opfType, solver)

    # model: The model created for optimization
    # genSolution: DataFrame containing the generators' solution
    # flowSolution: DataFrame containing the flow solution
    # voltageSolution: DataFrame containing the voltage solution (magnitude and angle)
    # solCosts: DataFrame containing solutions costs per hour and total sum of it
    # solCurt: Dataframe containing renewable curtailment per hour

    # Clear the terminal
    #clearTerminal()

    if model === nothing
        println("ERROR: El modelo es 'nothing', no se ha generado correctamente.")
        return
    end

    # Display results if the optimization was successful (globally, locally, or iteration limit reached)
    if termination_status(model) == OPTIMAL || termination_status(model) == LOCALLY_SOLVED || termination_status(model) == ITERATION_LIMIT

        # In case of a global solution
        if termination_status(model) == OPTIMAL
            println("Optimal solution found")

        # In case of a local solution
        elseif termination_status(model) == LOCALLY_SOLVED
            println("Local solution found")

        # In case of reaching the iteration limit
        elseif termination_status(model) == ITERATION_LIMIT
            println("Iteration limit reached")
        end

        # ##### VS Code terminal does not display the plot, but the Julia terminal does
        solution = 0
        if mFilePath != "None"
            # If an LP_OPF is to be solved
            if opfType == "LP-OPF"
                if solver == "Gurobi"
                    solution = solve_opf(mFilePath, DCMPPowerModel, Gurobi.Optimizer)
                elseif solver == "HiGHS"
                    solution = solve_opf(mFilePath, DCMPPowerModel, HiGHS.Optimizer)
                elseif solver == "Ipopt"
                    solution = solve_opf(mFilePath, DCMPPowerModel, Ipopt.Optimizer)
                else
                    print("Error loading DC resolution with PowerModels")
                end
            elseif opfType == "AC-OPF"
                if solver == "Ipopt"
                    solution = solve_opf(mFilePath, ACRPowerModel, Ipopt.Optimizer)
                else
                    print("Error loading AC resolution with PowerModels")
                end
            else
                println("Error loading solver type in PowerModels")
            end

            clearTerminal()

        else
            println(".m file for the case not found\n")
        end

        # Check the number of rows in the solution DataFrames
        genRows = DataFrames.nrow(genSolution)
        flowRows = DataFrames.nrow(flowSolution)
        voltageRows = DataFrames.nrow(voltageSolution)

        # Set the maximum number of rows to display
        maxRows = 10

        if genRows <= maxRows && flowRows <= maxRows && voltageRows <= maxRows
            println("\nDo you want to display the results in the terminal?")
            println("Press ENTER to confirm or any other input to deny")
            displayTerminal = readline(stdin)

            clearTerminal()

            if displayTerminal == ""
                println("Generator solution:")
                DataFrames.show(genSolution, allrows = true, allcols = true)
                println("\n\nFlow solution:")
                DataFrames.show(flowSolution, allrows = true, allcols = true)
                println("\n\nAngle solution:")
                DataFrames.show(voltageSolution, allrows = true, allcols = true)
                println("\n\nOperation Costs:")
                DataFrames.show(solCosts, allrows = true, allcols = true)
                println("") 
                println("\n\nCurtailment Results:")
                DataFrames.show(solCurt, allrows = true, allcols = true)
                println("")
                println("\n\nStorage Results:")
                DataFrames.show(solStorage, allrows = true, allcols = true)
            else
                println("\nResults will not be displayed")
            end
        else
            println("Tables are too large to display in the terminal")
        end

        # Print the solution obtained with PowerModels.jl if applicable (for las hour)
        if solution != 0
            println("Final cost obtained with PowerModels: ", round(solution["objective"], digits = 2), "€/h")
            println("Program execution time: ", solution["solve_time"] * 1000, " ms")
        end

        println("\nFinal cost with the program: ", round(solCosts.operation_cost[end], digits = 2), " €/h")
        println("Program execution time: ", solve_time(model) * 1000, " ms")
        println("\nDo you want to save the results in a CSV file?")
        println("Press ENTER to confirm or any other input to deny")
        saveCSV = readline(stdin)

        if saveCSV == ""
            println("Saving in CSV will overwrite any existing data")
            println("Are you sure you want to save?")
            println("\nPress ENTER to confirm or any other input to deny")
            confirmSaveCSV = readline(stdin)

            if confirmSaveCSV == ""
                CSV.write("./Results/voltageSolution.csv", voltageSolution, delim = ";")
                CSV.write("./Results/lineFlowSolution.csv", flowSolution, delim = ";")
                CSV.write("./Results/generatorSolution.csv", genSolution, delim = ";")
                CSV.write("./Results/costsSolution.csv", solCosts, delim = ";")
                CSV.write("./Results/curtailmentSolution.csv", solCurt, delim = ";")
                CSV.write("./Results/storageSolution.csv", solStorage, delim = ";")
                println("\nThe results have been saved in ./Results")
            else
                println("\nResults will not be saved")
            end
        else
            println("\nResults will not be saved")
        end
    else
        println("ERROR: ", termination_status(model))
    end

end