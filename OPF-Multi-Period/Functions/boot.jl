# This function precompiles all solvers so that they are already loaded during case execution

function boot()

    clearTerminal()

    println("Starting tests...")
    # Extract data from the test system
    test_line = CSV.read("Functions/test_system/lineData.csv", DataFrame)
    test_generator = CSV.read("Functions/test_system/generatorData.csv", DataFrame)
    test_nodes = CSV.read("Functions/test_system/nodeData.csv", DataFrame)
    # With this simple network, the different OPFs are generated so they are already loaded when the user uses them
    
    println("Test 1...")
    LP_OPF(test_line, test_generator, test_nodes, 2, 1, 100, "Gurobi",1)

    clearTerminal()

    println("Test 1 - Completed")
    println("Test 2...")
    LP_OPF(test_line, test_generator, test_nodes, 2, 1, 100, "HiGHS",1)

    clearTerminal()

    println("Test 1 - Completed")
    println("Test 2 - Completed")
    println("Test 3...")
    LP_OPF(test_line, test_generator, test_nodes, 2, 1, 100, "Ipopt",1)
    mPath = "Functions/test_system/test_system.m"
    solve_opf(mPath, DCMPPowerModel, Ipopt.Optimizer)

    clearTerminal()

    println("Test 1 - Completed")
    println("Test 2 - Completed")
    println("Test 3 - Completed")
    println("Test 4...")
    AC_OPF(test_line, test_generator, test_nodes, 2, 1, 100, "Ipopt")

    clearTerminal()

    println("Test 1 - Completed")
    println("Test 2 - Completed")
    println("Test 3 - Completed")
    println("Test 4 - Completed")
    println("Test 5...")
    AC_OPF(test_line, test_generator, test_nodes, 2, 1, 100, "Couenne")

    clearTerminal()
    
    println("Test 1 - Completed")
    println("Test 2 - Completed")
    println("Test 3 - Completed")
    println("Test 4 - Completed")
    println("Test 5 - Completed")
    sleep(1)
    
end
