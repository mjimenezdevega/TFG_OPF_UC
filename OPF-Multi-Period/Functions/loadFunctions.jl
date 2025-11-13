# Load all created functions

include("./boot.jl")                    # Initializer to load all functions

include("./clearTerminal.jl")           # Terminal cleaning
include("./chooseOption.jl")            # User choice and confirmation of the option
include("./extractData.jl")             # Extracts data from the selected case and stores it in SparseArrays
include("./selectStudyCase.jl")         # Choice of the case to study
include("./resultManager.jl")           # Manages the optimization result

include("../OPF/DC_OPF/DC_OPF.jl")      # Linear Programming - Optimal Power Flow function
include("../OPF/AC_OPF/AC_OPF.jl")      # Alternating Current - Optimal Power Flow function