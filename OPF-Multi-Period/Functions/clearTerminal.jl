# This function is to clear the terminal where the code is running
function clearTerminal()

    # If the terminal is Windows
    if Sys.iswindows()
        Base.run(`cmd /c cls`)

    # For other Unix-based terminals
    else
        Base.run(`clear`)
    end

end