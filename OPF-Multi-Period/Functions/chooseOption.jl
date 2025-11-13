# Loop that returns the option the user selects once confirmed

function chooseOption(o::Vector{String}, type::String)

    # Initialize variables
    valid = false  # This variable checks if the choice is valid; it becomes true if so
    selection = 0  # The option selected by the user

    # The loop continues until the response is valid
    while !valid 

        # Enter a try-catch block to handle inputs that cause exceptions
        try

            # Clear the terminal
            clearTerminal()

            # Print the possible enumerated options in the terminal
            for (i, k) in enumerate(o)
                println("$i. $k")
            end

            # Ask the user to enter their choice in the terminal
            println("\nChoose the number of the ", type, " you want to use: ")
            selection = parse(Int, readline())

            # If the input is a number and within the range of possible options
            if selection >= 1 && selection <= length(o)

                # Clear the terminal
                clearTerminal()

                # Show the selected option in the terminal
                println("You have selected:\n", selection, ". ", o[selection])

                # Ask for confirmation
                println("\nPress ENTER to confirm or any other input to reselect.")
                confirm = readline()
                
                # If the input is the ENTER key
                if confirm == ""
                    # Update "valid" to exit the loop
                    valid = true

                # Otherwise
                else
                    # Ignore the input and restart the loop
                    continue

                end

            # If the entered number is out of range
            else

                # Clear the terminal
                clearTerminal()

                # Show a message in the terminal indicating the range
                println("Please enter a number between 1 and $(length(o)).")

                # Display the message for 2 seconds
                sleep(1)
                continue

            end

        # If the input causes an exception, 
        # for example, entering a letter that cannot be converted to an int
        catch

            # Clear the terminal
            clearTerminal()

            # Display the message in the terminal for 2 seconds
            println("Invalid input. Please enter a number.")
            sleep(1)
            continue

        end

    end

    # Return the option selected by the user
    return o[selection]
    
end