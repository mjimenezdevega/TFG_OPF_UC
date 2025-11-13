function extractData(c::String)
    # Line data
    lineData = CSV.read("Cases/$c/lineData.csv", DataFrame)

    # Generator data
    generatorData = CSV.read("Cases/$c/generatorData.csv", DataFrame)

    
    # Demand data
    demandPath = "Cases/$c/demandData/"

    if isdir(demandPath)
        #If the directory demandData exists, read all files nodeData_X.csv
        nodeFiles = filter(x -> endswith(x, ".csv"), readdir(demandPath))  # CSV files list
        nodeDataList = [CSV.read(joinpath(demandPath, file), DataFrame) for file in nodeFiles] # Here we create a list of DataFrames with each hourly file.
    else
        #If demanData doesnt exist, we read a single demand file
        nodeDataList = CSV.read("Cases/$c/nodeData.csv", DataFrame)
    end

    # Solar Gen data
    solarData = CSV.read("Cases/$c/solarData.csv", DataFrame)

    # Wind Gen data
    windData = CSV.read("Cases/$c/windData.csv", DataFrame)

    # Battery data
    storageData = CSV.read("Cases/$c/storage.csv", DataFrame)

    # Number of nodes
    nNodes = maximum([lineData.fbus; lineData.tbus])

    # Number of lines
    nLines = size(lineData, 1)

    # Base power
    baseMVA = 100

    # Path to the .m file
    mFilePath = "Cases/$c/$c.m"

    if isfile(mFilePath)
        path = mFilePath
    else
        path = "None"
    end

    #Analyzed hours
    hours = length(nodeDataList)

    # Return all generated DataFrames and variables
    return(lineData, generatorData, nodeDataList, nNodes, nLines, baseMVA, path, hours, solarData, windData, storageData)
end