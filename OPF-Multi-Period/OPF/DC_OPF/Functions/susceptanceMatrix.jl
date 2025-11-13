# This function creates the susceptance matrix B using the line data

function susceptanceMatrix(data::DataFrame, nn::Int, nl::Int)

    # data:    DataFrame containing line data
    # nn:      Number of nodes
    # nl:      Number of lines

    # Using the data of the endpoints of each line (fbus and tbus), the incidence matrix is created,
    # where we assign 1 to the fbus nodes and -1 to the tbus nodes.
    # For more information, refer to: https://en.wikipedia.org/wiki/Incidence_matrix
    # For the SparseArrays' sparse function, the arguments are:
    # sparse([Row Indices], [Column Indices], [Value], [Total Number of Rows], [Total Number of Columns])
    A = SparseArrays.sparse(data.fbus, 1:nl, 1, nn, nl) + SparseArrays.sparse(data.tbus, 1:nl, -1, nn, nl)

    # A vector is created with the susceptance values of each line: B = -1/x
    B = -1 ./ (data.x) .* data.status

    # Once we have the Incidence Matrix "A" and the Susceptance Vector "B",
    # the Susceptance Matrix "B_0" can be created:
    B_0 = A * SparseArrays.spdiagm(B) * A'
    # Here, spdiagm creates a matrix without elements (SparseArray) and assigns the elements of vector B to the main diagonal.
    
    # The susceptance matrix is returned
    return B_0

end