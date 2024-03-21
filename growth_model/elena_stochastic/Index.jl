#### Author: Elena Pascual Garcia 
#### Date: 14th Dec 2023
#### Project: Hybrid Stochastic Bacterial Growth & Division Model + Antibiotics


### Indexing structure: Define the state and parameter names and assign an 
### indexing number (position in vector) for each name. This will help with 
### consistency across the rest of the model definition. (see example of use
### at the end of the script)

struct MyIndexDict
    indices::Dict{Symbol, Int}
    nrOfItems::Int
    allItems::UnitRange{Int}
end

# Function to create an IndexDict with given names
function create_index_dict(names)
    indices = Dict{Symbol, Int}()
    for (i, name) in enumerate(names)
        indices[Symbol(name)] = i
    end
    return MyIndexDict(indices, length(names), 1:length(names))
end

# Helper function for compact indexing
function index(dict::MyIndexDict, item::String)
    sym_item = Symbol(item)
    haskey(dict.indices, sym_item) || throw(ArgumentError("Item $item not found in the index dictionary"))
    return dict.indices[sym_item]
end

function indexing()
    stateNames = ["NUTR", "ATP",
                  "mENZ", "mTRP", "mHKP", "mRIB",
                  "RIB_mENZ", "RIB_mTRP", "RIB_mHKP", "RIB_mRIB",
                  "AB_RIB_mENZ", "AB_RIB_mTRP", "AB_RIB_mHKP", "AB_RIB_mRIB",
                  "ENZ", "TRP", "HKP", "RIB", "AB", "TotProp"]

    I = create_index_dict(stateNames)

    parNames = ["NUTR_ext", "ns", "vt", "vm", "Km", "Kt", "thetar",
                "thetax", "we", "wr", "wq", "Kq", "nq", "gmax", "Kp", "nx",
                "nr", "kb", "ku", "dm", "Mref",
                "AB_ext", "kon", "KD", "pin", "pout",
                "NA"]

    J = create_index_dict(parNames)

    return (I, J)
end

    global (I, J) = indexing()

    # Compact indexing function
    sidx(item) = index(I, item)
    pidx(item) = index(J, item)



# ####### =======================================================================
# ####### ======================= Example of usage ==============================

#     # Create the indexing structure
#     (I, J) = indexing()

#     # Example vector of states
#     states = rand(20)

#     # Access the state corresponding to "NUTR"
#     nutr_state = states[sidx("NUTR")]

#     # Example vector of parameters
#     parameters = rand(30)

#     # Access the parameter corresponding to "NUTR_ext"
#     nutr_param = parameters[pidx("NUTR_ext")]

#     # Display the results
#     println("State value for NUTR: ", nutr_state)
#     println("Parameter value for NUTR_ext: ", nutr_param)

#     # Access the nrOfItems for states
#     nr_of_states = I.nrOfItems

#     # Access the nrOfItems for parameters
#     nr_of_params = J.nrOfItems

#     # Display the results
#     println("Number of states: ", nr_of_states)
#     println("Number of parameters: ", nr_of_params)