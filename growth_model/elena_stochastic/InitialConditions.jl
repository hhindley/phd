#### Author: Elena Pascual Garcia 
#### Date: 14th Dec 2023
#### Project: Hybrid Stochastic Bacterial Growth & Division Model + Antibiotics


### Initial conditions: Define the starting value for each species through the 
### indexing structure defined in Index.jl. Input parameter is "I", the indexing 
### structure for the state variables. Output = vector with all initial conditions. 
### Usage example can be found at the end of the script

function get_X0(I)
   
    X0 = zeros(I.nrOfItems)

    # Manually set initial conditions
    X0[sidx("NUTR")] = 0
    X0[sidx("ATP")] = 10
    X0[sidx("mENZ")] = 0
    X0[sidx("mTRP")] = 0
    X0[sidx("mHKP")] = 0
    X0[sidx("mRIB")] = 0
    X0[sidx("RIB_mENZ")] = 0
    X0[sidx("RIB_mTRP")] = 0
    X0[sidx("RIB_mHKP")] = 0
    X0[sidx("RIB_mRIB")] = 0
    X0[sidx("AB")] = 0
    X0[sidx("AB_RIB_mENZ")] = 0
    X0[sidx("AB_RIB_mTRP")] = 0
    X0[sidx("AB_RIB_mHKP")] = 0
    X0[sidx("AB_RIB_mRIB")] = 0
    X0[sidx("ENZ")] = 0
    X0[sidx("TRP")] = 0
    X0[sidx("HKP")] = 0
    X0[sidx("RIB")] = 10
    X0[sidx("TotProp")] = 0

    return X0
end

# # ##### ============ Usage =================
# X0 = get_X0(J)
# X0[sidx("RIB")]

