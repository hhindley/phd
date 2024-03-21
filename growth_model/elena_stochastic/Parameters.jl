#### Author: Elena Pascual Garcia 
#### Date: 14th Dec 2023
#### Project: Hybrid Stochastic Bacterial Growth & Division Model + Antibiotics


### Parameter values: Define the value for each parameter name through the 
### indexing structure defined in Index.jl. Input parameter is "J", the indexing 
### structure for the parameters. Output = vector with all parameter values. 
### Usage example can be found at the end of the script


function get_parameters(J)
    par = fill(NaN, 1, J.nrOfItems)

    SF = 22000
    # Nutrient import and metabolism
    par[pidx("NUTR_ext")] = 1.0e6           # Amount of substrate [#]
    par[pidx("ns")]       = 0.5             # Nutrient quality [AU]
    par[pidx("vt")]       = 660             # max. nutrient import rate [1/min]
    par[pidx("vm")]       = 5800.0          # max. enzymatic rate [1/min]
    par[pidx("Km")]       = 1.0e3           # enzymatic threshold [1/(10^8 aa)]
    par[pidx("Kt")]       = 1.0e3           # nutrient import threshold [1/(10^8 aa)]

    # Transcription
    par[pidx("thetar")]   = 426.8693338968694 * SF      # ribosome transcription threshold [1/(10^8 aa)]
    par[pidx("thetax")]   = 4.379733394834643 * SF      # non-ribosome transcription threshold [1/(10^8 aa)]
    par[pidx("we")]       = 0.003 * 4.139172187824451   # max. enzyme transcription rate [1/(aa min)]
    par[pidx("wr")]       = 0.008 * 929.9678874564831   # max ribosome transcription rate [1/(aa min)]
    par[pidx("wq")]       = 0.02 * 948.9349882947897    # max HKP transcription rate [1/(aa min)]

    # Autoinhibition for housekeeping genes
    par[pidx("Kq")]       = 1.522190403737490e5         # q-autoinhibition threshold [1/(10^8 aa)]
    par[pidx("nq")]       = 4                           # q-autoinhibition Hill. coeff [AU]

    # Translation
    par[pidx("gmax")]     = 12600.0                     # max. transl. elongation rate [aa/min]
    par[pidx("Kp")]       = 1801.378030928276
    par[pidx("Kp")]       = (par[pidx("gmax")] / par[pidx("Kp")]) * SF  # transl. elongation threshold [1/(10^8 aa)]
    par[pidx("nx")]       = 300.0                       # length of non-ribosomal proteins [aa]
    par[pidx("nr")]       = 7549.0                      # ribosome length [aa]

    par[pidx("kb")]       = 1.0                         # mRNA-ribosome binding rate [aa/min]
    par[pidx("ku")]       = 1.0                         # mRNA-ribosome unbinding rate [1/min]

    par[pidx("dm")]       = 0.1                         # mRNA-degradation rate [1/min]

    par[pidx("Mref")]     = 1.0e8                       # Typical cell total mass [aa]

    # Chloramphenicol effect
    par[pidx("AB_ext")]   = 0                           # External antibiotic concentration [uM]
    par[pidx("kon")]      = 19.3                        # Antibiotic binding rate 
    par[pidx("KD")]       = 0.03                        # Antibiotic unbinding/binding ratio
    par[pidx("pin")]      = 0.09                        # Antibiotic import rate
    par[pidx("pout")]     = 0.02                        # Antibiotic output rate
    
    par[pidx("NA")]       = 602214085700000000519317.96307935    # Avogadro's number [molec/mol]

    return par
end



# ##### ============ Usage =================
# par = get_parameters(J)
# par[pidx("C")]