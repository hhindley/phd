using LabelledArrays, DataFrames, CSV

conv = 1e9
avo = 6.022e23
da = 650

copies = dna*avo/(rtca_length*conv*da)

function calc_conc_copyno(rtc_length, copies)
    dna = rtc_length*conv*da*copies/(avo)
    return dna
end

function calc_copyno(rtc_length, dna_conc)
    return dna_conc*avo/(rtc_length*conv*da)
end

copies = [1, 10, 100, 1000, 10000]

dna_concs = [0.0001, 1.038e-5, 1.038e-6, 1.038e-7, 1.038e-8] 
rtc_lengths = @LArray [7025, 7235, 8222] (:rtca, :rtcb, :rtcr)

df_res = DataFrame(:RtcA=>[], :RtcB=>[], :RtcR=>[])

for (rtc, res) in zip(rtc_lengths, eachcol(df_res))
    for copy in dna_concs 
        println(copy)
        dna = calc_copyno(rtc,copy)
        push!(res, dna)
    end   
end

insertcols!(df_res, 1, :Copies=>copies)
df_res
CSV.write("/home/holliehindley/phd/data/plasmid_concs.csv", df_res)
