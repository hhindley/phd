using DataFrames, CSV, InteractiveViz, GLMakie, ProgressBars

include("/home/hollie_hindley/Documents/stochastic_hybrid/analysis_funcs.jl")

kdam_vals = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]

dfs = [DataFrame(CSV.File(joinpath("/home/hollie_hindley/Documents/stochastic_hybrid/kdam_test1", file), header=["event", "time", "rm_a", "rtca", "rm_b", "rtcb", "rm_r", "rtcr", "rh", "rd", "rt", "volume", "totprop"])) for file in readdir("/home/hollie_hindley/Documents/stochastic_hybrid/kdam_test1")]

props, df_rs, df_ps = df_sort_all(dfs)

f = Figure()
ilines(f[1,1], range(0,df_ps[1].time[end], length=length(df_ps[1].rtca)), df_ps[1][:,:rh]./df_ps[1].volume);
ilines!(f[1,1], range(0,df_ps[1].time[end], length=length(df_ps[1].rtca)), df_ps[1][:,:rt]./df_ps[1].volume);
ilines!(f[1,1], range(0,df_ps[1].time[end], length=length(df_ps[1].rtca)), df_ps[1][:,:rd]./df_ps[1].volume);

ilines(f[2,1], range(0,df_ps[2].time[end], length=length(df_ps[2].rtca)), df_ps[2][:,:rh]./df_ps[2].volume);
ilines!(f[2,1], range(0,df_ps[2].time[end], length=length(df_ps[2].rtca)), df_ps[2][:,:rt]./df_ps[2].volume);
ilines!(f[2,1], range(0,df_ps[2].time[end], length=length(df_ps[2].rtca)), df_ps[2][:,:rd]./df_ps[2].volume);

ilines(f[3,1], range(0,df_ps[3].time[end], length=length(df_ps[3].rtca)), df_ps[3][:,:rh]./df_ps[3].volume);
ilines!(f[3,1], range(0,df_ps[3].time[end], length=length(df_ps[3].rtca)), df_ps[3][:,:rt]./df_ps[3].volume);
ilines!(f[3,1], range(0,df_ps[3].time[end], length=length(df_ps[3].rtca)), df_ps[3][:,:rd]./df_ps[3].volume);

ilines(f[4,1], range(0,df_ps[4].time[end], length=length(df_ps[4].rtca)), df_ps[4][:,:rh]./df_ps[4].volume);
ilines!(f[4,1], range(0,df_ps[4].time[end], length=length(df_ps[4].rtca)), df_ps[4][:,:rt]./df_ps[4].volume);
ilines!(f[4,1], range(0,df_ps[4].time[end], length=length(df_ps[4].rtca)), df_ps[4][:,:rd]./df_ps[4].volume);

ilines(f[5,1], range(0,df_ps[5].time[end], length=length(df_ps[5].rtca)), df_ps[5][:,:rh]./df_ps[5].volume);
ilines!(f[5,1], range(0,df_ps[5].time[end], length=length(df_ps[5].rtca)), df_ps[5][:,:rt]./df_ps[5].volume);
ilines!(f[5,1], range(0,df_ps[5].time[end], length=length(df_ps[5].rtca)), df_ps[5][:,:rd]./df_ps[5].volume);

ilines(f[6,1], range(0,df_ps[6].time[end], length=length(df_ps[6].rtca)), df_ps[6][:,:rh]./df_ps[6].volume);
ilines!(f[6,1], range(0,df_ps[6].time[end], length=length(df_ps[6].rtca)), df_ps[6][:,:rt]./df_ps[6].volume);
ilines!(f[6,1], range(0,df_ps[6].time[end], length=length(df_ps[6].rtca)), df_ps[6][:,:rd]./df_ps[6].volume);

ilines(f[7,1], range(0,df_ps[7].time[end], length=length(df_ps[7].rtca)), df_ps[7][:,:rh]./df_ps[7].volume);
ilines!(f[7,1], range(0,df_ps[7].time[end], length=length(df_ps[7].rtca)), df_ps[7][:,:rt]./df_ps[7].volume);
ilines!(f[7,1], range(0,df_ps[7].time[end], length=length(df_ps[7].rtca)), df_ps[7][:,:rd]./df_ps[7].volume);

ilines(f[8,1], range(0,df_ps[8].time[end], length=length(df_ps[8].rtca)), df_ps[8][:,:rh]./df_ps[8].volume);
ilines!(f[8,1], range(0,df_ps[8].time[end], length=length(df_ps[8].rtca)), df_ps[8][:,:rt]./df_ps[8].volume);
ilines!(f[8,1], range(0,df_ps[8].time[end], length=length(df_ps[8].rtca)), df_ps[8][:,:rd]./df_ps[8].volume);

Label(f[1,1,Top()], "rh, rt, rd")









f = Figure()
ilines(f[1,1], range(0,df_ps[3].time[end], length=length(df_ps[3].rtca)), df_ps[3][:,:rh]./df_ps[3].volume);
ilines!(f[1,1], range(0,df_ps[3].time[end], length=length(df_ps[3].rtca)), df_ps[3][:,:rt]./df_ps[3].volume);
ilines!(f[1,1], range(0,df_ps[3].time[end], length=length(df_ps[3].rtca)), df_ps[3][:,:rd]./df_ps[3].volume);

ilines(f[2,1], range(0,df_ps[4].time[end], length=length(df_ps[4].rtca)), df_ps[4][:,:rh]./df_ps[4].volume);
ilines!(f[2,1], range(0,df_ps[4].time[end], length=length(df_ps[4].rtca)), df_ps[4][:,:rt]./df_ps[4].volume);
ilines!(f[2,1], range(0,df_ps[4].time[end], length=length(df_ps[4].rtca)), df_ps[4][:,:rd]./df_ps[4].volume);

ilines(f[3,1], range(0,df_ps[6].time[end], length=length(df_ps[6].rtca)), df_ps[6][:,:rh]./df_ps[6].volume);
ilines!(f[3,1], range(0,df_ps[6].time[end], length=length(df_ps[6].rtca)), df_ps[6][:,:rt]./df_ps[6].volume);
ilines!(f[3,1], range(0,df_ps[6].time[end], length=length(df_ps[6].rtca)), df_ps[6][:,:rd]./df_ps[6].volume);

ilines(f[4,1], range(0,df_ps[8].time[end], length=length(df_ps[8].rtca)), df_ps[8][:,:rh]./df_ps[8].volume);
ilines!(f[4,1], range(0,df_ps[8].time[end], length=length(df_ps[8].rtca)), df_ps[8][:,:rt]./df_ps[8].volume);
ilines!(f[4,1], range(0,df_ps[8].time[end], length=length(df_ps[8].rtca)), df_ps[8][:,:rd]./df_ps[8].volume);

ilines(f[5,1], range(0,df_ps[9].time[end], length=length(df_ps[9].rtca)), df_ps[9][:,:rh]./df_ps[9].volume);
ilines!(f[5,1], range(0,df_ps[9].time[end], length=length(df_ps[9].rtca)), df_ps[9][:,:rt]./df_ps[9].volume);
ilines!(f[5,1], range(0,df_ps[9].time[end], length=length(df_ps[9].rtca)), df_ps[9][:,:rd]./df_ps[9].volume);

ilines(f[6,1], range(0,df_ps[11].time[end], length=length(df_ps[11].rtca)), df_ps[11][:,:rh]./df_ps[11].volume);
ilines!(f[6,1], range(0,df_ps[11].time[end], length=length(df_ps[11].rtca)), df_ps[11][:,:rt]./df_ps[11].volume);
ilines!(f[6,1], range(0,df_ps[11].time[end], length=length(df_ps[11].rtca)), df_ps[11][:,:rd]./df_ps[11].volume);

ilines(f[7,1], range(0,df_ps[13].time[end], length=length(df_ps[13].rtca)), df_ps[13][:,:rh]./df_ps[13].volume);
ilines!(f[7,1], range(0,df_ps[13].time[end], length=length(df_ps[13].rtca)), df_ps[13][:,:rt]./df_ps[13].volume);
ilines!(f[7,1], range(0,df_ps[13].time[end], length=length(df_ps[13].rtca)), df_ps[13][:,:rd]./df_ps[13].volume);

ilines(f[8,1], range(0,df_ps[14].time[end], length=length(df_ps[14].rtca)), df_ps[14][:,:rh]./df_ps[14].volume);
ilines!(f[8,1], range(0,df_ps[14].time[end], length=length(df_ps[14].rtca)), df_ps[14][:,:rt]./df_ps[14].volume);
ilines!(f[8,1], range(0,df_ps[14].time[end], length=length(df_ps[14].rtca)), df_ps[14][:,:rd]./df_ps[14].volume);

Label(f[1,1,Top()], "rh, rt, rd")






f = Figure()
ilines(f[1,1], range(0,df_ps[3].time[end], length=length(df_ps[3].rtca)), df_ps[3][:,x]./df_ps[3].volume);

ilines(f[2,1], range(0,df_ps[4].time[end], length=length(df_ps[4].rtca)), df_ps[4][:,x]./df_ps[4].volume);

ilines(f[3,1], range(0,df_ps[6].time[end], length=length(df_ps[6].rtca)), df_ps[6][:,x]./df_ps[6].volume);

ilines(f[4,1], range(0,df_ps[8].time[end], length=length(df_ps[8].rtca)), df_ps[8][:,x]./df_ps[8].volume);

ilines(f[5,1], range(0,df_ps[9].time[end], length=length(df_ps[9].rtca)), df_ps[9][:,x]./df_ps[9].volume);

ilines(f[6,1], range(0,df_ps[11].time[end], length=length(df_ps[11].rtca)), df_ps[11][:,x]./df_ps[11].volume);

ilines(f[7,1], range(0,df_ps[13].time[end], length=length(df_ps[13].rtca)), df_ps[13][:,x]./df_ps[13].volume);

ilines(f[8,1], range(0,df_ps[14].time[end], length=length(df_ps[14].rtca)), df_ps[14][:,x]./df_ps[14].volume);

Label(f[1,1,Top()], "$x")
