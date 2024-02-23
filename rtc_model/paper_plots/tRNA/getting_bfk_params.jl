ds_range = 10 .^ range(log10(1e-5),log10(1),length=10)
dsmax_range = 10 .^ range(log10(1e-5),log10(1),length=10)
θ_range = 10 .^ range(log10(1e-5),log10(1),length=10)
dsmin_range = 10 .^ range(log10(1e-5),log10(1),length=10)
brs=[]
vals=[]
for ds in ds_range
    for dsmin in dsmin_range
        for dsmax in dsmax_range
            for θ in θ_range
                try
                    prob = ODEProblem(rtc_trna_model, ssvals_trna, tspan, params_trna; jac=true)
                    odefun = prob.f
                    F = (u,p) -> odefun(u,p,0)
                    J = (u,p) -> odefun.jac(u,p,0)
                    id_kdam = indexof(kdam, parameters(rtc_trna_model))
                    par_tm = prob.p
                    prob_bf = BifurcationProblem(F, prob.u0, setproperties(par_tm), (@lens _[id_kdam]); J = J,
                    record_from_solution = (x, p) -> (rm_a = x[1], rtca = x[2], rm_b = x[3], rtcb = x[4], rm_r = x[5], rtcr = x[6], trna = x[7], rd = x[8], rt = x[9]),)
                    opts_br = ContinuationPar(p_min = 0., p_max = 200., ds = ds,# a=a,
                    dsmax = dsmax, dsmin = dsmin,# 0.15
                    # options to detect bifurcations
                    detect_bifurcation = 3, n_inversion = 2, max_bisection_steps = 10, #3,2,10
                    # number of eigenvalues
                    # nev =100, #tolParamBisectionEvent=1e-30, 
                    # maximum number of continuation steps
                    max_steps = 1000,)# dsmin_bisection=1e-30)#, tol_bisection_eigenvalue=1e-10)# a=0.9, )
                    # tolStability=1e-10, tolBisectionEigenvalue=1e-10)#,tolParamBisectionEvent=1e-1)
                    # only using parameters that make a difference to solution
                    # continuation of equilibria
                    
                    br = continuation(prob_bf, PALC(θ=θ), opts_br; plot = false, bothside=true, normC = norminf)
                catch
                    println("skipping error")
                end
                if br.specialpoint[1].param == 0. && br.specialpoint[end].param == 1. && all(x->x==1, ([br.specialpoint[i].status == :converged for i in range(2,length(br.specialpoint)-1)])) == true && length([br.specialpoint[i].status == :converged for i in range(2,length(br.specialpoint)-1)]) > 1 && all(x->x>0, br.rm_a[end-5:end])
                 
                    push!(brs, br)
                    push!(vals, [ds, dsmin, dsmax, θ])
                end
            end
        end
    end
end
brs
vals

p=Plots.plot()
for i in range(1,28)
    Plots.plot!(p,brs[i])
end
p

vals[1:28]

df = DataFrame(ds=[val[1] for val in vals[1:28]], dsmin=[val[2] for val in vals[1:28]], dsmax=[val[3] for val in vals[1:28]], θ=[val[4] for val in vals[1:28]])
df1 = DataFrame(CSV.File("/home/holliehindley/phd/rtc_model/paper_plots/tRNA/bf_params.csv"))

