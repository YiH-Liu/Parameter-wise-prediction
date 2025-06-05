using Plots, DifferentialEquations
using CSV
using DataFrames
using LaTeXStrings
using Interpolations
using Distributions
using Measures
using NLopt
using FilePathsBase
using SpecialFunctions
using Dierckx


# In this script, we generate the result for Case 1

# Get the path of the current script.
path = dirname(@__FILE__)

# Read the CSV file back into a DataFrame
df = CSV.read("$path\\data.csv", DataFrame)

# Convert the DataFrame back into a matrix
step = 18
data1 = df.a[1:step:end]
data2 = df.b[1:step:end]
cmax=100
# Fixed parameters
# Data is collected at time t = t1.
t = 100
I = 200; J = 20; tfs=8
# union of prediction interval
ub_union1= zeros(I,4); ub_union2= zeros(I,4) 
lb_union1= zeros(I,4); lb_union2= zeros(I,4)


# solution of PDE
function solution(x,t,P,ρ,k1,k2,c10,c20,h,L) 
    D = P/4
    v = (P*ρ)/2
    x0 = L/2
    Δ = k1^2 + k2^2 + 6*k1*k2
    b1_a = c10*(-k1/(2*sqrt(Δ))) + c20*((-k1+k2+sqrt(Δ))/(4*sqrt(Δ)))
    b2_a = c10*(k1/(2*sqrt(Δ))) + c20*((k1-k2+sqrt(Δ))/(4*sqrt(Δ)))
    b1_b = exp(((-k1-k2-sqrt(Δ))/2)*t)
    b2_b = exp(((-k1-k2+sqrt(Δ))/2)*t)
    b1 = b1_a.*(erf.((h.-(x.-v*t.-x0))./(2*sqrt(D*t))) + erf.((h.+(x.-v*t.-x0))./(2*sqrt(D*t)))).*b1_b
    b2 = b2_a.*(erf.((h.-(x.-v*t.-x0))./(2*sqrt(D*t))) + erf.((h.+(x.-v*t.-x0))./(2*sqrt(D*t)))).*b2_b
    c1a = (-k1+k2-sqrt(Δ))/(2*k1); c1b = (-k1+k2+sqrt(Δ))/(2*k1)
    c1 = c1a.*b1 .+ c1b.*b2
    c2 = b1 .+ b2 
    return c1,c2
end 

# model(): The continuum model
function model(t,a) 
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # t: (number) Analytic solution of PDE at t.
    # a: (vector) Parameter vector.
    #------------------------------------------------------------------------------------------------------------------
    # Evaluate the analytic solution
    x = 0:1:199
    c1,c2 = solution(x,t,a[1],a[2],a[3],a[4],0.5,0.2,20,199)
    return c1,c2
end

# error(): The loglikelihood function for multinomial measurement error model for Case 1
function error(data1,data2,t,a)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # data1: (vector) Count data for subpopulation 1.
    # data2: (vector) Count data for subpopulation 2.
    # t: (number) Data are simulated at time t.
    # a: (vector) Parameter vector.
    # Output:
    # e: (number) Log-likelihood for the Additive Gaussian Measurement Error Model (Case 1) with parameter vector a.
    #------------------------------------------------------------------------------------------------------------------
    # Sample from continuum model at time t with parameter vector a.
    c1,c2 =model(t,a)
    c1 = c1[1:step:end]
    c2 = c2[1:step:end]
    # Estimate the loglikelihood function.
    e=0;
    I = length(data1)
    for i = 1:I
        e = e + loglikelihood(Poisson(c1[i]*J),data1[i]) + loglikelihood(Poisson(c2[i]*J),data2[i])
    end
    return e
end

# fun(): Evaluate the log-likelihood function as a function of the parameters a.
function fun(a)
    # ------------------------------------------------------------------------------------------------------------------
    # Input:
    # a: (vector) [D,v,k1,k2].
    # Output:
    # error(data1,data2,t,a): (number) Log-likelihood for the Additive Gaussian Measurement Error Model (Case 1) with parameter vector a.
    #------------------------------------------------------------------------------------------------------------------
    return error(data1,data2,t,a)
end

# optimise(): NLopt routine to maximise the function fun, with parameter estimates θ₀ subject to bound constraints lb, ub
function optimise(fun,θ₀,lb,ub;
    dv = false,
    method = dv ? :LD_LBFGS : :LN_BOBYQA,
)

    if dv || String(method)[2] == 'D'
        tomax = fun
        
    else
        tomax = (θ,∂θ) -> fun(θ)
        
    end
    
    opt = Opt(method,length(θ₀))
    opt.max_objective = tomax
    opt.lower_bounds = lb       # Lower bound
    opt.upper_bounds = ub       # Upper bound
    opt.local_optimizer = Opt(:LN_NELDERMEAD, length(θ₀))
  #=   opt.local_optimizer = Opt(:GN_CRS2_LM, length(θ₀)) =#
    opt.maxtime = 10*60
    res = optimize(opt,θ₀)
    return res[[2,1]]
end


#Expected Parameters
P = 1; ρ = 0.5; k1 = 0.02; k2 = 0.03;
θ = [P,ρ,k1,k2]
l_theta = fun(θ)
#MLE----------------------------------------------------------------------------
# Inital guess.
Pg = 1; ρg = 0.5; k1g = 0.025; k2g = 0.035
lb_P = 0.001   ; lb_ρ = -1 ; lb_k1 = 0.001; lb_k2 = 0.001
ub_P = 2; ub_ρ = 1; ub_k1 = 0.06; ub_k2 = 0.06
θG = [Pg,ρg,k1g,k2g] # inital guess
lb=[lb_P,lb_ρ,lb_k1,lb_k2] # lower bound
ub=[ub_P,ub_ρ,ub_k1,ub_k2] # upper bound
# Call numerical optimization routine to give the vector of parameters xopt, and the maximum loglikelihood fopt.
@time (xopt,fopt)  = optimise(fun,θG,lb,ub)
fmle=fopt
# Print MLE parameters
Pmle=xopt[1]; 
ρmle=xopt[2]; 
k1mle=xopt[3];
k2mle=xopt[4] 

println("Pmle: ", Pmle)
println("ρmle: ", ρmle)
println("k1mle: ", k1mle)
println("k2mle: ", k2mle)
println("likelihood at θ: ",l_theta)
println("likelihood at MLE: ",fmle)
c1mle,c2mle= model(t,[Pmle,ρmle,k1mle,k2mle])
c1,c2= model(t,θ)
scatter(data2)
#Profile P-------------------------------------------------------------------
df = 1
llstar = -quantile(Chisq(df),0.95)/2 # 95% asymptotic threshold
nptss=40 # points on the profile
Pmin=0.001 # lower bound for the profile  
Pmax=2# upper bound for the profile 
Prange=LinRange(Pmin,Pmax,nptss) # vector of D1 values along the profile
nrange=zeros(3,nptss) # matrix to store the nuisance parameters once optimized out
llP=zeros(nptss) # loglikelihood at each point along the profile
nllP=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun1(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data1,data2,t,[Prange[i],aa[1],aa[2],aa[3]])
    end
    local lb1=[lb_ρ,lb_k1,lb_k2] # lower bound 
    local ub1=[ub_ρ,ub_k1,ub_k2] # upper bound
    local θG1=[ρmle,k1mle,k2mle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun1,θG1,lb1,ub1)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llP[i]=fo[1] # store the loglikelihood 
end
nllP=llP.-maximum(llP); # calculate normalised loglikelihood
println("Profile P Complete-----------------------------------------------------------------------------")

# Plot the univariate profile likelihood for D.
s1a=plot(Prange,nllP,lw=2,xlim=(0,2),ylim=(-3,0.1),legend=false,tickfontsize=tfs,labelfontsize=tfs,grid=false,
xticks=([0,0.5,1,1.5,2], [latexstring("0"),"", latexstring("1"),"",latexstring("2")]), yticks=([-3,-2,-1,0],[latexstring("-3"),latexstring("-2"),latexstring("-1"),latexstring("0")]),framestyle=:box,left_margin = 15mm)
s1a=hline!([llstar],lw=2)
s1a=vline!([Pmle],lw=2)
s1a = annotate!(1, -4, text(L"P", :left, tfs))
s1a = annotate!(-0.8, -1.5, text(L"\bar{\ell}_{p}", :left, tfs,rotation=90))
s1a = annotate!(-1, 0.41, text(L"(\mathrm{a})", :left, tfs))


function predicrealP(lsp,nrange,t)
    # ------------------------------------------------------------------------------------------------------------------
       # Input:
       # lsp: (vector) Parameter sets within the 95% log-likelihood threshold.
       # t: (number) Prediction interval constructed at time t.
       # Output:
       # lb1, ub1: (vectors) Lower and upper bounds of prediction interval for subpopulation 1.
       # ------------------------------------------------------------------------------------------------------------------
       # Generate empty lower and upper bounds of prediction intervals for subpopulations 1.
       xx = 0:1:199
       gp=length(xx)
       lb1 = ones(200,1).*cmax; ub1 = zeros(200,1)
       lb2 = ones(200,1).*cmax; ub2 = zeros(200,1)
       for k = 1:length(lsp) 
           # Solve the PDE with parameter set a.
           a = lsp[k]
           c1,c2 = model(t,[a,nrange[1,k],nrange[2,k],nrange[3,k]])
           # Construct prediction interval for data realizations.
           for i = 1:I
              
                c1_05 = (quantile(Poisson(c1[i]*J),[.05,.95])[1])
                c1_95 = (quantile(Poisson(c1[i]*J),[.05,.95])[2])
                c2_05 = (quantile(Poisson(c2[i]*J),[.05,.95])[1])
                c2_95 = (quantile(Poisson(c2[i]*J),[.05,.95])[2])
                if c1_05 < lb1[i] 
                    lb1[i] = c1_05
                end
        
                if c1_95 > ub1[i] 
                    ub1[i] = c1_95
                end
                if c2_05 < lb2[i] 
                    lb2[i] = c2_05
                end
        
                if c2_95 > ub2[i] 
                    ub2[i] = c2_95
                end
           end
       end
    return lb1,ub1,lb2,ub2
end

function FindinterceptP()
    UnivariateP = linear_interpolation(Prange,vec(nllP));
    g(x)=UnivariateP(x)-llstar
    ϵ=(Pmax-Pmin)/10^6
    x0=Pmle
    x1=Pmin
    x2=(x1+x0)/2;
    while abs(x1-x0) > ϵ
        x2=(x1+x0)/2;
        if g(x0)*g(x2) < 0 
            x1=x2
        else
            x0=x2
        end
    end
    PPmin = x2
    x0=Pmle
    x1=Pmax
    x2=(x1+x0)/2;
    while abs(x1-x0) > ϵ
        x2=(x1+x0)/2;
        if g(x0)*g(x2) < 0 
            x1=x2
        else
            x0=x2
        end
    end
    PPmax = x2
    return PPmin,PPmax
end

PPmin,PPmax = FindinterceptP()
println(PPmin,PPmax)
nptss=100 # points on the profile
Pmin=PPmin # lower bound for the profile  
Pmax=PPmax# upper bound for the profile 
Prange=LinRange(Pmin,Pmax,nptss) # vector of D1 values along the profile
nrange=zeros(3,nptss) # matrix to store the nuisance parameters once optimized out
llP=zeros(nptss) # loglikelihood at each point along the profile
nllP=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun1b(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data1,data2,t,[Prange[i],aa[1],aa[2],aa[3]])
    end
    local lb1b=[lb_ρ,lb_k1,lb_k2] # lower bound 
    local ub1b=[ub_ρ,ub_k1,ub_k2] # upper bound
    local θG1b=[ρmle,k1mle,k2mle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun1b,θG1b,lb1b,ub1b)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llP[i]=fo[1] # store the loglikelihood 
end
nllP=llP.-maximum(llP); # calculate normalised loglikelihood
println("Profile P in confidence set Complete-----------------------------------------------------------------------------")
# Construct prediction interval.
@time lb1,ub1,lb2,ub2=predicrealP(Prange,nrange,t)
println("profile realisation Prediction for P Complete---------------------------------------------------------------") 
lb1 = Float64.(lb1); ub1 = Float64.(ub1)
lb2 = Float64.(lb2); ub2 = Float64.(ub2)
ub_union1[:,1] .= ub1; lb_union1[:,1] .= lb1
ub_union2[:,1] .= ub2; lb_union2[:,1] .= lb2

x_range = 1:1:200
x = 1:1:200
s1b1a = plot(x_range, lb1, lw=0, fillrange=ub1, fillalpha=0.40, xlims=(1,200), ylims=(-10,50),
         color=:red, label=false, grid=true, tickfontsize=tfs, margin=10mm, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-10,0,10,20,30,40,50], ["",latexstring("0"),"","","","",latexstring("50")]),linecolor=:transparent)
s1b1a = plot!(x,c1mle.*J,color=:red,legend=false,lw=2,ls=:dash)
s1b1a = hline!([0],lw=2,ls=:dash,color=:black)
a1b1a = scatter!(1:step:200,data1,color=:red,markerstrokecolor=:red,ms=1.5)
s1b1a = annotate!(100, -27.5,text(L"i", :left, tfs))
s1b1a = annotate!(-100*0.6, 20, text(L"C^{\mathrm{o}}_{1} ", :left, tfs,rotation=90))
s1b1a = annotate!(-200*0.5, 55.5, text(L"(\mathrm{b})", :left, tfs))


s1b1b = plot(x_range, lb2, lw=0, fillrange=ub2, fillalpha=0.40, xlims=(1,200), ylims=(-10,50),
         color=:green, label=false, grid=true, tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-10,0,10,20,30,40,50], ["",latexstring("0"),"","","","",latexstring("50")]),linecolor=:transparent)
s1b1b = plot!(x,c2mle.*J,color=:green,legend=false,lw=2,ls=:dash)
s1b1b = hline!([0],lw=2,ls=:dash,color=:black)
a1b1b = scatter!(1:step:200,data2,color=:green,markerstrokecolor=:green,ms=1.5)
s1b1b = annotate!(100, -27.5,text(L"i", :left, tfs))
s1b1b = annotate!(-100*0.6, 20, text(L"C^{\mathrm{o}}_{2} ", :left, tfs,rotation=90))
s1b1b = annotate!(-200*0.5, 55.5, text(L"(\mathrm{c})", :left, tfs))

s1b2 = plot(x_range, lb1 .- c1mle.*J, lw=0, fillrange=ub1 .- c1mle.*J, fillalpha=0.40,
         color=:red, label=false, grid=true, tickfontsize=tfs, framestyle=:box,linecolor=:transparent)

s1b2 = plot!(x_range, lb2 .- c2mle.*J, lw=0, fillrange=ub2 .- c2mle.*J, fillalpha=0.40,
         color=:green, label=false, grid=true, xlims=(1,200), ylims=(-30,20), tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-30,-20,-10,0,10,20], [latexstring("-30"),"","",latexstring("0"),"",latexstring("20")]),linecolor=:transparent)

s1b2 = hline!([0],lw=2,ls=:dash,color=:black,legend=false)
s1b2 = annotate!(100, -44,text(L"i", :left, tfs))
s1b2 = annotate!(-100*0.65, -23, text(L"C_{s,0.95}^{\mathrm{o}} - \hat{C}_{s} ", :left, tfs,rotation=90))
s1b2 = annotate!(-200*0.5, 25, text(L"(\mathrm{d})", :left, tfs))

s1 = plot(s1a,s1b1a,s1b1b,s1b2,layout = (1,4), link = :y, 
bottom_margin = 6mm,top_margin=6mm, right_margin = 1mm,left_margin = 4mm,size=(580,150)) 
savefig(s1,"$path\\P.pdf")
display(s1) 


#Profile ρ-------------------------------------------------------------------
df = 1
llstar = -quantile(Chisq(df),0.95)/2 # 95% asymptotic threshold
nptss=40 # points on the profile
ρmin=0 # lower bound for the profile  
ρmax=1.5# upper bound for the profile 
ρrange=LinRange(ρmin,ρmax,nptss) # vector of D1 values along the profile
nrange=zeros(3,nptss) # matrix to store the nuisance parameters once optimized out
llρ=zeros(nptss) # loglikelihood at each point along the profile
nllρ=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun2(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data1,data2,t,[aa[1],ρrange[i],aa[2],aa[3]])
    end
    local lb2=[lb_P,lb_k1,lb_k2] # lower bound 
    local ub2=[ub_P,ub_k1,ub_k2] # upper bound
    local θG2=[Pmle,k1mle,k2mle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun2,θG2,lb2,ub2)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llρ[i]=fo[1] # store the loglikelihood 
end
nllρ=llρ.-maximum(llρ); # calculate normalised loglikelihood
println("Profile ρ Complete-----------------------------------------------------------------------------")

# Plot the univariate profile likelihood for ρ.
s2a=plot(ρrange,nllρ,lw=2,xlim=(0,2),ylim=(-3,0.1),legend=false,tickfontsize=tfs,labelfontsize=tfs,grid=false,
xticks=([0,0.5,1,1.5,2], [latexstring("0"),"", latexstring("1"),"",latexstring("2")]),framestyle=:box,left_margin = 15mm)
s2a=hline!([llstar],lw=2)
s2a=vline!([ρmle],lw=2)
s2a = annotate!(1, -4, text(L"ρ", :left, tfs))
s2a = annotate!(-0.8, -1.5, text(L"\bar{\ell}_{p}", :left, tfs,rotation=90))
s2a = annotate!(-1, 0.41, text(L"(\mathrm{e})", :left, tfs))




function predicrealρ(lsp,nrange,t)
    # ------------------------------------------------------------------------------------------------------------------
       # Input:
       # lsp: (vector) Parameter sets within the 95% log-likelihood threshold.
       # t: (number) Prediction interval constructed at time t.
       # Output:
       # lb1, ub1: (vectors) Lower and upper bounds of prediction interval for subpopulation 1.
       # ------------------------------------------------------------------------------------------------------------------
       # Generate empty lower and upper bounds of prediction intervals for subpopulations 1.
       xx = 0:1:199
       gp=length(xx)
       lb1 = ones(200,1).*cmax; ub1 = zeros(200,1)
       lb2 = ones(200,1).*cmax; ub2 = zeros(200,1)
       for k = 1:length(lsp) 
           # Solve the PDE with parameter set a.
           a = lsp[k]
           c1,c2 = model(t,[nrange[1,k],a,nrange[2,k],nrange[3,k]])
           # Construct prediction interval for data realizations.
           for i = 1:I
              
                c1_05 = (quantile(Poisson(c1[i]*J),[.05,.95])[1])
                c1_95 = (quantile(Poisson(c1[i]*J),[.05,.95])[2])
                c2_05 = (quantile(Poisson(c2[i]*J),[.05,.95])[1])
                c2_95 = (quantile(Poisson(c2[i]*J),[.05,.95])[2])
                if c1_05 < lb1[i] 
                    lb1[i] = c1_05
                end
        
                if c1_95 > ub1[i] 
                    ub1[i] = c1_95
                end
                if c2_05 < lb2[i] 
                    lb2[i] = c2_05
                end
        
                if c2_95 > ub2[i] 
                    ub2[i] = c2_95
                end
           end
       end
    return lb1,ub1,lb2,ub2
end
#Calculate location where the profile intersects the threshold log-likelihood
function Findinterceptρ()
    Univariateρ = linear_interpolation(ρrange,vec(nllρ));
    g(x)=Univariateρ(x)-llstar
    ϵ=(ρmax-ρmin)/10^6
    x0=ρmle
    x1=ρmin
    x2=(x1+x0)/2;
    while abs(x1-x0) > ϵ
        x2=(x1+x0)/2;
        if g(x0)*g(x2) < 0 
            x1=x2
        else
            x0=x2
        end
    end
    ρρmin = x2
    x0=ρmle
    x1=ρmax
    x2=(x1+x0)/2;
    while abs(x1-x0) > ϵ
        x2=(x1+x0)/2;
        if g(x0)*g(x2) < 0 
            x1=x2
        else
            x0=x2
        end
    end
    ρρmax = x2
    return ρρmin,ρρmax
end

ρρmin,ρρmax = Findinterceptρ()
println(ρρmin,ρρmax)
llstar = -quantile(Chisq(df),0.95)/2 # 95% asymptotic threshold
nptss=100 # points on the profile
ρmin=ρρmin # lower bound for the profile  
ρmax=ρρmax# upper bound for the profile 
ρrange=LinRange(ρmin,ρmax,nptss) # vector of D1 values along the profile
nrange=zeros(3,nptss) # matrix to store the nuisance parameters once optimized out
llρ=zeros(nptss) # loglikelihood at each point along the profile
nllρ=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun2b(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data1,data2,t,[aa[1],ρrange[i],aa[2],aa[3]])
    end
    local lb2b=[lb_P,lb_k1,lb_k2] # lower bound 
    local ub2b=[ub_P,ub_k1,ub_k2] # upper bound
    local θG2b=[Pmle,k1mle,k2mle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun2b,θG2b,lb2b,ub2b)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llρ[i]=fo[1] # store the loglikelihood 
end
nllρ=llρ.-maximum(llρ); # calculate normalised loglikelihood
println("Profile ρ in confidence set Complete-----------------------------------------------------------------------------")
# Construct prediction interval.
@time lb1,ub1,lb2,ub2=predicrealρ(ρrange,nrange,t)
lb1 = Float64.(lb1); ub1 = Float64.(ub1)
lb2 = Float64.(lb2); ub2 = Float64.(ub2)
ub_union1[:,2] .= ub1; lb_union1[:,2] .= lb1
ub_union2[:,2] .= ub2; lb_union2[:,2] .= lb2

x_range = 1:1:200
x = 1:1:200
s2b1a = plot(x_range, lb1, lw=0, fillrange=ub1, fillalpha=0.40, xlims=(1,200), ylims=(-10,50),
         color=:red, label=false, grid=true, tickfontsize=tfs, margin=10mm, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-10,0,10,20,30,40,50], ["",latexstring("0"),"","","","",latexstring("50")]),linecolor=:transparent)
s2b1a = plot!(x,c1mle.*J,color=:red,legend=false,lw=2,ls=:dash)
s2b1a = hline!([0],lw=2,ls=:dash,color=:black)
a2b1a = scatter!(1:step:200,data1,color=:red,markerstrokecolor=:red,ms=1.5)
s2b1a = annotate!(100, -27.5,text(L"i", :left, tfs))
s2b1a = annotate!(-100*0.6, 20, text(L"C^{\mathrm{o}}_{1} ", :left, tfs,rotation=90))
s2b1a = annotate!(-200*0.5, 55.5, text(L"(\mathrm{f})", :left, tfs))


s2b1b = plot(x_range, lb2, lw=0, fillrange=ub2, fillalpha=0.40, xlims=(1,200), ylims=(-10,50),
         color=:green, label=false, grid=true, tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-10,0,10,20,30,40,50], ["",latexstring("0"),"","","","",latexstring("50")]),linecolor=:transparent)
s2b1b = plot!(x,c2mle.*J,color=:green,legend=false,lw=2,ls=:dash)
s2b1b = hline!([0],lw=2,ls=:dash,color=:black)
a2b1b = scatter!(1:step:200,data2,color=:green,markerstrokecolor=:green,ms=1.5)
s2b1b = annotate!(100, -27.5,text(L"i", :left, tfs))
s2b1b = annotate!(-100*0.6, 20, text(L"C^{\mathrm{o}}_{2} ", :left, tfs,rotation=90))
s2b1b = annotate!(-200*0.5, 55.5, text(L"(\mathrm{g})", :left, tfs))

s2b2 = plot(x_range, lb1 .- c1mle.*J, lw=0, fillrange=ub1 .- c1mle.*J, fillalpha=0.40,
         color=:red, label=false, grid=true, tickfontsize=tfs, framestyle=:box,linecolor=:transparent)

s2b2 = plot!(x_range, lb2 .- c2mle.*J, lw=0, fillrange=ub2 .- c2mle.*J, fillalpha=0.40,
         color=:green, label=false, grid=true, xlims=(1,200), ylims=(-30,20), tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-30,-20,-10,0,10,20], [latexstring("-30"),"","",latexstring("0"),"",latexstring("20")]),linecolor=:transparent)

s2b2 = hline!([0],lw=2,ls=:dash,color=:black,legend=false)
s2b2 = annotate!(100, -44,text(L"i", :left, tfs))
s2b2 = annotate!(-100*0.65, -23, text(L"C_{s,0.95}^{\mathrm{o}} - \hat{C}_{s} ", :left, tfs,rotation=90))
s2b2 = annotate!(-200*0.5, 25, text(L"(\mathrm{h})", :left, tfs))

s2 = plot(s2a,s2b1a,s2b1b,s2b2,layout = (1,4), link = :y, 
bottom_margin = 6mm,top_margin=6mm, right_margin = 1mm,left_margin = 4mm,size=(580,150))  # Here we set the top margin
savefig(s2,"$path\\ρ.pdf")
display(s2)


#Profile k1-------------------------------------------------------------------
df = 1
llstar = -quantile(Chisq(df),0.95)/2 # 95% asymptotic threshold
nptss=40 # points on the profile
k1min=0.001 # lower bound for the profile  
k1max=0.035# upper bound for the profile 
k1range=LinRange(k1min,k1max,nptss) # vector of k1 values along the profile
nrange=zeros(3,nptss) # matrix to store the nuisance parameters once optimized out
llk1=zeros(nptss) # loglikelihood at each point along the profile
nllk1=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun3(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data1,data2,t,[aa[1],aa[2],k1range[i],aa[3]])
    end
    local lb3=[lb_P,lb_ρ,lb_k2] # lower bound 
    local ub3=[ub_P,ub_ρ,ub_k2] # upper bound
    local θG3=[Pmle,ρmle,k2mle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun3,θG3,lb3,ub3)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llk1[i]=fo[1] # store the loglikelihood 
end
nllk1=llk1.-maximum(llk1); # calculate normalised loglikelihood
println("Profile λ1 Complete-----------------------------------------------------------------------------")

# Plot the univariate profile likelihood for k1.
s3a=plot(k1range,nllk1,lw=2,xlim=(0,0.1),ylim=(-3,0.1),legend=false,tickfontsize=tfs,labelfontsize=tfs,grid=false,
xticks=([0,0.025,0.05,0.075,0.1], [latexstring("0"),"", latexstring("0.05"),"",latexstring("0.1")]), yticks=([-3,-2,-1,0],[latexstring("-3"),latexstring("-2"),latexstring("-1"),latexstring("0")]),framestyle=:box,left_margin = 15mm)
s3a=hline!([llstar],lw=2)
s3a=vline!([k1mle],lw=2)
s3a = annotate!(0.05, -4, text(L"λ_{1}", :left, tfs))
s3a = annotate!(-0.04, -1.5, text(L"\bar{\ell}_{p}", :left, tfs,rotation=90))
s3a = annotate!(-0.05, 0.41, text(L"(\mathrm{i})", :left, tfs))


function predicrealk1(lsp,nrange,t)
    # ------------------------------------------------------------------------------------------------------------------
       # Input:
       # lsp: (vector) Parameter sets within the 95% log-likelihood threshold.
       # t: (number) Prediction interval constructed at time t.
       # Output:
       # lb1, ub1: (vectors) Lower and upper bounds of prediction interval for subpopulation 1.
       # ------------------------------------------------------------------------------------------------------------------
       # Generate empty lower and upper bounds of prediction intervals for subpopulations 1.
       xx = 0:1:199
       gp=length(xx)
       lb1 = ones(200,1).*cmax; ub1 = zeros(200,1)
       lb2 = ones(200,1).*cmax; ub2 = zeros(200,1)
       for k = 1:length(lsp) 
           # Solve the PDE with parameter set a.
           a = lsp[k]
           c1,c2 = model(t,[nrange[1,k],nrange[2,k],a,nrange[3,k]])
           # Construct prediction interval for data realizations.
           for i = 1:I
              
                c1_05 = (quantile(Poisson(c1[i]*J),[.05,.95])[1])
                c1_95 = (quantile(Poisson(c1[i]*J),[.05,.95])[2])
                c2_05 = (quantile(Poisson(c2[i]*J),[.05,.95])[1])
                c2_95 = (quantile(Poisson(c2[i]*J),[.05,.95])[2])
                if c1_05 < lb1[i] 
                    lb1[i] = c1_05
                end
        
                if c1_95 > ub1[i] 
                    ub1[i] = c1_95
                end
                if c2_05 < lb2[i] 
                    lb2[i] = c2_05
                end
        
                if c2_95 > ub2[i] 
                    ub2[i] = c2_95
                end
           end
       end
    return lb1,ub1,lb2,ub2
end
#Calculate location where the profile intersects the threshold log-likelihood
function Findinterceptk1()
    Univariatek1 = linear_interpolation(k1range,vec(nllk1));
    g(x)=Univariatek1(x)-llstar
    ϵ=(k1max-k1min)/10^6
    x0=k1mle
    x1=k1min
    x2=(x1+x0)/2;
    while abs(x1-x0) > ϵ
        x2=(x1+x0)/2;
        if g(x0)*g(x2) < 0 
            x1=x2
        else
            x0=x2
        end
    end
    k1k1min = x2
    x0=k1mle
    x1=k1max
    x2=(x1+x0)/2;
    while abs(x1-x0) > ϵ
        x2=(x1+x0)/2;
        if g(x0)*g(x2) < 0 
            x1=x2
        else
            x0=x2
        end
    end
    k1k1max = x2
    return k1k1min,k1k1max
end

k1k1min,k1k1max = Findinterceptk1()
println(k1k1min,k1k1max)
nptss=100 # points on the profile
k1min=k1k1min # lower bound for the profile  
k1max=k1k1max # upper bound for the profile 
k1range=LinRange(k1min,k1max,nptss) # vector of k1 values along the profile
nrange=zeros(3,nptss) # matrix to store the nuisance parameters once optimized out
llk1=zeros(nptss) # loglikelihood at each point along the profile
nllk1=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun3b(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data1,data2,t,[aa[1],aa[2],k1range[i],aa[3]])
    end
    local lb3b=[lb_P,lb_ρ,lb_k2] # lower bound 
    local ub3b=[ub_P,ub_ρ,ub_k2] # upper bound
    local θG3b=[Pmle,ρmle,k2mle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun3b,θG3b,lb3b,ub3b)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llk1[i]=fo[1] # store the loglikelihood 
end
nllk1=llk1.-maximum(llk1); # calculate normalised loglikelihood
println("Profile λ1 in confidence set Complete-----------------------------------------------------------------------------")
# Construct prediction interval.
@time lb1,ub1,lb2,ub2=predicrealk1(k1range,nrange,t)
lb1 = Float64.(lb1); ub1 = Float64.(ub1)
lb2 = Float64.(lb2); ub2 = Float64.(ub2)
ub_union1[:,3] .= ub1; lb_union1[:,3] .= lb1
ub_union2[:,3] .= ub2; lb_union2[:,3] .= lb2

x_range = 1:1:200
x = 1:1:200
s3b1a = plot(x_range, lb1, lw=0, fillrange=ub1, fillalpha=0.40, xlims=(1,200), ylims=(-10,50),
         color=:red, label=false, grid=true, tickfontsize=tfs, margin=10mm, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-10,0,10,20,30,40,50], ["",latexstring("0"),"","","","",latexstring("50")]),linecolor=:transparent)
s3b1a = plot!(x,c1mle.*J,color=:red,legend=false,lw=2,ls=:dash)
s3b1a = hline!([0],lw=2,ls=:dash,color=:black)
a3b1a = scatter!(1:step:200,data1,color=:red,markerstrokecolor=:red,ms=1.5)
s3b1a = annotate!(100, -27.5,text(L"i", :left, tfs))
s3b1a = annotate!(-100*0.6, 20, text(L"C^{\mathrm{o}}_{1} ", :left, tfs,rotation=90))
s3b1a = annotate!(-200*0.5, 55.5, text(L"(\mathrm{j})", :left, tfs))


s3b1b = plot(x_range, lb2, lw=0, fillrange=ub2, fillalpha=0.40, xlims=(1,200), ylims=(-10,50),
         color=:green, label=false, grid=true, tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-10,0,10,20,30,40,50], ["",latexstring("0"),"","","","",latexstring("50")]),linecolor=:transparent)
s3b1b = plot!(x,c2mle.*J,color=:green,legend=false,lw=2,ls=:dash)
s3b1b = hline!([0],lw=2,ls=:dash,color=:black)
a3b1b = scatter!(1:step:200,data2,color=:green,markerstrokecolor=:green,ms=1.5)
s3b1b = annotate!(100, -27.5,text(L"i", :left, tfs))
s3b1b = annotate!(-100*0.6, 20, text(L"C^{\mathrm{o}}_{2} ", :left, tfs,rotation=90))
s3b1b = annotate!(-200*0.5, 55.5, text(L"(\mathrm{k})", :left, tfs))

s3b2 = plot(x_range, lb1 .- c1mle.*J, lw=0, fillrange=ub1 .- c1mle.*J, fillalpha=0.40,
         color=:red, label=false, grid=true, tickfontsize=tfs, framestyle=:box,linecolor=:transparent)

s3b2 = plot!(x_range, lb2 .- c2mle.*J, lw=0, fillrange=ub2 .- c2mle.*J, fillalpha=0.40,
         color=:green, label=false, grid=true, xlims=(1,200), ylims=(-30,20), tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-30,-20,-10,0,10,20], [latexstring("-30"),"","",latexstring("0"),"",latexstring("20")]),linecolor=:transparent)

s3b2 = hline!([0],lw=2,ls=:dash,color=:black,legend=false)
s3b2 = annotate!(100, -44,text(L"i", :left, tfs))
s3b2 = annotate!(-100*0.65, -23, text(L"C_{s,0.95}^{\mathrm{o}} - \hat{C}_{s} ", :left, tfs,rotation=90))
s3b2 = annotate!(-200*0.5, 25, text(L"(\mathrm{l})", :left, tfs))

s3 = plot(s3a,s3b1a,s3b1b,s3b2,layout = (1,4), link = :y, 
bottom_margin = 6mm,top_margin=6mm, right_margin = 1mm,left_margin = 4mm,size=(580,150))  # Here we set the top margin
savefig(s3,"$path\\λ1.pdf")
display(s3)


#Profile λ2-------------------------------------------------------------------
df = 1
llstar = -quantile(Chisq(df),0.95)/2 # 95% asymptotic threshold
nptss=40 # points on the profile
k2min=0.001 # lower bound for the profile  
k2max=0.08# upper bound for the profile 
k2range=LinRange(k2min,k2max,nptss) # vector of k2 values along the profile
nrange=zeros(3,nptss) # matrix to store the nuisance parameters once optimized out
llk2=zeros(nptss) # loglikelihood at each point along the profile
nllk2=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun4(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data1,data2,t,[aa[1],aa[2],aa[3],k2range[i]])
    end
    local lb4=[lb_P,lb_ρ,lb_k1] # lower bound 
    local ub4=[ub_P,ub_ρ,ub_k1] # upper bound
    local θG4=[Pmle,ρmle,k1mle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun4,θG4,lb4,ub4)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llk2[i]=fo[1] # store the loglikelihood 
end
nllk2=llk2.-maximum(llk2); # calculate normalised loglikelihood
println("Profile λ2 Complete-----------------------------------------------------------------------------")

# Plot the univariate profile likelihood for k2.
s4a=plot(k2range,nllk2,lw=2,xlim=(0,0.1),ylim=(-3,0.1),legend=false,tickfontsize=tfs,labelfontsize=tfs,grid=false,
xticks=([0,0.025,0.05,0.075,0.1], [latexstring("0"),"", latexstring("0.05"),"",latexstring("0.1")]), yticks=([-3,-2,-1,0],[latexstring("-3"),latexstring("-2"),latexstring("-1"),latexstring("0")]),framestyle=:box,left_margin = 15mm)
s4a=hline!([llstar],lw=2)
s4a=vline!([k2mle],lw=2)
s4a = annotate!(0.05, -4, text(L"λ_{2}", :left, tfs))
s4a = annotate!(-0.04, -1.5, text(L"\bar{\ell}_{p}", :left, tfs,rotation=90))
s4a = annotate!(-0.05, 0.41, text(L"(\mathrm{m})", :left, tfs))

function predicrealk2(lsp,nrange,t)
    # ------------------------------------------------------------------------------------------------------------------
       # Input:
       # lsp: (vector) Parameter sets within the 95% log-likelihood threshold.
       # t: (number) Prediction interval constructed at time t.
       # Output:
       # lb1, ub1: (vectors) Lower and upper bounds of prediction interval for subpopulation 1.
       # ------------------------------------------------------------------------------------------------------------------
       # Generate empty lower and upper bounds of prediction intervals for subpopulations 1.
       xx = 0:1:199
       gp=length(xx)
       lb1 = ones(200,1).*cmax; ub1 = zeros(200,1)
       lb2 = ones(200,1).*cmax; ub2 = zeros(200,1)
       for k = 1:length(lsp) 
           # Solve the PDE with parameter set a.
           a = lsp[k]
           c1,c2 = model(t,[nrange[1,k],nrange[2,k],nrange[3,k],a])
           # Construct prediction interval for data realizations.
           for i = 1:I
              
                c1_05 = (quantile(Poisson(c1[i]*J),[.05,.95])[1])
                c1_95 = (quantile(Poisson(c1[i]*J),[.05,.95])[2])
                c2_05 = (quantile(Poisson(c2[i]*J),[.05,.95])[1])
                c2_95 = (quantile(Poisson(c2[i]*J),[.05,.95])[2])
                if c1_05 < lb1[i] 
                    lb1[i] = c1_05
                end
        
                if c1_95 > ub1[i] 
                    ub1[i] = c1_95
                end
                if c2_05 < lb2[i] 
                    lb2[i] = c2_05
                end
        
                if c2_95 > ub2[i] 
                    ub2[i] = c2_95
                end
           end
       end
    return lb1,ub1,lb2,ub2
end
#Calculate location where the profile intersects the threshold log-likelihood
function Findinterceptk2()
    Univariatek2 = linear_interpolation(k2range,vec(nllk2));
    g(x)=Univariatek2(x)-llstar
    ϵ=(k2max-k2min)/10^6
    x0=k2mle
    x1=k2min
    x2=(x1+x0)/2;
    while abs(x1-x0) > ϵ
        x2=(x1+x0)/2;
        if g(x0)*g(x2) < 0 
            x1=x2
        else
            x0=x2
        end
    end
    k2k2min = x2
    x0=k2mle
    x1=k2max
    x2=(x1+x0)/2;
    while abs(x1-x0) > ϵ
        x2=(x1+x0)/2;
        if g(x0)*g(x2) < 0 
            x1=x2
        else
            x0=x2
        end
    end
    k2k2max = x2
    return k2k2min,k2k2max
end

k2k2min,k2k2max = Findinterceptk2()
println(k2k2min,k2k2max)
nptss=100 # points on the profile
k2min=k2k2min # lower bound for the profile  
k2max=k2k2max# upper bound for the profile 
k2range=LinRange(k2min,k2max,nptss) # vector of k2 values along the profile
nrange=zeros(3,nptss) # matrix to store the nuisance parameters once optimized out
llk2=zeros(nptss) # loglikelihood at each point along the profile
nllk2=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun4b(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data1,data2,t,[aa[1],aa[2],aa[3],k2range[i]])
    end
    local lb4b=[lb_P,lb_ρ,lb_k1] # lower bound 
    local ub4b=[ub_P,ub_ρ,ub_k1] # upper bound
    local θG4b=[Pmle,ρmle,k1mle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun4b,θG4b,lb4b,ub4b)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llk2[i]=fo[1] # store the loglikelihood 
end
nllk2=llk2.-maximum(llk2); # calculate normalised loglikelihood
println("Profile λ2 in confidence set Complete-----------------------------------------------------------------------------")
# Construct prediction interval.
@time lb1,ub1,lb2,ub2=predicrealk2(k2range,nrange,t)
lb1 = Float64.(lb1); ub1 = Float64.(ub1)
lb2 = Float64.(lb2); ub2 = Float64.(ub2)
ub_union1[:,4] .= ub1; lb_union1[:,4] .= lb1
ub_union2[:,4] .= ub2; lb_union2[:,4] .= lb2

x_range = 1:1:200
x = 1:1:200
s4b1a = plot(x_range, lb1, lw=0, fillrange=ub1, fillalpha=0.40, xlims=(1,200), ylims=(-10,50),
         color=:red, label=false, grid=true, tickfontsize=tfs, margin=10mm, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-10,0,10,20,30,40,50], ["",latexstring("0"),"","","","",latexstring("50")]),linecolor=:transparent)
s4b1a = plot!(x,c1mle.*J,color=:red,legend=false,lw=2,ls=:dash)
s4b1a = hline!([0],lw=2,ls=:dash,color=:black)
a4b1a = scatter!(1:step:200,data1,color=:red,markerstrokecolor=:red,ms=1.5)
s4b1a = annotate!(100, -27.5,text(L"i", :left, tfs))
s4b1a = annotate!(-100*0.6, 20, text(L"C^{\mathrm{o}}_{1} ", :left, tfs,rotation=90))
s4b1a = annotate!(-200*0.5, 55.5, text(L"(\mathrm{n})", :left, tfs))


s4b1b = plot(x_range, lb2, lw=0, fillrange=ub2, fillalpha=0.40, xlims=(1,200), ylims=(-10,50),
         color=:green, label=false, grid=true, tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-10,0,10,20,30,40,50], ["",latexstring("0"),"","","","",latexstring("50")]),linecolor=:transparent)
s4b1b = plot!(x,c2mle.*J,color=:green,legend=false,lw=2,ls=:dash)
s4b1b = hline!([0],lw=2,ls=:dash,color=:black)
a4b1b = scatter!(1:step:200,data2,color=:green,markerstrokecolor=:green,ms=1.5)
s4b1b = annotate!(100, -27.5,text(L"i", :left, tfs))
s4b1b = annotate!(-100*0.6, 20, text(L"C^{\mathrm{o}}_{2} ", :left, tfs,rotation=90))
s4b1b = annotate!(-200*0.5, 55.5, text(L"(\mathrm{o})", :left, tfs))

s4b2 = plot(x_range, lb1 .- c1mle.*J, lw=0, fillrange=ub1 .- c1mle.*J, fillalpha=0.40,
         color=:red, label=false, grid=true, tickfontsize=tfs, framestyle=:box,linecolor=:transparent)

s4b2 = plot!(x_range, lb2 .- c2mle.*J, lw=0, fillrange=ub2 .- c2mle.*J, fillalpha=0.40,
         color=:green, label=false, grid=true, xlims=(1,200), ylims=(-30,20), tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-30,-20,-10,0,10,20], [latexstring("-30"),"","",latexstring("0"),"",latexstring("20")]),linecolor=:transparent)

s4b2 = hline!([0],lw=2,ls=:dash,color=:black,legend=false)
s4b2 = annotate!(100, -44,text(L"i", :left, tfs))
s4b2 = annotate!(-100*0.65, -23, text(L"C_{s,0.95}^{\mathrm{o}} - \hat{C}_{s} ", :left, tfs,rotation=90))
s4b2 = annotate!(-200*0.5, 25, text(L"(\mathrm{p})", :left, tfs))

s4 = plot(s4a,s4b1a,s4b1b,s4b2,layout = (1,4), link = :y, 
bottom_margin = 6mm,top_margin=6mm, right_margin = 1mm,left_margin = 4mm,size=(580,150))  # Here we set the top margin
savefig(s4,"$path\\λ2.pdf")
display(s4)