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


# In this script, we generate the result for Case 2b
# Variable P indicate P^\star in the main document
# Variable k indicate λ^\star in the main document
# Variable c0 indicate c_0 in the main document

# Get the path of the current script.
path = dirname(@__FILE__)

# Read the CSV file back into a DataFrame
df = CSV.read("$path\\data.csv", DataFrame)

# Convert the DataFrame back into a matrix
step =  10
tfs = 8
data = df.a[1:step:end]
length(data)

# Fixed parameters
r0 = 20; R=99.5; δ= 0.5
# Data is collected at time t = t1.
t = 600
I = 200; J = 20
# union of prediction interval
ub_union1= zeros(I,4)
lb_union1= zeros(I,4)

# diff!(): Discretized PDE 
function diff!(d_u, u, p, t)

    (D,k,δ) = p
    R = length(u) - 1
    r = 0:δ:R
    N = R + 1

    # Boundary at x = 0, associated with Equation (24).
    d_u[1] = (D/(δ^2))*(u[2]-u[1]) + k*u[1]*(1-u[1])
    

    # Associated with Equation (25).
    for n in 2:N-1
        d_u[n] = (D/(δ^2))*(u[n+1]-2*u[n]+u[n-1]) + (D/(2*δ*r[n]))*(u[n+1]-u[n-1]) + k*u[n]*(1-u[n])
    end
    
    # Boundary at x = 199, associated with Equation (26).
    d_u[N] =  (D/(δ^2))*(-u[N]+u[N-1]) + k*u[N]*(1-u[N])
   

    return d_u
end

# pdesolver(): Solves the ODE in diff!().
function pdesolver(t,P,k,r0,c_0,δ,R)
    # ------------------------------------------------------------------------------------------------------------------
    # Input:
    # Time: (number) Solve ODE/PDE at t = Time.
    # D: (number) Diffusivity.
    # k: (number) proliferation rate.
    # r0: (number) intial radius.
    # c0: (number) intial density.
    # δ: (number) Grid spacing.
    # R: (number) ODE/PDE is solved on 0 <= x <= R.
    # Output:
    # c1, c2: Solution of ODE/PDE at time t = Time.
    # ------------------------------------------------------------------------------------------------------------------
    # Construct initial condition for PDE.
    D = P/4
    # Total number of mesh points N.
    N = Int(R/δ + 1)
    # Initial condition u0.
    c0 = zeros(1,N)
    # Update initial condition u0: c1(x) = 1 at 79 <= x <= 119.
    for n=1:N
       xx = (n-1)*δ
       if xx <= r0
          c0[n] = c_0
       end
    end


    # Return the inital conditon if ODE/PDE is solved at t=0
    if t == 0
        c0 = vec(c0)
        return c0
    end

    # Solve the PDE using Heun's method at t=Time 
    p=(D,k,δ)
    tspan=(0,t)
    prob=ODEProblem(diff!,c0,tspan,p)
    alg=Heun() 
    sol=solve(prob,alg,saveat=t);
    sol = sol[end]
    c = vec(sol)
    return c
end 

# model(): The continuum model
function model(t,a) 
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # t: (number) Analytic solution of PDE at t.
    # a: (vector) Parameter vector.
    #------------------------------------------------------------------------------------------------------------------
    # Evaluate the analytic solution
    
    c05 = pdesolver(t,a[1],a[2],r0,a[3],δ,R)
    c50 = reverse(c05)
    l = length(c05)
    c = zeros(l+l-1,1)
    c[1:l-1] .= c50[1:l-1]
    c[l:end] .= c05[1:l]
    return c
end

# error(): The loglikelihood function for multinomial measurement error model for Case 1
function error(data,t,a)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # data1: (vector) Count data for subpopulation 1.
    # δ: (number) Grid spacing.
    # L: (number) PDE solved on 0 <= x <= L.
    # a: (vector) Parameter vector.
    # Output:
    # e: (number) Log-likelihood for the Additive Gaussian Measurement Error Model (Case 1) with parameter vector a.
    #------------------------------------------------------------------------------------------------------------------
    # Solve the PDE at time t1 with parameter vector a.
    c=model(t,a) 
    xlocdata = 0:δ:199
    interpr = linear_interpolation(xlocdata,vec(c));
    x = 0:step:199
    c = interpr(x)
    # Estimate the loglikelihood function.
    e=0.0;
    for i = 1:length(data)
        e = e + data[i]*log(c[i]) + (J-data[i])*log((1-c[i]))
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
    return error(data,t,a)
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
    opt.maxtime = 10*600
    res = optimize(opt,θ₀)
    return res[[2,1]]
end






#Expected Parameters
P = 1; k = 0.01; c_0 = 0.8;
θ = [P,k,c_0]
l_theta = fun(θ)
#MLE----------------------------------------------------------------------------
# Inital guess.
Pg = 1; kg = 0.01; c0g = 0.8;
lb_P = 0.001; lb_k = 0 ; lb_c0 = 0;
ub_P = 2; ub_k = 0.05; ub_c0 = 1;
θG = [Pg,kg,c0g] # inital guess
lb=[lb_P,lb_k,lb_c0] # lower bound
ub=[ub_P,ub_k,ub_c0] # upper bound
# Call numerical optimization routine to give the vector of parameters xopt, and the maximum loglikelihood fopt.
@time (xopt,fopt)  = optimise(fun,θG,lb,ub)
fmle=fopt
# Print MLE parameters
Pmle=xopt[1]; 
kmle=xopt[2]; 
c0mle=xopt[3];


println("Pmle: ", Pmle)
println("λmle: ", kmle)
println("c0mle: ", c0mle)
println("likelihood at θ: ",l_theta)
println("likelihood at MLE: ",fmle)
cmle = model(t,xopt) 
fcmle = linear_interpolation(1:δ:200,vec(cmle)); 
scatter(1:step:200, data)
plot!(1:δ:200,cmle.*J)


#Profile P-------------------------------------------------------------------
df = 1
llstar = -quantile(Chisq(df),0.95)/2 # 95% asymptotic threshold
nptss=40 # points on the profile
Pmin=0.8 # lower bound for the profile  
Pmax=2# upper bound for the profile 
Prange=LinRange(Pmin,Pmax,nptss) # vector of D1 values along the profile
nrange=zeros(2,nptss) # matrix to store the nuisance parameters once optimized out
llP=zeros(nptss) # loglikelihood at each point along the profile
nllP=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun1(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data,t,[Prange[i],aa[1],aa[2]])
    end
    local lb1=[lb_k,lb_c0] # lower bound 
    local ub1=[ub_k,ub_c0] # upper bound
    local θG1=[kmle,c0mle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun1,θG1,lb1,ub1)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llP[i]=fo[1] # store the loglikelihood 
end
nllP=llP.-maximum(llP); # calculate normalised loglikelihood
println("Profile P Complete-----------------------------------------------------------------------------")

# Plot the univariate profile likelihood for D.
spl=Spline1D(Prange,nllP.-maximum(nllP),w=ones(length(Prange)),k=3,bc="nearest",s=0)
yyP=evaluate(spl,Prange)
s1a=plot(Prange,yyP,lw=2,xlim=(0,2),ylim=(-3,0.1),legend=false,tickfontsize=tfs,labelfontsize=tfs,grid=false,
xticks=([0,0.5,1,1.5,2], [latexstring("0"),"", latexstring("1"),"",latexstring("2")]), yticks=([-3,-2,-1,0],[latexstring("-3"),latexstring("-2"),latexstring("-1"),latexstring("0")]),framestyle=:box,left_margin = 15mm)
s1a=hline!([llstar],lw=2)
s1a=vline!([Pmle],lw=2)
s1a = annotate!(1, -4, text(L"P^\star", :left, tfs))
s1a = annotate!(-0.548, -1.5, text(L"\bar{\ell}_{p}", :left, tfs,rotation=90))
s1a = annotate!(-0.8, 0.41, text(L"(\mathrm{a})", :left, tfs))

function predicrealP(lsp,nrange,t)
    # ------------------------------------------------------------------------------------------------------------------
       # Input:
       # lsp: (vector) Parameter sets within the 95% log-likelihood threshold.
       # t: (number) Prediction interval constructed at time t.
       # Output:
       # lb1, ub1: (vectors) Lower and upper bounds of prediction interval for subpopulation 1.
       # ------------------------------------------------------------------------------------------------------------------
       # Generate empty lower and upper bounds of prediction intervals for subpopulations 1.
       xx = 0:δ:199
       x = 0:1:199
       gp=length(x)
       lb = ones(gp,1); ub = zeros(gp,1)
       for k = 1:length(lsp) 
           # Solve the PDE with parameter set a.
           a = lsp[k]
           c = model(t,[a,nrange[1,k],nrange[2,k]])
           p = linear_interpolation(xx,vec(c));
           # Construct prediction interval for data realizations.
           for i = 1:gp
                c_05 = (quantile(Binomial(J,p(x[i])),[.05,.95])[1])/J
                c_95 = (quantile(Binomial(J,p(x[i])),[.05,.95])[2])/J
                if c_05 < lb[i] 
                    lb[i] = c_05
                end
        
                if c_95 > ub[i] 
                    ub[i] = c_95
                end
           end
       end
    return lb,ub
end

function FindinterceptP(nllP)
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

PPmin,PPmax = FindinterceptP(yyP)
println(PPmin,PPmax)
nptss=100 # points on the profile
Pmin= PPmin# lower bound for the profile  
Pmax= PPmax# upper bound for the profile 
Prange=LinRange(Pmin,Pmax,nptss) # vector of D1 values along the profile
nrange=zeros(2,nptss) # matrix to store the nuisance parameters once optimized out
llP=zeros(nptss) # loglikelihood at each point along the profile
nllP=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun1b(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data,t,[Prange[i],aa[1],aa[2]])
    end
    local lb1b=[lb_k,lb_c0] # lower bound 
    local ub1b=[ub_k,ub_c0] # upper bound
    local θG1b=[kmle,c0mle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun1b,θG1b,lb1b,ub1b)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llP[i]=fo[1] # store the loglikelihood 
end
nllP=llP.-maximum(llP); # calculate normalised loglikelihood
println("Profile P in confidence set Complete-----------------------------------------------------------------------------")
# Construct prediction interval.
@time lb,ub=predicrealP(Prange,nrange,t)
println("profile realisation Prediction for P Complete---------------------------------------------------------------") 
lb = Float64.(lb); ub = Float64.(ub)



x_range = 1:1:200
x = 1:δ:200
s1b1 = plot(x,cmle.*J,color=:red,legend=false,lw=2,ls=:dash)
s1b1 = plot!(x_range, lb.*J, lw=0, fillrange=ub.*J, fillalpha=0.40, xlims=(1,200), ylims=(-5,25),
         color=:red, label=false, grid=true, tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-5,0,10,20,25], ["",latexstring("0"),"",latexstring("20")]),linecolor=:transparent)
s1b1 = hline!([0],lw=2,ls=:dash,color=:black)
s1b1 = hline!([20],lw=2,ls=:dash,color=:black)
s1b1 = scatter!(1:step:200,data,color=:red,markerstrokecolor=:red,ms=1.5)
s1b1 = annotate!(100, -13.7,text(L"i", :left, tfs))
s1b1 = annotate!(-100*0.6, 10, text(L"C^{\mathrm{o}} ", :left, tfs,rotation=90))
s1b1 = annotate!(-200*0.4, 28, text(L"(\mathrm{b})", :left, tfs))



s1b2 = plot(x_range, lb.*J .- fcmle(x_range).*J, lw=0, fillrange=ub.*J .- fcmle(x_range).*J, fillalpha=0.40,
         color=:red, label=false, grid=true, xlims=(1,200), ylims=(-10,10), tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-10,-5,0,5,10], [latexstring("-10"),"",latexstring("0"),"",latexstring("10")]),linecolor=:transparent)
s1b2 = hline!([0],lw=2,ls=:dash,color=:black,legend=false)
s1b2 = annotate!(100, -16,text(L"i", :left, tfs))
s1b2 = annotate!(-100*0.6, -7, text(L"C_{0.9}^{\mathrm{o}} - \hat{C} ", :left, tfs,rotation=90))
s1b2 = annotate!(-200*0.4, 12, text(L"(\mathrm{c})", :left, tfs))


s1 = plot(s1a,s1b1,s1b2,layout = (1,3), link = :y, 
             bottom_margin = 6mm,top_margin=6mm, right_margin = 2mm,left_margin = 5mm,size=(580,150))  # Here we set the top margin
savefig(s1,"$path\\p.pdf")
display(s1)



#Profile λ-------------------------------------------------------------------
df = 1
llstar = -quantile(Chisq(df),0.95)/2 # 95% asymptotic threshold
nptss=40 # points on the profile
kmin=0 # lower bound for the profile  
kmax=0.015# upper bound for the profile 
krange=LinRange(kmin,kmax,nptss) # vector of D1 values along the profile
nrange=zeros(2,nptss) # matrix to store the nuisance parameters once optimized out
llk=zeros(nptss) # loglikelihood at each point along the profile
nllk=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun2(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data,t,[aa[1],krange[i],aa[2]])
    end
    local lb2=[lb_P,lb_c0] # lower bound 
    local ub2=[ub_P,ub_c0] # upper bound
    local θG2=[Pmle,c0mle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun2,θG2,lb2,ub2)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llk[i]=fo[1] # store the loglikelihood 
end
nllk=llk.-maximum(llk); # calculate normalised loglikelihood
println("Profile λ Complete-----------------------------------------------------------------------------")

# Plot the univariate profile likelihood for k.
spl=Spline1D(krange,nllk.-maximum(nllk),w=ones(length(krange)),k=3,bc="nearest",s=0.4)
yyk=evaluate(spl,krange)
s2a=plot(krange,yyk,lw=2,xlim=(0,0.02),ylim=(-3,0.1),legend=false,tickfontsize=tfs,labelfontsize=tfs,grid=false,
xticks=([0,0.005,0.01,0.015,0.02], [latexstring("0"),"", latexstring("0.01"),"",latexstring("0.02")]), yticks=([-3,-2,-1,0],[latexstring("-3"),latexstring("-2"),latexstring("-1"),latexstring("0")]),framestyle=:box,left_margin = 15mm)
s2a=hline!([llstar],lw=2)
s2a=vline!([kmle],lw=2)
s2a = annotate!(0.01, -4, text(L"λ^\star", :left, tfs))
s2a = annotate!(-0.0055, -1.5, text(L"\bar{\ell}_{p}", :left, tfs,rotation=90))
#= s2a = annotate!(-0.008, 0.41, text(L"(\mathrm{d})", :left, tfs))
 =#


function predicrealk(lsp,nrange,t)
    # ------------------------------------------------------------------------------------------------------------------
       # Input:
       # lsp: (vector) Parameter sets within the 95% log-likelihood threshold.
       # t: (number) Prediction interval constructed at time t.
       # Output:
       # lb1, ub1: (vectors) Lower and upper bounds of prediction interval for subpopulation 1.
       # ------------------------------------------------------------------------------------------------------------------
       # Generate empty lower and upper bounds of prediction intervals for subpopulations 1.
       xx = 0:δ:199
       x = 0:1:199
       gp=length(x)
       lb = ones(gp,1); ub = zeros(gp,1)
       for k = 1:length(lsp) 
           # Solve the PDE with parameter set a.
           a = lsp[k]
           c = model(t,[nrange[1,k],a,nrange[2,k]])
           p = linear_interpolation(xx,vec(c));
           # Construct prediction interval for data realizations.
           for i = 1:gp
                c_05 = (quantile(Binomial(J,p(x[i])),[.05,.95])[1])/J
                c_95 = (quantile(Binomial(J,p(x[i])),[.05,.95])[2])/J
                if c_05 < lb[i] 
                    lb[i] = c_05
                end
        
                if c_95 > ub[i] 
                    ub[i] = c_95
                end
           end
       end
    return lb,ub
end

function Findinterceptk(nllk)
    Univariatek = linear_interpolation(krange,vec(nllk));
    g(x)=Univariatek(x)-llstar
    ϵ=(kmax-kmin)/10^6
    x0=kmle
    x1=kmin
    x2=(x1+x0)/2;
    while abs(x1-x0) > ϵ
        x2=(x1+x0)/2;
        if g(x0)*g(x2) < 0 
            x1=x2
        else
            x0=x2
        end
    end
    kkmin = x2
    x0=kmle
    x1=kmax
    x2=(x1+x0)/2;
    while abs(x1-x0) > ϵ
        x2=(x1+x0)/2;
        if g(x0)*g(x2) < 0 
            x1=x2
        else
            x0=x2
        end
    end
    kkmax = x2
    return kkmin,kkmax
end

kkmin,kkmax = Findinterceptk(yyk)
println(kkmin,kkmax)
nptss=100 # points on the profile
kmin=kkmin # lower bound for the profile  
kmax=kkmax# upper bound for the profile 
krange=LinRange(kmin,kmax,nptss) # vector of D1 values along the profile
nrange=zeros(2,nptss) # matrix to store the nuisance parameters once optimized out
llk=zeros(nptss) # loglikelihood at each point along the profile
nllk=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun2b(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data,t,[aa[1],krange[i],aa[2]])
    end
    local lb2b=[lb_P,lb_c0] # lower bound 
    local ub2b=[ub_P,ub_c0] # upper bound
    local θG2b=[Pmle,c0mle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun2b,θG2b,lb2b,ub2b)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llk[i]=fo[1] # store the loglikelihood 
end
nllk=llk.-maximum(llk); # calculate normalised loglikelihood

println("Profile λ in confidence set Complete-----------------------------------------------------------------------------")
# Construct prediction interval.
@time lb,ub=predicrealk(krange,nrange,t)
println("profile realisation Prediction for λ Complete---------------------------------------------------------------") 
lb = Float64.(lb); ub = Float64.(ub)



x_range = 1:1:200
x = 1:δ:200
s2b1 = plot(x,cmle.*J,color=:red,legend=false,lw=2,ls=:dash)
s2b1 = plot!(x_range, lb.*J, lw=0, fillrange=ub.*J, fillalpha=0.40, xlims=(1,200), ylims=(-5,25),
         color=:red, label=false, grid=true, tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-5,0,10,20,25], ["",latexstring("0"),"",latexstring("20")]),linecolor=:transparent)
s2b1 = hline!([0],lw=2,ls=:dash,color=:black)
s2b1 = hline!([20],lw=2,ls=:dash,color=:black)
s2b1 = scatter!(1:step:200,data,color=:red,markerstrokecolor=:red,ms=1.5)
s2b1 = annotate!(100, -13.7,text(L"i", :left, tfs))
s2b1 = annotate!(-100*0.6, 10, text(L"C^{\mathrm{o}} ", :left, tfs,rotation=90))
#= s2b1 = annotate!(-200*0.4, 28, text(L"(\mathrm{e})", :left, tfs))
 =#


s2b2 = plot(x_range, lb.*J .- fcmle(x_range).*J, lw=0, fillrange=ub.*J .- fcmle(x_range).*J, fillalpha=0.40,
         color=:red, label=false, grid=true, xlims=(1,200), ylims=(-10,10), tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-10,-5,0,5,10], [latexstring("-10"),"",latexstring("0"),"",latexstring("10")]),linecolor=:transparent)
s2b2 = hline!([0],lw=2,ls=:dash,color=:black,legend=false)
s2b2 = annotate!(100, -16,text(L"i", :left, tfs))
s2b2 = annotate!(-100*0.6, -7, text(L"C_{0.9}^{\mathrm{o}} - \hat{C} ", :left, tfs,rotation=90))
#= s2b2 = annotate!(-200*0.4, 12, text(L"(\mathrm{f})", :left, tfs)) =#


s2 = plot(s2a,s2b1,s2b2,layout = (1,3), link = :y, 
             bottom_margin = 6mm,top_margin=6mm, right_margin = 2mm,left_margin = 5mm,size=(580,150))  # Here we set the top margin
savefig(s2,"$path\\λ.pdf")
display(s2)



#Profile c0-------------------------------------------------------------------
df = 1
llstar = -quantile(Chisq(df),0.95)/2 # 95% asymptotic threshold
nptss=40 # points on the profile
c0min=0.000000704 # lower bound for the profile  
c0max=1# upper bound for the profile 
c0range=LinRange(c0min,c0max,nptss) # vector of D1 values along the profile
nrange=zeros(2,nptss) # matrix to store the nuisance parameters once optimized out
llc0=zeros(nptss) # loglikelihood at each point along the profile
nllc0=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun3(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data,t,[aa[1],aa[2],c0range[i]])
    end
    local lb3=[lb_P,lb_k] # lower bound 
    local ub3=[ub_P,ub_k] # upper bound
    local θG3=[Pmle,kmle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun3,θG3,lb3,ub3)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llc0[i]=fo[1] # store the loglikelihood 
end
nllc0=llc0.-maximum(llc0); # calculate normalised loglikelihood
println("Profile c0 Complete-----------------------------------------------------------------------------")

# Plot the univariate profile likelihood for k.
spl=Spline1D(c0range,nllc0.-maximum(nllc0),w=ones(length(c0range)),k=3,bc="nearest",s=0.1)
yyc0=evaluate(spl,c0range)
s3a=plot(c0range,yyc0,lw=2,xlim=(0,1),ylim=(-3,0.1),legend=false,tickfontsize=tfs,labelfontsize=tfs,grid=false,
xticks=([0,0.25,0.5,0.75,1], [latexstring("0"),"", latexstring("0.5"),"",latexstring("1")]), yticks=([-3,-2,-1,0],[latexstring("-3"),latexstring("-2"),latexstring("-1"),latexstring("0")]),framestyle=:box,left_margin = 15mm)
s3a=hline!([llstar],lw=2)
s3a=vline!([c0mle],lw=2)
s3a = annotate!(0.5, -3.9, text(L"c_{0}", :left, tfs))
s3a = annotate!(-0.27, -1.5, text(L"\bar{\ell}_{p}", :left, tfs,rotation=90))
#= s3a = annotate!(-0.4, 0.41, text(L"(\mathrm{g})", :left, tfs)) =#

function predicrealc0(lsp,nrange,t)
    # ------------------------------------------------------------------------------------------------------------------
       # Input:
       # lsp: (vector) Parameter sets within the 95% log-likelihood threshold.
       # t: (number) Prediction interval constructed at time t.
       # Output:
       # lb1, ub1: (vectors) Lower and upper bounds of prediction interval for subpopulation 1.
       # ------------------------------------------------------------------------------------------------------------------
       # Generate empty lower and upper bounds of prediction intervals for subpopulations 1.
       xx = 0:δ:199
       x = 0:1:199
       gp=length(x)
       lb = ones(gp,1); ub = zeros(gp,1)
       for k = 1:length(lsp) 
           # Solve the PDE with parameter set a.
           a = lsp[k]
           c = model(t,[nrange[1,k],nrange[2,k],a])
           p = linear_interpolation(xx,vec(c));
           # Construct prediction interval for data realizations.
           for i = 1:gp
                c_05 = (quantile(Binomial(J,p(x[i])),[.05,.95])[1])/J
                c_95 = (quantile(Binomial(J,p(x[i])),[.05,.95])[2])/J
                if c_05 < lb[i] 
                    lb[i] = c_05
                end
        
                if c_95 > ub[i] 
                    ub[i] = c_95
                end
           end
       end
    return lb,ub
end

function Findinterceptc0(nllc0)
    Univariatec0 = linear_interpolation(c0range,vec(nllc0));
    g(x)=Univariatec0(x)-llstar
    ϵ=(c0max-c0min)/10^6
    x0=c0mle
    x1=c0min
    x2=(x1+x0)/2;
    while abs(x1-x0) > ϵ
        x2=(x1+x0)/2;
        if g(x0)*g(x2) < 0 
            x1=x2
        else
            x0=x2
        end
    end
    c0c0min = x2
    x0=c0mle
    x1=c0max
    x2=(x1+x0)/2;
    while abs(x1-x0) > ϵ
        x2=(x1+x0)/2;
        if g(x0)*g(x2) < 0 
            x1=x2
        else
            x0=x2
        end
    end
    c0c0max = x2
    return c0c0min,c0c0max
end

c0c0min,c0c0max = Findinterceptc0(yyc0)

nptss=100 # points on the profile
c0min=c0c0min # lower bound for the profile  
c0max=c0c0max# upper bound for the profile 
c0range=LinRange(c0min,c0max,nptss) # vector of D1 values along the profile
nrange=zeros(2,nptss) # matrix to store the nuisance parameters once optimized out
llc0=zeros(nptss) # loglikelihood at each point along the profile
nllc0=zeros(nptss) # normalised loglikelihood at each point along the profile

@time for i in 1:nptss
    function fun3b(aa) # function to return loglikelihood by fixing the interest parameter along the profile
        return error(data,t,[aa[1],aa[2],c0range[i]])
    end
    local lb3b=[lb_P,lb_k] # lower bound 
    local ub3b=[ub_P,ub_k] # upper bound
    local θG3b=[Pmle,kmle] # initial estimate - take the MLE 
    local (xo,fo)=optimise(fun3b,θG3b,lb3b,ub3b)
    nrange[:,i]=xo[:] # store the nuisance parameters
    llc0[i]=fo[1] # store the loglikelihood 
end
nllc0=llc0.-maximum(llc0); # calculate normalised loglikelihood
println("Profile c0 in confidence set Complete-----------------------------------------------------------------------------")
# Construct prediction interval.
@time lb,ub=predicrealc0(c0range,nrange,t)
println("profile realisation Prediction for c0 Complete---------------------------------------------------------------") 
lb = Float64.(lb); ub = Float64.(ub)



x_range = 1:1:200
x = 1:δ:200
s3b1 = plot(x,cmle.*J,color=:red,legend=false,lw=2,ls=:dash)
s3b1 = plot!(x_range, lb.*J, lw=0, fillrange=ub.*J, fillalpha=0.40, xlims=(1,200), ylims=(-5,25),
         color=:red, label=false, grid=true, tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-5,0,10,20,25], ["",latexstring("0"),"",latexstring("20")]),linecolor=:transparent)
s3b1 = hline!([0],lw=2,ls=:dash,color=:black)
s3b1 = hline!([20],lw=2,ls=:dash,color=:black)
s3b1 = scatter!(1:step:200,data,color=:red,markerstrokecolor=:red,ms=1.5)
s3b1 = annotate!(100, -13.7,text(L"i", :left, tfs))
s3b1 = annotate!(-100*0.6, 10, text(L"C^{\mathrm{o}} ", :left, tfs,rotation=90))
#= s3b1 = annotate!(-200*0.4, 28, text(L"(\mathrm{h})", :left, tfs)) =#



s3b2 = plot(x_range, lb.*J .- fcmle(x_range).*J, lw=0, fillrange=ub.*J .- fcmle(x_range).*J, fillalpha=0.40,
         color=:red, label=false, grid=true, xlims=(1,200), ylims=(-10,10), tickfontsize=tfs, framestyle=:box,
         xticks=([1,50,100,150,200], [latexstring("1"),"",latexstring("100"),"",latexstring("200")]),
         yticks=([-10,-5,0,5,10], [latexstring("-10"),"",latexstring("0"),"",latexstring("10")]),linecolor=:transparent)
s3b2 = hline!([0],lw=2,ls=:dash,color=:black,legend=false)
s3b2 = annotate!(100, -16,text(L"i", :left, tfs))
s3b2 = annotate!(-100*0.6, -7, text(L"C_{0.9}^{\mathrm{o}} - \hat{C} ", :left, tfs,rotation=90))
#= s3b2 = annotate!(-200*0.4, 12, text(L"(\mathrm{i})", :left, tfs)) =#


s3 = plot(s3a,s3b1,s3b2,layout = (1,3), link = :y, 
             bottom_margin = 6mm,top_margin=6mm, right_margin = 2mm,left_margin = 5mm,size=(580,150))  # Here we set the top margin
savefig(s3,"$path\\c0.pdf")
display(s3)
