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

# In this script, We solve the PDE in Equaiton (2) and (3) snalytically and numerically at k = 200

# Get the path of the current script
path = dirname(@__FILE__)

# Analytic solution
function solution(x,t,D,v,k1,k2,c10,c20,h,L)
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
# diff!(): Discretized PDE 
function diff!(d_u, u, p, t)

    (D,v,k1,k2,δ) = p
    (S,N) = size(u)

    # Boundary at x = 0, associated with Equation (24).
    d_u[1,1] = (D/(δ^2))*(u[1,2]-u[1,1]) - k1*u[1,1] + 2*k2*u[2,1]
    d_u[2,1] = (D/(δ^2))*(u[2,2]-u[2,1]) + k1*u[1,1] - k2*u[2,1]

    # Associated with Equation (25).
    for n in 2:N-1
        d_u[1,n] = (D/(δ^2))*(u[1,n+1]-2*u[1,n]+u[1,n-1]) - (v/(2*δ))*(u[1,n+1]-u[1,n-1]) - k1*u[1,n] + 2*k2*u[2,n]
        d_u[2,n] = (D/(δ^2))*(u[2,n+1]-2*u[2,n]+u[2,n-1]) - (v/(2*δ))*(u[2,n+1]-u[2,n-1]) + k1*u[1,n] - k2*u[2,n]
    end
    
    # Boundary at x = 199, associated with Equation (26).
    d_u[1,N] =  (D/(δ^2))*(-u[1,N]+u[1,N-1]) - k1*u[1,N] + 2*k2*u[2,N]
    d_u[2,N] =  (D/(δ^2))*(-u[2,N]+u[2,N-1]) + k1*u[1,N] - k2*u[2,N]

    return d_u
end

# pdesolver(): Solves the ODE in diff!().
function pdesolver(time,D,v,k1,k2,h,c10,c20,δ,L)
    # ------------------------------------------------------------------------------------------------------------------
    # Input:
    # Time: (number) Solve ODE/PDE at t = Time.
    # D1: (number) Diffusivity for agents in subpopulation 1.
    # D2: (number) Diffusivity for agents in subpopulation 2.
    # v1: (number) Drift velocity for agents in subpopulation 1.
    # v2: (number) Drift velocity for agents in subpopulation 2.
    # δ: (number) Grid spacing.
    # L: (number) ODE/PDE is solved on 0 <= x <= L.
    # Output:
    # c1, c2: Solution of ODE/PDE at time t = Time.
    # ------------------------------------------------------------------------------------------------------------------
    # Construct initial condition for PDE.
    # Total number of mesh points N.
    N = Int(L/δ + 1)
    # Initial condition u0.
    c0 = zeros(2,N)
    # Update initial condition u0: c1(x) = 1 at 79 <= x <= 119.
    for n=1:N
       xx = (n-1)*δ
       if abs(xx-(L/2)) < h
          c0[1,n] = c10
          c0[2,n] = c20
       end
    end


    # Return the inital conditon if ODE/PDE is solved at t=0
    if time == 0
        c1_0 = vec(c0[1,:]); c2_0 = vec(c0[2,:])
        return c1_0,c2_0
    end

    # Solve the PDE using Heun's method at t=Time 
    p=(D,v,k1,k2,δ)
    tspan=(0,time)
    prob=ODEProblem(diff!,c0,tspan,p)
    alg=Heun() 
    sol=solve(prob,alg,saveat=time);
    sol = sol[end]
    c1 = sol[1,:];  c2= sol[2,:]
    c1 = vec(c1); c2 = vec(c2)
    return c1,c2
end 


D = 0.25;
v = 0.25
k1 = 0.02;
k2 = 0.03;
h = 20;
c10 = 0.5;
c20 = 0.2;
L = 199;
Time = 200
# test grid independent.
δ = 1
@time c1,c2=pdesolver(Time,D,v,k1,k2,h,c10,c20,δ,L)
f = plot(0:δ:L,c1, color=:red, ls=:dash,lw=6,label=L"δ = 1")
f = plot!(0:δ:L,c2, color=:green, ls=:dash,lw=6,label=L"δ = 1")

δ = 0.5
@time c1,c2=pdesolver(Time,D,v,k1,k2,h,c10,c20,δ,L)
f = plot!(0:δ:L,c1, color=:red,lw=3,label=L"δ = 0.5")
f = plot!(0:δ:L,c2, color=:green,lw=3,label=L"δ = 0.5")
display(f)

# Compare numerical and Analytic
x = 0:0.5:L
@time ac1,ac2 = solution(x,Time,D,v,k1,k2,c10,c20,h,L)
f1 = plot(0:0.5:L,ac1, color=:red, ls=:dash,lw=6,label="Analytic")
f1 = plot!(0:0.5:L,ac2, color=:green, ls=:dash,lw=6,label="Analytic")
f1 = plot!(0:δ:L,c1, color=:red,lw=3,label="Numerical")
f1 = plot!(0:δ:L,c2, color=:green,lw=3,label="Numerical")
display(f1)

sol = DataFrame(a = vec(ac1),b = vec(ac2))
CSV.write("$path\\solution.csv", sol)


ac10,ac20 = solution(x,0.00001,D,v,k1,k2,c10,c20,h,L)
sol = DataFrame(a = vec(ac10),b = vec(ac20))
CSV.write("$path\\solution0.csv", sol)