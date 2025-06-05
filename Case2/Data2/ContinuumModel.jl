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

# In this script, We solve the PDE in Equaiton (9) numerically at k = 100.

# Get the path of the current script
path = dirname(@__FILE__)



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
function pdesolver(t,D,k,r0,c_0,δ,R)
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

# Uxy(): convert solution to Cartesian coordinates.
function Uxy(I,J,u,r)
    ur = linear_interpolation(r,u)
    U = zeros(J,I)
    for i in 1:I
        for j in 1:J
            x_i = i-1
            y_j = J-j
            r_xy = sqrt((x_i-99.5)^2 + (y_j-99.5)^2)
            if r_xy < 99
                U[j,i] = ur(r_xy)
            end
        end
    end
    return U
end


D = 0.25;
k = 0.01;
r0 = 20;
c_0 = 0.8;

R = 99.5;
Time = 600
# test grid independent.
δ = 0.25
@time c=pdesolver(Time,D,k,r0,c_0,δ,R)
f = plot(0:δ:R,c, color=:red, ls=:dash,lw=6,label=L"δ = 0.25")
δ = 0.5
@time c0=pdesolver(0,D,k,r0,c_0,δ,R)
@time c=pdesolver(Time,D,k,r0,c_0,δ,R)
f = plot!(0:δ:R,c, color=:red,lw=3,label=L"δ = 0.5")
display(f)
# Save radial Density
ContinuumDens = DataFrame(r = vec(1:δ:R+1), u = vec(c))
CSV.write("$path\\ContinuumDens.csv", ContinuumDens)


#Save solution in Cartesian coordinates
I = 200; J =200;
A0Continuum = Uxy(I,J,c0,0:δ:R)
AContinuum = Uxy(I,J,c,0:δ:R)
# Reverse each column of the matrix
reversed_matrix0 =reverse(A0Continuum, dims=1)
reversed_matrix =reverse(AContinuum, dims=1) 
# Define the x and y coordinates
x = 1:size(reversed_matrix, 2)
y = 1:size(reversed_matrix, 1)

# Convert the reversed matrix to a DataFrame
df0 = DataFrame(reversed_matrix0,:auto)
df = DataFrame(reversed_matrix,:auto)

#Save solution at inital condition in Cartesian coordinates
CSV.write("$path\\A0Continuum.csv", df0)
#Save solution at k=600 condition in Cartesian coordinates
CSV.write("$path\\AContinuum.csv", df)