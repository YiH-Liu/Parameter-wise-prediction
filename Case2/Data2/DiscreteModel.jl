using Random
using CSV
using DataFrames
using StatsBase
using FilePathsBase
using LaTeXStrings
using Plots

# In this script, we generate the data required to create Figure 3, Figure S2 and Figure S3.
path = dirname(@__FILE__)

# Init(): Set up the Initial condition for stochstic model.
function init(W,H,Δ,r_0,c_0)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # W: (number) The width of the lattice.
    # H: (number) The height of the lattice.
    # Δ: (number) The lattice size.
    # r_0: (number) Initial radius of circular barrier assay.
    # c_0: (number) Initial cells density.
    # Output:
    # A: (matrix) Lattice in matrix form, A[j,i] is associate with site (J+1-j,i) of the lattice.
    # N: (number) Total number of agents.
    #------------------------------------------------------------------------------------------------------------------
    # Genrate the empty lattice based on W and H.
    I = Int(W/Δ + 1); J = Int(H/Δ + 1); 
    A = zeros(J,I);
    # Find the center of the lattice (Cartesian coordinate).
    x_ctr = W/2; y_ctr = H/2
    # Find the site in the circle center at A[j_ctr,i_ctr] with radius r_0.
    site_available = []
    for i = 1:I
        for j = 1:J
            # Find the corresponding Cartesian coordinate.
            x_i = (i-1)*Δ; y_i = (J - j)*Δ
            # Find the distance from site to the center of circle.
            r_i = sqrt(((x_i-x_ctr)^2) + ((y_i-y_ctr)^2))
            if r_i <= r_0
                if rand() <= c_0
                    A[j,i] = 1
                end
            end
        end
    end

    N = sum(A)
    return A,N
end 

# realisation(): update the lattice matrix A from time t to time t + τ.
function realisation(A,P_m,k,I,J,N)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # A: (matrix) Lattice in matrix form, A[j,i] is associate with site (J+1-j,i) of the lattice.
    # P_m: (number) Probability of move.
    # I: (number) Number of colomns in the lattice.
    # J: (number) Number of rows in the lattice.
    # N: (number) Totol number of cells on the lattice.
    # Output:
    # Ac: (matrix) Updated lattice in matrix form, A[j,i] is associate with site (J+1-j,i) of the lattice.
    # Nc: (number) Total number of agents.
    #------------------------------------------------------------------------------------------------------------------
    # Number of select cells
    n_m = 0;
    # Random sequential update- select until we have total N cells.
    while n_m < N
         # Select a random ith colomn and jth row of the lattice matrix.
         i = rand(1:I); j = rand(1:J)
         if A[j,i] == 1
            # A non empty site is selected. 
            n_m = n_m + 1
            if rand() <= P_m
                # If the cell move decide the direction of movement.
                prob = rand()
                if prob <= 1/4   # move up
                    if j!=1 # not at boundary
                        if A[j-1,i] == 0 
                            A[j-1,i] = 1
                            A[j,i] = 0
                        end
                    end 
                elseif prob <= 1/2 # move down
                    if j!=J # not at boundary
                        if A[j+1,i] == 0 
                            A[j+1,i] = 1
                            A[j,i] = 0
                        end
                    end
                elseif prob <= 3/4 # move right
                    if i != I # not at boundary
                        if A[j,i+1] == 0
                            A[j,i+1] = 1
                            A[j,i] = 0
                        end
                    end
                else # move left
                    if i != 1 # not at boundary
                        if A[j,i-1] == 0
                            A[j,i-1] = 1
                            A[j,i] = 0
                        end
                    end
                end
            end
         end
    end
    n_p = 0
    Ac = A
    while n_p < N
        # Select a random ith colomn and jth row of the lattice matrix.
        i = rand(1:I); j = rand(1:J)
        if A[j,i] == 1
            n_p = n_p + 1
            if rand() <= k
                # If the agent proliferate decide the direction of proliferation.
                prob = rand()
                if prob <= 1/4   # move up
                    if j!=1 # not at boundary
                        if Ac[j-1,i] == 0 
                            Ac[j-1,i] = 1
                        end
                    end 
                elseif prob <= 1/2 # move down
                    if j!=J # not at boundary
                        if Ac[j+1,i] == 0 
                            Ac[j+1,i] = 1
                        end
                    end
                elseif prob <= 3/4 # move right
                    if i != I # not at boundary
                        if Ac[j,i+1] == 0
                            Ac[j,i+1] = 1
                        end
                    end
                else # move left
                    if i != 1 # not at boundary
                        if Ac[j,i-1] == 0
                            Ac[j,i-1] = 1
                        end
                    end
                end
            end
        end
    end
    Nc = sum(Ac)
    return Ac,Nc
end


function discrete_simulation(A,N,P_m,k,W,H,Δ,τ,t)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # A: (matrix) Lattice in matrix form at initial condition, A[j,i] is associate with site (J+1-j,i) of the lattice.
    # N: (number) Totol number of agents on the lattice.
    # P_m: (number) Probability of move.
    # k: (number) Probability of proliferate.
    # W: (number) The width of the lattice.
    # H: (number) The height of the lattice.
    # Δ: (number) The lattice site size.
    # τ: (number) The discrete time step duration.
    # t: (number) Time of Simulation.
    # Output:
    # A: (matrix) Updated lattice in matrix form at time t, A[j,i] is associate with site (J+1-j,i) of the lattice.
    #------------------------------------------------------------------------------------------------------------------
    I = Int(W/Δ + 1); J = Int(H/Δ + 1); 
    for i = τ:τ:t
        A,N = realisation(A,P_m,k,I,J,N)
    end 
    return A
end

# ctrfer(): Transform jth row and ith colomn of the lattice matrix A into lattice site (i,j).
function ctrfer(i,j,J)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # i: (number) i-th colomn of the lattice matrix A.
    # j: (number) j-th row of the lattice matrix A.
    # J: (number) Number of rows in the lattice.
    # Output:
    # (site_i,site_j): (tuple) Lattice site (site_i,site_j) correspond to jth row and ith colomn of the lattice matrix A.
    #------------------------------------------------------------------------------------------------------------------
    site_i = i
    site_j = J+1 - j
    return (site_i,site_j)
end

# Indices(): Extract the Indices of each agent in lattice matrix A and store them in dataframe.
function Indices(A)  
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # A: (matrix) lattice in matrix form, A[j,i] is associate with site (J+1-j,i) of the lattice.
    # Output:
    # index: (dataframe) index.i stores all the i-indices of agents on the lattice and the corresponding j-indices are stored in index.j.
    #------------------------------------------------------------------------------------------------------------------
    (J,I) = size(A)
    #Indices for subpopulation 1.
    indices = []
    # Extract the indices of each agent on lattice and store them in as list of tuple.
    for j in 1:J
        for i in 1:I
            if A[j,i] >0.025
                site = ctrfer(i,j,J)
                push!(indices, site)
            end
        end
    end
    # i j index for subpopulation 1.
    i_indices = [index[1] for index in indices]
    j_indices = [index[2] for index in indices]
    # Store in a dataframe
    index = DataFrame(i = vec(i_indices),j = vec(j_indices))
    return index
end



# Lattice size.
W = 199; H = 199; Δ = 1
# Discrete model parameters.
P_m = 1; k =0.01
# Simuation time.
t = 600; τ = 1

# Inital condition
#---------------------------
r_0 = 20; c_0 = 0.8


# Genrate the initial condition.
A0,N0= init(W,H,Δ,r_0,c_0)
# Extract indices at initial condition.
index_0 = Indices(A0) 
count0 = sum(A0[91:110,:], dims=1)
# save the indices at initial condition for snapshot.
CSV.write("$path\\index_0.csv", index_0)

# Simulate the discrete model at time t.
@time A = discrete_simulation(A0,N0,P_m,k,W,H,Δ,τ,t)

# Extract indices at time t.
index = Indices(A) 

# Save the indices at time t.
CSV.write("$path\\index.csv", index)

# Save the count data at initial condition
count1 = sum(A[91:110,:], dims=1)
countdf0 = DataFrame(a = vec(count0))
countdf = DataFrame(a = vec(count1))
CSV.write("$path\\data0.csv", countdf0)
CSV.write("$path\\data.csv", countdf)


# Simulate the discrete model at time t for m times and average it.
(r,c)=size(A0)
m = 100
global A_m = zeros(r,c)
global A0_m = zeros(r,c)
for i = 1:m
    local A0,N0= init(W,H,Δ,r_0,c_0)
    global A0_m = A0_m .+ A0
    local A = discrete_simulation(A0,N0,P_m,k,W,H,Δ,τ,t)
    global A_m = A_m .+ A
    display(i)
end
A0_average = A0_m./m
A_average = A_m./m
function densR(A,m,Δ,R,bn)
    bs = R/bn
    (J,I) = size(A)
    # Find the center of the lattice (Cartesian coordinate).
    x_ctr = (I-1)/2; y_ctr = (J-1)/2

    dens = zeros(bn,1)
    for i = 1:I
        for j = 1:J
            if A[j,i] != 0
                # Find the corresponding Cartesian coordinate.
                x_i = (i-1)*Δ; y_i = (J - j)*Δ
                # Find the distance from site to the center of circle.
                r_i = sqrt(((x_i-x_ctr)^2) + ((y_i-y_ctr)^2))
                if r_i <= R
                    index = Int(floor((r_i/bs) + 1))
                    dens[index] = dens[index] +(A[j,i]./m)
                end
            end
        end
    end
    Area = zeros(Int(round((R/bs))),1)
    N = length(Area)
    r = 0:bs:R
    for i = 1:N
        Area[i] = (pi*(r[i+1]/Δ)^2) - (pi*(r[i]/Δ)^2)
    end
    
    dens = (dens./Area)
    r = bs:bs:R
    return r,dens
end
# Find the density along radius
R = 99.5;bn = 21
r,dens = densR(A_m,m,Δ,R,bn)
DiscreteDens = DataFrame(r = vec(r), u = vec(dens))
CSV.write("$path\\DiscreteDens.csv", DiscreteDens)



# store averaged density in matrix form
# Reverse each column of the matrix
reversed_matrix0 = reverse(A0_average, dims=1)
reversed_matrix  = reverse(A_average, dims=1)
# Define the x and y coordinates
x = 1:size(reversed_matrix, 2)
y = 1:size(reversed_matrix, 1)

# Convert the reversed matrix to a DataFrame
df0 = DataFrame(reversed_matrix0,:auto)
df = DataFrame(reversed_matrix,:auto)
# Save the DataFrame to a CSV file
CSV.write("$path\\A0Discrete.csv", df0)
CSV.write("$path\\ADiscrete.csv", df)