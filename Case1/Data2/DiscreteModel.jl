using Random
using CSV
using DataFrames
using StatsBase
using FilePathsBase
using LaTeXStrings
using Plots
# In this script, we generate the data required to create Supplementary Information (SI) Figure S4.

# Get the path of the current script.
path = dirname(@__FILE__)

# Init(): Set up the Initial condition for stochstic model.
function init(W,H,Δ,h_0,d1_0,d2_0)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # W: (number) The width of the lattice.
    # H: (number) The height of the lattice.
    # Δ: (number) The lattice size.
    # h_0: (number) Initial width of occupancy.
    # d1_0: (number) Initial density for agent from subpopulation 1.
    # d2_0: (number) Initial density for agent from subpopulation 2.
    # Output:
    # A1: (matrix) Lattice in matrix form contain occupancy information for agent from subpopulation 1, A1[j,i] is associate with site (J+1-j,i) of the lattice.
    # A2: (matrix) Lattice in matrix form contain occupancy information for agent from subpopulation 2, A2[j,i] is associate with site (J+1-j,i) of the lattice.
    # N1: (number) Total number of agent from subpopulation 1.
    # N2: (number) Total number of agent from subpopulation 2.
    #------------------------------------------------------------------------------------------------------------------
    # Genrate the empty lattice based on W and H.
    I = Int(W/Δ + 1); J = Int(H/Δ + 1); 
    A1 = zeros(J,I);  A2 = zeros(J,I);
    # Find the center of the lattice (site index).
    i_ctr = Int(round((I+1)/2));
    # Find the center of the lattice (Cartesian coordinate).
    x_ctr = W/2
    # Number of cells in each colomn.
    N1_0 = Int(round(d1_0*J)); N2_0 = Int(round(d2_0*J))
    # Total number of cells.
    N1 = 0; N2 =0
    # Find the colomn in the initial condition.
    for i = 1:I
        # Find the corresponding Cartesian coordinate.
        x_i = (i-1)*Δ
        # Find the distance from site to the center of circle.
        h_i = abs(x_i-x_ctr)
        if h_i <= h_0
            # occupy the colomn with given initial density
            index1 = sample(1:J, N1_0, replace=true)
            index2 = sample(1:J, N2_0, replace=true)
            for j in index1
                A1[Int(j),i] = A1[Int(j),i] + 1
                N1 = N1 + 1
            end

            for j in index2
                A2[Int(j),i] = A2[Int(j),i] + 1
                N2 = N2 + 1
            end
        end
        
    end
    return A1,A2,N1,N2
end 
    
# realisation(): update the lattice matrix A from time t to time t + τ.
function realisation(A1,A2,P_m,ρ,Rg,Rr,I,J,N1,N2)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # A1: (matrix) Lattice in matrix form contain occupancy information for agent from subpopulation 1, A1[j,i] is associate with site (J+1-j,i) of the lattice.
    # A2: (matrix) Lattice in matrix form contain occupancy information for agent from subpopulation 2, A2[j,i] is associate with site (J+1-j,i) of the lattice.
    # P_m: (number) Probability of move.
    # ρ: (number) Bias parameter.
    # Rg (number) Red to green transition rate.
    # Rr (number) Green to red transition rate.
    # I: (number) Number of colomns in the lattice.
    # J: (number) Number of rows in the lattice.
    # N1: (number) Total number of agent from subpopulation 1.
    # N2: (number) Total number of agent from subpopulation 2.
    # Output:
    # A1_c: (matrix) Lattice in matrix form contain occupancy information for agent from subpopulation 1 after 1 realisation, A1[j,i] is associate with site (J+1-j,i) of the lattice.
    # A2_c: (matrix) Lattice in matrix form contain occupancy information for agent from subpopulation 2 after 1 realisation, A2[j,i] is associate with site (J+1-j,i) of the lattice.
    # N1_c: (number) Total number of agent from subpopulation 1 after 1 realisation.
    # N2_c: (number) Total number of agent from subpopulation 2 after 1 realisation.
    #------------------------------------------------------------------------------------------------------------------
    # Total number of cell
    N = N1 + N2
    # Number of select cells
    n1_m = 0; n1_g = 0; n1_r = 0
    n2_m = 0; n2_g = 0; n2_r = 0
    # Probability to move to the right for subpopulations 1.
    prob_right = (1 + ρ) / 4
    # Random sequential update- select until we have total N cells.
    while n1_m < N1
         # Select a random ith colomn and jth row of the lattice matrix.
         i = rand(1:I); j = rand(1:J)
         if A1[j,i] != 0
            # A1 non empty site is selected. 
            n1_m = n1_m + A1[j,i]
            for k = 1:A1[j,i] 
                if rand() <= P_m
                    # If the cell move decide the direction of movement.
                    prob = rand()
                    if prob <= 1/4   # move up
                        if j!=1 # not at boundary
                            A1[j-1,i] = A1[j-1,i] + 1
                            A1[j,i] = A1[j,i] - 1
                        else
                            A1[J,i] = A1[J,i] + 1
                            A1[j,i] = A1[j,i] - 1
                        end
                    elseif prob <= 1/2 # move down
                        if j!=J # not at boundary
                            A1[j+1,i] = A1[j+1,i] + 1
                            A1[j,i] = A1[j,i] - 1
                        else
                            A1[1,i] = A1[1,i] + 1
                            A1[j,i] = A1[j,i] - 1
                        end
                    elseif prob <= 1/2 +  prob_right# move right
                        if i != I # not at boundary  
                            A1[j,i+1] = A1[j,i+1] + 1
                            A1[j,i] = A1[j,i] - 1
                        end
                    else # move left
                        if i != 1 # not at boundary      
                            A1[j,i-1] = A1[j,i-1] + 1
                            A1[j,i] = A1[j,i] - 1
                        end
                    end
                end
            end
         end
    end
    while n2_m < N2
        # Select a random ith colomn and jth row of the lattice matrix.
        i = rand(1:I); j = rand(1:J)
        if A2[j,i] != 0
           # A2 non empty site is selected. 
           n2_m = n2_m + A2[j,i]
           for k = 1:A2[j,i] 
               if rand() <= P_m
                   # If the cell move decide the direction of movement.
                   prob = rand()
                   if prob <= 1/4   # move up
                       if j!=1 # not at boundary
                           A2[j-1,i] = A2[j-1,i] + 1
                           A2[j,i] = A2[j,i] - 1
                       else
                           A2[J,i] = A2[J,i] + 1
                           A2[j,i] = A2[j,i] - 1
                       end
                   elseif prob <= 1/2 # move down
                       if j!=J # not at boundary
                           A2[j+1,i] = A2[j+1,i] + 1
                           A2[j,i] = A2[j,i] - 1
                       else
                           A2[1,i] = A2[1,i] + 1
                           A2[j,i] = A2[j,i] - 1
                       end
                   elseif prob <= 1/2 +  prob_right# move right
                       if i != I # not at boundary  
                           A2[j,i+1] = A2[j,i+1] + 1
                           A2[j,i] = A2[j,i] - 1
                       end
                   else # move left
                       if i != 1 # not at boundary      
                           A2[j,i-1] = A2[j,i-1] + 1
                           A2[j,i] = A2[j,i] - 1
                       end
                   end
               end
           end
        end
    end
    N1_c = N1; N2_c = N2; 
    A1_c = A1; A2_c = A2;
    while n1_g < N1 
        # Select a random ith colomn and jth row of the lattice matrix.
        i = rand(1:I); j = rand(1:J)
        if A1[j,i] != 0
           # A non empty site is selected. 
           n1_g = n1_g + A1[j,i]
           for k = 1:A1[j,i] 
            if rand() <= Rg
                # If the cell die remove the cell.
                A1_c[j,i] = A1_c[j,i] - 1
                N1_c = N1_c - 1

                A2_c[j,i] = A2_c[j,i] + 1
                N2_c = N2_c + 1
            end
           end
        end
    end

    while n2_r < N2
        # Select a random ith colomn and jth row of the lattice matrix.
        i = rand(1:I); j = rand(1:J)
        if A2[j,i] != 0
           # A non empty site is selected. 
           n2_r = n2_r + A2[j,i]
           for k = 1:A2[j,i] 
            if rand() <= Rr
                # If the cell die remove the cell.
                A1_c[j,i] = A1_c[j,i] + 2
                N1_c = N1_c + 2

                A2_c[j,i] = A2_c[j,i] - 1
                N2_c = N2_c - 1
            end
           end
        end
    end
    return A1_c,A2_c,N1_c,N2_c
end




function discrete_simulation(A1,A2,N1,N2,P_m,ρ,Rg,Rr,W,H,Δ,τ,t)
    #------------------------------------------------------------------------------------------------------------------
    # Input:
    # A1: (matrix) Lattice in matrix form contain occupancy information for agent from subpopulation 1, A1[j,i] is associate with site (J+1-j,i) of the lattice.
    # A2: (matrix) Lattice in matrix form contain occupancy information for agent from subpopulation 2, A2[j,i] is associate with site (J+1-j,i) of the lattice.
    # P_m: (number) Probability of move.
    # ρ: (number) Bias parameter.
    # Rg (number) Red to green transition rate.
    # Rr (number) Green to red transition rate.
    # N1: (number) Total number of agent from subpopulation 1.
    # N2: (number) Total number of agent from subpopulation 2.
    # W: (number) The width of the lattice.
    # H: (number) The height of the lattice.
    # Δ: (number) The lattice site size.
    # τ: (number) The discrete time step duration.
    # t: (number) Time of Simulation.
    # Output:
    # A1: (matrix) Updated lattice in matrix form contain occupancy information for agent from subpopulation 1 at time t, A1[j,i] is associate with site (J+1-j,i) of the lattice.
    # A2: (matrix) Updated lattice in matrix form contain occupancy information for agent from subpopulation 2 at time t, A2[j,i] is associate with site (J+1-j,i) of the lattice.
    #------------------------------------------------------------------------------------------------------------------   
    I = Int(W/Δ + 1); J = Int(H/Δ + 1); 
    for i = τ:τ:t
        A1,A2,N1,N2 = realisation(A1,A2,P_m,ρ,Rg,Rr,I,J,N1,N2)
    end 
    return A1,A2
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
            for n = 1:A[j,i]
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
W = 199; H = 19; Δ = 1
# Time step duration.

τ= 1
# Discrete model parameters.
P_m = 1; Rg = 0.02;
Rr = 0.03; ρ=0.5
D = P_m/4
# Simuation time.
t=200
# Inital condition
#---------------------------
h_0 = 20; d1_0 = 0.5; d2_0 = 0.2;


# Genrate the initial condition.
A1_0,A2_0,N1_0,N2_0= init(W,H,Δ,h_0,d1_0,d2_0)

# Extract indices at initial condition.
index1_0 = Indices(A1_0) 
index2_0 = Indices(A2_0) 
# save the indices at initial condition for snapshot.
CSV.write("$path\\index1_0.csv", index1_0)
CSV.write("$path\\index2_0.csv", index2_0)
# initial count data and save it
count1_0 = sum(A1_0, dims=1)
count2_0 = sum(A2_0, dims=1)
count_0 = DataFrame(a = vec(count1_0),b = vec(count2_0))
CSV.write("$path\\data_0.csv", count_0)

# model simulation
A1_0,A2_0,N1_0,N2_0 = init(W,H,Δ,h_0,d1_0,d2_0)
@time A1,A2 = discrete_simulation(A1_0,A2_0,N1_0,N2_0,P_m,ρ,Rg,Rr,W,H,Δ,τ,t)
# Extract indices at time t.
index1 = Indices(A1) 
index2 = Indices(A2) 
# Save the indices at time t.
CSV.write("$path\\index1.csv", index1)
CSV.write("$path\\index2.csv", index2)

# count data and save it
count1 = sum(A1, dims=1)
count2 = sum(A2, dims=1)
count = DataFrame(a = vec(count1),b = vec(count2))
CSV.write("$path\\data.csv", count)

