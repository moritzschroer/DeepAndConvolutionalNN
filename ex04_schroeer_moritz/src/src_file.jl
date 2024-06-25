module WolfModule
using ProgressMeter
import Statistics: mean, std
import Pkg; 
Pkg.add("StatsBase")
using StatsBase: autocor

function flipRandomCluster(σ, J, flip_prob, M, E)
    nx, ny, nz = size(σ)
    i = rand(1:nx)
    j = rand(1:ny)
    k = rand(1:nz)
    visited =  fill(false, nx, ny, nz)
    σ, visited, M, E = add_neighbors_to_cluster(σ, J, i, j, k, flip_prob, nx, ny, nz, visited, M, E)
   
    return σ, M, E 
end

function add_neighbors_to_cluster(σ, J, i, j, k, flip_prob, nx, ny, nz, visited, M, E)
    #println("Hallo")
    #println(σ)
    visited[i,j,k] = true
    
    σ_i_j_k = σ[i,j,k]
    M -= 2*σ_i_j_k
    
    σ[i,j,k] = -σ_i_j_k
    E += 2*J*sumNearestNeighbors(σ, i, j, k, nx, ny, nz)
    
    if rand() < flip_prob 
        if i == nx  
            if visited[1,j,k] != true && σ_i_j_k == σ[1,j,k]
                σ, visited, M, E = add_neighbors_to_cluster(σ, J, 1, j, k, flip_prob, nx, ny, nz, visited, M, E)
            end
        elseif visited[i+1,j,k] != true && σ_i_j_k == σ[i+1,j,k]
            σ, visited, M, E = add_neighbors_to_cluster(σ, J, i+1, j, k, flip_prob, nx, ny, nz, visited, M, E)
        end
    end
    if rand() < flip_prob 
        if i == 1 
            if visited[nx,j,k] != true && σ_i_j_k == σ[nx,j,k]
                σ, visited, M, E = add_neighbors_to_cluster(σ, J, nx, j, k, flip_prob, nx, ny, nz, visited, M, E)
            end 
        elseif visited[i-1,j,k] != true && σ_i_j_k == σ[i-1,j,k]
            σ, visited, M, E = add_neighbors_to_cluster(σ, J, i-1, j, k, flip_prob, nx, ny, nz, visited, M, E)
        end
    end

    if rand() < flip_prob 
        if j == 1 
            if visited[i,ny,k] != true && σ_i_j_k == σ[i,ny,k]
                σ, visited, M, E = add_neighbors_to_cluster(σ, J, i, ny, k, flip_prob, nx, ny, nz, visited, M, E)
            end
        elseif visited[i,j-1,k] != true && σ_i_j_k == σ[i,j-1,k]
            σ, visited, M, E = add_neighbors_to_cluster(σ, J, i, j-1, k, flip_prob, nx, ny, nz, visited, M, E)
        end
    end
    if rand() < flip_prob 
        if j == ny 
            if visited[i,1,k] != true && σ_i_j_k == σ[i,1,k]
                σ, visited, M, E = add_neighbors_to_cluster(σ, J, i, 1, k, flip_prob, nx, ny, nz, visited, M, E)
            end
        elseif visited[i,j+1,k] != true && σ_i_j_k == σ[i,j+1,k]
            σ, visited, M, E = add_neighbors_to_cluster(σ, J, i, j+1, k, flip_prob, nx, ny, nz, visited, M, E)
        end
    end

    if rand() < flip_prob 
        if k == 1 
            if visited[i,j,nz] != true && σ_i_j_k == σ[i,j,nz]
                σ, visited, M, E = add_neighbors_to_cluster(σ, J, i, j, nz, flip_prob, nx, ny, nz, visited, M, E)
            end
        elseif visited[i,j,k-1] != true && σ_i_j_k == σ[i,j,k-1]
            σ, visited, M, E = add_neighbors_to_cluster(σ, J, i, j, k-1, flip_prob, nx, ny, nz, visited, M, E)
        end
    end
    if rand() < flip_prob 
        if k == nz 
            if visited[i,j,1] != true && σ_i_j_k == σ[i,j,1]
                σ, visited, M, E = add_neighbors_to_cluster(σ, J, i, j, 1, flip_prob, nx, ny, nz, visited, M, E)
            end
        elseif visited[i,j,k+1] != true && σ_i_j_k == σ[i,j,k+1]
            σ, visited, M, E = add_neighbors_to_cluster(σ, J, i, j, k+1, flip_prob, nx, ny, nz, visited, M, E)
        end
    end
    return σ, visited, M, E 
end


function init_random(nx, ny, nz)
    σ = rand([-1,1], nx, ny, nz)
    return σ
end


function sumNearestNeighbors(σ,i,j,k,nx,ny,nz)
    sum_nn = 0
    #nearest neighbors in x direction
    if i == nx
        sum_nn += σ[1,j,k]*σ[i,j,k]
    else
        sum_nn += σ[i+1,j,k]*σ[i,j,k]
    end
    if i == 1
        sum_nn += σ[nx,j,k]*σ[i,j,k]
    else
        sum_nn += σ[i-1,j,k]*σ[i,j,k]
    end

    #nearest neighbors in y direction
    if j == ny
        sum_nn += σ[i,1,k]*σ[i,j,k]
    else
        sum_nn += σ[i,j+1,k]*σ[i,j,k]
    end
    if j == 1
        sum_nn += σ[i,ny,k]*σ[i,j,k]
    else
        sum_nn += σ[i,j-1,k]*σ[i,j,k]
    end

    #nearest neighbors in z direction
    if k == nz
        sum_nn += σ[i,j,1]*σ[i,j,k]
    else
        sum_nn += σ[i,j,k+1]*σ[i,j,k]
    end
    if k == 1
        sum_nn += σ[i,j,nz]*σ[i,j,k]
    else
        sum_nn += σ[i,j,k-1]*σ[i,j,k]
    end
    
    return sum_nn
end

function sumSpinsAndNNInteraction(σ,nx,ny,nz)
    sumSpins, interactionSum = 0, 0   
     for i in 1:nx
            for j in 1:ny
                for k in 1:nz
                    sumSpins += σ[i,j,k] 
                    interactionSum += sumNearestNeighbors(σ,i,j,k,nx,ny,nz)/2
                end
            end
        end
    return sumSpins, interactionSum
end

function simulateWolf(Temperatures, J::Float64, nx::Int64, ny::Int64, nz::Int64, thermalization_steps::Int64, sample_steps::Int64)
    
    t_steps = length(Temperatures)
    E_list, M_list, T_list = zeros(t_steps), zeros(t_steps), zeros(t_steps)
    M_std, E_std, susceptibility, heatcapacity = zeros(t_steps), zeros(t_steps), zeros(t_steps), zeros(t_steps)
    binderCumulant, autocorM, autocorE = zeros(t_steps), zeros(sample_steps, t_steps), zeros(sample_steps, t_steps)
    runtimes = zeros(t_steps)
    for j in 1:t_steps
        t = @elapsed begin
            Energies, Magnetizations = zeros(sample_steps), zeros(sample_steps)
            T = Temperatures[j]
            σ = init_random(nx,ny,nz)

            sumSpins, interactionSum = sumSpinsAndNNInteraction(σ,nx,ny,nz)
            E = -J*interactionSum 
            M = sumSpins
            flip_prob = 1 - exp(-2/T*J)
            
            @showprogress for i in 1:thermalization_steps
                σ, M, E = flipRandomCluster(σ, J, flip_prob, M, E)
            end
            @showprogress for i in 1:sample_steps
                σ, M, E = flipRandomCluster(σ, J, flip_prob, M, E)
                Energies[i] = E
                Magnetizations[i] = abs(M)/(nx*ny*nz)
            end
            M_list[j]  = mean(Magnetizations)
            E_list[j]         = mean(Energies)
            M_std[j]          = std(Magnetizations)
            E_std[j]          = std(Energies)
            susceptibility[j] = nx*ny*nz*1/T*M_std[j]^2
            heatcapacity[j]   = (1/T)^2*E_std[j]^2
            binderCumulant[j] = 1 - mean(Magnetizations.^4)/(3*mean(Magnetizations.^2)^2)
            autocorM[:, j] = autocor(Magnetizations, 0:sample_steps-1)
            autocorE[:, j] = autocor(Energies, 0:sample_steps-1)
        end
        runtimes[j] = t 
    end
    return Temperatures, E_list, M_list, susceptibility, heatcapacity, binderCumulant, autocorM, autocorE, runtimes

end

#sort unsorted list to sorted list, sorted from lowest to highest temperature
function sort(tempNotSorted, EList, MList, susceptibility, heat_capacity, binder_cumulant, AutocorrM_list, AutocorrE_list, runtimes)
    n = length(tempNotSorted);
    indices = [];
    for i in 1:n
        push!(indices,i);
    end
    for i in 1:n-1
        min_index = i;
        for j in i+1:n
            if tempNotSorted[indices[j]]<= tempNotSorted[indices[min_index]]
                min_index = j;
            end
        end
        indices[i], indices[min_index] = indices[min_index], indices[i];
    end
    M_sorted, E_sorted, susceptibility_sorted, heat_capacity_sorted, temp_sort, binder_cumulant_sort, AutocorrM_list_sorted, AutocorrE_list_sorted, runtimes_sorted = [], [], [], [], [], [], [], [], []
    for i in 1:n
        push!(E_sorted, EList[indices[i]])
        push!(M_sorted, MList[indices[i]])
        push!(susceptibility_sorted, susceptibility[indices[i]])
        push!(heat_capacity_sorted, heat_capacity[indices[i]])
        push!(binder_cumulant_sort, binder_cumulant[indices[i]])
        push!(AutocorrM_list_sorted, AutocorrM_list[indices[i]])
        push!(AutocorrE_list_sorted, AutocorrE_list[indices[i]])
        push!(runtimes_sorted, runtimes[indices[i]])
        
    end
    return E_sorted, M_sorted, susceptibility_sorted, heat_capacity_sorted, binder_cumulant_sort, AutocorrM_list_sorted, AutocorrE_list_sorted, runtimes_sorted    
end

#initilaize random spins
function initializeRandom(σ)
    n_x, n_y, n_z = size(σ)
    for i in 1:n_x
        for j in 1:n_y
            for k in 1:n_z
                σ[i,j,k] = rand([0,1])
            end
        end
    end
    return σ
end

#initialize all spins in spin up state
function initializeSpinsUp(σ)
    n_x, n_y, n_z = size(σ)
    for i in 1:n_x
        for j in 1:n_y
            for k in 1:n_z
                σ[i,j,k] = 1
            end
        end
    end
    return σ
end

#calculate the energy exponentials and return the values in an array
function initializeEnergieExponentials(β, J, H)
    #initialize energies such that for the field term the case that
    #σ_i=-1->σ_j=1 , i.e. dE = 2*H is stored in expValsExtFieldTerm[5]
    #and dE = -2*H is stored in expValsExtFieldTerm[1], because that makes
    #it possible to acess the values of the exponential without additional 
    #if statement in the method accept (compared to the case where expValsExtFieldTerm
    #is just an array with len = 2
    #the interaction term exponentials are stored in an array where the case, that
    #dE = -J*12 is stored in expValsInteractionTerm[1], dE = -J*11 is stored in 
    #expValsInteractionTerm[2] ...
    expValsExtFieldTerm, expValsInteractionTerm = zeros(5), zeros(25)

    expValsInteractionTerm[13] = 1
    for i in 1:12
        expValsInteractionTerm[i] = exp(-β*J*(13-i))
        expValsInteractionTerm[i+13] = exp(β*J*i)
    end
  
    expValsExtFieldTerm[1] = exp(-β*H*2)
    expValsExtFieldTerm[5] = exp(β*H*2)
    return expValsExtFieldTerm, expValsInteractionTerm
end

#sum over the nearest neighbours of a a spin
function sumNearestNeighbors(σ,i,j,k,nx,ny,nz)
    sum_nn = 0
    #nearest neighbors in x direction
    if i == nx
        sum_nn += σ[1,j,k]*σ[i,j,k]
    else
        sum_nn += σ[i+1,j,k]*σ[i,j,k]
    end
    if i == 1
        sum_nn += σ[nx,j,k]*σ[i,j,k]
    else
        sum_nn += σ[i-1,j,k]*σ[i,j,k]
    end

    #nearest neighbors in y direction
    if j == ny
        sum_nn += σ[i,1,k]*σ[i,j,k]
    else
        sum_nn += σ[i,j+1,k]*σ[i,j,k]
    end
    if j == 1
        sum_nn += σ[i,ny,k]*σ[i,j,k]
    else
        sum_nn += σ[i,j-1,k]*σ[i,j,k]
    end

    #nearest neighbors in z direction
    if k == nz
        sum_nn += σ[i,j,1]*σ[i,j,k]
    else
        sum_nn += σ[i,j,k+1]*σ[i,j,k]
    end
    if k == 1
        sum_nn += σ[i,j,nz]*σ[i,j,k]
    else
        sum_nn += σ[i,j,k-1]*σ[i,j,k]
    end
    
    return sum_nn
end

#sum over all spins  and over all nearest neighbours of all spins
function sumSpinsAndNNInteraction(σ,nx,ny,nz)
sumSpins, interactionSum = 0, 0   
 for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                sumSpins += σ[i,j,k] 
                interactionSum += sumNearestNeighbors(σ,i,j,k,nx,ny,nz)/2
            end
        end
    end
    return sumSpins, interactionSum
end

#now the Metropolis algorithm

#condition to accept new configuration
function acceptMetropolis(expValsExtFieldTerm, expValsInteractionTerm, dSumSpins, dSumNNSpins, dE)
    if dE <= 0.
        return true
    elseif rand() <= expValsExtFieldTerm[Int(dSumSpins)+3]*expValsInteractionTerm[Int(dSumNNSpins)+13]
        return true
    else
        return false
    end
end

#thermalization
function thermalizeMetropolis(σ, J, H, thermalization_steps, nx, ny, nz, expValsExtFieldTerm, expValsInteractionTerm)
    for i in 1:thermalization_steps
        i,j,k = rand(1:nx), rand(1:ny), rand(1:nz)
        σ[i,j,k] = (-1)*σ[i,j,k]

        dSumNNSpins = 2*sumNearestNeighbors(σ, i, j, k, nx, ny, nz)
        dSumSpins = 2*σ[i,j,k]
         
        dE = -J*dSumNNSpins - H*dSumSpins

        if  acceptMetropolis(expValsExtFieldTerm, expValsInteractionTerm, dSumSpins, dSumNNSpins, dE) == false
             σ[i,j,k] = (-1)*σ[i,j,k]
        end
    end
    return σ
end

#metropolis algorithm
function metropolis(σ, β, J, H, MCSweepSteps, MCSampleSteps, nx, ny, nz, expValsExtFieldTerm, expValsInteractionTerm)

    sumSpins, interactionSum = sumSpinsAndNNInteraction(σ,nx,ny,nz)
    N = nx*ny*nz
    E = (-J*interactionSum -H*sumSpins)
    M = sumSpins/N

    E_values, M_values = zeros(MCSampleSteps), zeros(MCSampleSteps)
    
    push!(E_values, E)
    push!(M_values, abs(M))
    
    for l in 1:MCSampleSteps
        for m in 1:MCSweepSteps
            i,j,k = rand(1:nx), rand(1:ny), rand(1:nz)
            σ[i,j,k] = (-1)*σ[i,j,k]

            dSumNNSpins = 2*sumNearestNeighbors(σ, i, j, k, nx, ny, nz)
            dSumSpins = 2*σ[i,j,k]  
            
            dE = -J*dSumNNSpins - H*dSumSpins

            if  acceptMetropolis(expValsExtFieldTerm, expValsInteractionTerm, dSumSpins, dSumNNSpins, dE) == true
                E += dE
                M += 2*σ[i,j,k]/(nx*ny*nz)
                
            else
                σ[i,j,k] = (-1)*σ[i,j,k]
            end
        end
        
        E_values[l] = E
        M_values[l] = abs(M)
    end

    E_mean = mean(E_values)
    M_mean = mean(M_values)
    autocorM = autocor(M_values, 0:MCSampleSteps-1)
    autocorE = autocor(E_values, 0:MCSampleSteps-1)
    
    M_2 = mean(M_values .^ 2)
    M_4 = mean(M_values .^ 4)
    M_std = std(M_values)
    E_std = std(E_values)
    
    χ = N*β*M_std^2
    hc = β^2*E_std^2
    U = 1 - M_4/(3*M_2^2)

    
    return E_mean, M_mean, χ, hc , σ, U, autocorM, autocorE
end

#simulation using metropolis algorithm
function simulateIsingMetropolis(temperature, J, H, MCSweepSteps, MCSampleSteps, thermalization_steps, nx, ny, nz)
    E_list, M_list, susceptibility, heat_capacity, temp_not_sorted, binder_cumulant_not_sorted, AutocorrM_list, AutocorrE_list = [], [], [], [], [], [], [], []
    temp_steps = length(temperature)
    runtimes = zeros(temp_steps)
    Threads.@threads for i in 1:temp_steps
        t = @elapsed begin
            β = 1/temperature[i]
            
            expValsExtFieldTerm, expValsInteractionTerm = initializeEnergieExponentials(β, J, H)
            
            σ = zeros(Int8, nx, ny, nz)
            σ  = initializeRandom(σ)
            σ  = initializeSpinsUp(σ)
            σ = thermalizeMetropolis(σ, J, H, thermalization_steps, nx, ny, nz, expValsExtFieldTerm, expValsInteractionTerm)
            E_value, M_value, χ, hc , σ, U, M_Autocorr, E_Autocorr  = metropolis(σ, β, J, H, MCSweepSteps, MCSampleSteps, nx, ny, nz, expValsExtFieldTerm, expValsInteractionTerm)
            push!(E_list, E_value)
            push!(M_list, M_value)
            push!(susceptibility, χ)
            push!(heat_capacity, hc)
            push!(binder_cumulant_not_sorted, U)
            push!(temp_not_sorted, temperature[i])
            push!(AutocorrM_list, M_Autocorr)
            push!(AutocorrE_list, E_Autocorr)

        end
        runtimes[i] = t  
    end

    E_sorted, M_sorted, susceptibility_sorted, heat_capacity_sorted, binder_cumulant_sorted, AutocorrM_sorted, AutocorrE_sorted, runtimes_sorted = sort(temp_not_sorted, E_list, M_list, susceptibility, heat_capacity, binder_cumulant_not_sorted, AutocorrM_list, AutocorrE_list, runtimes)
    autocorM_return, autocorE_return = zeros(MCSampleSteps, length(temperature)), zeros(MCSampleSteps, length(temperature))
    # return M corr as 2d Matrix
    for j in 1:length(AutocorrE_sorted)
        autocorE_return[:, j] = AutocorrE_sorted[j]
        autocorM_return[:, j] = AutocorrM_sorted[j] 
    end
    return temperature, E_sorted, M_sorted, susceptibility_sorted, heat_capacity_sorted, binder_cumulant_sorted, autocorM_return, autocorE_return, runtimes_sorted
end


function calculateDerivativeNumerically(x_data, y_data)
    
    n = length(y_data)
    derivative = zeros(n)

    # boundary values
    derivative[1] = (y_data[2] - y_data[1]) / (x_data[2] - x_data[1])
    derivative[n] = (y_data[n] - y_data[n-1]) / (x_data[n] - x_data[n-1])
    
    # average values for the other data points
    for i in 2:n-1
        derivative[i] = (y_data[i+1] - y_data[i-1]) / (x_data[i+1] - x_data[i-1])
    end

    return derivative
end

export simulateIsingCreutz, simulateIsingMetropolis, calculateDerivativeNumerically



export simulateWolf, simulateMC

end


