# Define an "ActiveBrownianParticle" Type
abstract type ActiveBrownianParticle end

# Define a specific type for 2D ABPs
struct ABP2 <: ActiveBrownianParticle
	x::Float64 		# x position (μm)
	y::Float64 		# y position (μm)
	θ::Float64  	# orientation (rad)
	R::Float64 		# Radius (μm)
	v::Float64 		# velocity (μm/s)
	DT::Float64 	# translational diffusion coefficient (μm^2/s)
	DR::Float64 	# rotational diffusion coefficient (rad^2/s)
end

## Get position and orientation of the particle or ensemble (CURRENTLY ONLY 2D)
position(abp::ABP2) = ( abp.x, abp.y, abp.θ )

#position(abpe::ABPE2) = [ abpe.x abpe.y ]
#orientation(abpe::ABPE2) = abpe.θ

## Initialize ABP particle (CURRENTLY ONLY 2D)
function initABP(position::NTuple, R::Float64, v::Float64; T::Float64=300.0, η::Float64=1e-3)
    # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
    DT, DR = diffusion_coeff(1e-6R)

    if length(position) == 3
        abp = ABP2( position..., R, v, 1e12DT, DR)
    else
        println("No init method available")
    end

    return abp
end

## Calculate the step position (and orientation) for the particle (CURRENTLY ONLY 2D)
function step(abp::ABP, δt::Float64) where {ABP <: ActiveBrownianParticle}
    if length(position(abp)) == 3
        δx = sqrt(2*abp.DT*δt)*randn() + abp.v*δt*cos(abp.θ)
        δy = sqrt(2*abp.DT*δt)*randn() + abp.v*δt*sin(abp.θ)
        δθ = sqrt(2*abp.DR*δt)*randn()
        δp = ( δx, δy, δθ )
    else
        println("No step method available")
    end
    return δp
end

## Update position and orientation of the particle 
update(abp::ABP, step) where {ABP <: ActiveBrownianParticle} = ABP( (position(abp) .+ step)..., abp.R, abp.v, abp.DT, abp.DR )

##Calculate diffusion coefficient
function diffusion_coeff(R::Float64, T::Float64=300.0, η::Float64=1e-3)
    # Boltzmann constant [J/K]
    kB = 1.38e-23
    # friction coefficient [Ns/m]
    γ = 6*pi*R*η
    # translational diffusion coefficient [m^2/s]
    DT = kB*T/γ
    # rotational diffusion coefficient [rad^2/s]
    DR = 6*DT/(8*R^2)
    return DT, DR
end;

## Create single particle trajectory
function trajectory(abp::ABP, N, δt::Float64) where {ABP <: ActiveBrownianParticle}
	p = [position(abp)]
	t = 0:δt:δt*N
    
	for i in 1:N
		abp = update(abp, step(abp,δt))
		push!(p, position(abp))
	end
	return p, t
end

## Plot trajectory
function plot_trajectory(x,y)
    plot(x, y, legend=false, aspect_ratio=:equal);
	scatter!([x[end]],[y[end]], legend=false);
    xlabel!("x [μm]");
    ylabel!("y [μm]");
end


function animate_trajectory(p, t, time_vec, tidx)
	push!(time_vec, time())
	if length(time_vec)>1
		framerate = round(1/(time_vec[end] - time_vec[end-1]),digits=1)
	else
		framerate = 0.0
	end
	nframe = length(time_vec)
	
	ti = tidx[1:nframe]
	x = first.(p)[ti]
	y = [ pi[2] for pi in p[ti] ]
	
	plot_trajectory(x, y)
	title!("t = $(round(t[ti[end]],digits=1)) s, $framerate fps")
end

#------------------------------
# MSD calculate (Estratta da file MSD di Gaia e riscritta per semplificarla)
function MSDcalculate(x,y,N_Max,N)
    ltrack= N+1
    msd=zeros(N_Max+1)
    for tau in 1:N_Max
        for i in tau+1:ltrack
            msd[tau+1]+=((x[i]-x[i-tau])^2+(y[i]-y[i-tau])^2)/(ltrack-tau)
        end
    end
    return msd
end

#--------------------------------
# Funzione che genera n traiettorie, mi plotta solo le prime 4 in un grafico e fa la media degli MSD di tutte 
function traj_and_MSD(x0, y0, R::Float64, v::Float64, num_traj::Int64, N, Delta_t::Float64, N_Max)
    graph = plot();
    matrMSD = fill(NaN, N_Max+1, num_traj)

    for i in 1:num_traj
        orientazione = rand()*2*pi
        orientazione = round(orientazione, digits=3)
        abp = initABP( (x0, y0, orientazione), R, v);

        p, t = trajectory( abp, N, Delta_t);

        x = [pi[1] for pi in p]
        y = [pi[2] for pi in p]

        if i <= 4 
            plot!(x,y, range=[-50,50],  title = "ActiveParticle (R=$R µm, v=$v µm/s)", aspect_ratio= :equal, legend=false)
            xlabel!("x [μm]")
            ylabel!("y [μm]")
        end

        matrMSD[1:N_Max+1, i] = MSDcalculate(x,y, N_Max, N)
    end

    return graph, matrMSD 

end

#------------------------------ 
#------------------------------ 
#------------------------------ 
#------------------------------ 

# INIZIALIZZAZIONE DI UN INSIEME DI PARTICELLE IN CUI CONSIDERIAMO LE INTERAZIONI STERICHE

# Define an "ABPsEnsemble" Type
abstract type ABPsEnsemble end

# Define a specific type for 2D ABPsEnsembles (CURRENTLY ASSUMING ALL PARTICLES ARE EQUAL)
struct ABPE2 <: ABPsEnsemble
    Np::Int64                      # number of particles
    L::Float64                      # size of observation space (μm)
	R::Float64  # Radius (μm)                                   --> Vector{Float64}(undef,Np)
	v::Float64 	# velocity (μm/s)                               --> Vector{Float64}(undef,Np)
	DT::Float64 # translational diffusion coefficient (μm^2/s)  --> Vector{Float64}(undef,Np)
	DR::Float64 # rotational diffusion coefficient (rad^2/s)    --> Vector{Float64}(undef,Np)
	x::Vector{Float64}    # x position (μm)
	y::Vector{Float64}    # y position (μm)
	θ::Vector{Float64}    # orientation (rad)
end

## Initialize ABP ensemble (CURRENTLY ONLY 2D)
function initABPE(Np::Int64, L::Float64, R::Float64, v::Float64; T::Float64=300.0, η::Float64=1e-3)
    # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
    DT, DR = diffusion_coeff(1e-6R)

    # ONLY 2D!
    xyθ = (rand(Np,3).-0.5).*repeat([L L 2π],Np)
    xyθ[:,1:2], dists, superpose, uptriang = hardsphere(xyθ[:,1:2],R)
    abpe = ABPE2( Np, L, R, v, 1e12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3])

    return abpe, (dists, superpose, uptriang)
end

#------------------------------ 
# FUNZIONI PER LA CORREZIONE DELLE POSIZIONI DI PARTICELLE SOVRAPPOSTE
function hardsphere_correction!(xy::Array{Float64,2}, dists::Array{Float64,2}, superpose::BitArray{2}, R::Float64; tol::Float64=1e-3)
    Np = size(superpose,1) #mi dà il numero di righe di superpose 
    for np1 in 1:Np
        if any(superpose[np1,:])  #se almeno un valore della riga np1 è true  
            np2 = findfirst(superpose[np1,:])
            Δp = (xy[np1,:] - xy[np2,:]) .* ( ( (1+tol)*2R / dists[np1,np2] - 1 ) / 2 )
            xy[np1,:] += Δp
            xy[np2,:] -= Δp
            dists[np2,np2+1:Np] = pairwise(Euclidean(), xy[np2:np2,:], xy[np2+1:Np,:], dims=1 )  #????
            superpose[np2,np2+1:Np] = (dists[np2,np2+1:Np] .< 2R*(1-tol))  #????
        end
    end
    return nothing
end

function hardsphere!(xy::Array{Float64,2}, dists::Array{Float64,2}, superpose::BitArray{2}, uptriang::BitArray{2}, R::Float64; tol::Float64=1e-3)
    superpositions = 1
    counter = 0
    # @time begin
    while superpositions > 0
        dists .= pairwise(Euclidean(),xy,xy,dims=1)
        superpose .= (dists .< 2R*(1-tol)).*uptriang
        # @show(findall(superpose))
        superpositions = sum(superpose)
        # @show(superpositions)
        if superpositions > 0
            hardsphere_correction!(xy,dists,superpose,R,tol=tol)
        end
        counter += 1
        # @show(counter)
        if counter >= 100
            println("$superpositions superpositions remaining after 100 cycles")
            break
        end
    end
    # end
    return nothing
end

function hardsphere(xy::Array{Float64,2}, R::Float64; tol::Float64=1e-3)
    Np = size(xy,1)
    dists = zeros(Np,Np)
    superpose = falses(Np,Np)
    uptriang = falses(Np,Np)
    for i = 1:Np-1
        uptriang[i,i+1:Np] .= true
    end
    hardsphere!(xy, dists, superpose, uptriang, R; tol=tol)
    return xy, dists, superpose, uptriang
end
#------------------------------ 

function multiparticleE(Np::Integer, L::Float64, R::Float64, v::Float64, Nt::Int64=2, δt::Float64=1e-3)
    (Nt isa Int64) ? Nt : Nt=convert(Int64,Nt)
    
    ABPE = Vector{ABPE2}(undef,Nt+1)
    ABPE[1], matrices = initABPE( Np, L, R, v ) # including initial hardsphere correction
    
    simulate!(ABPE, matrices, Nt, δt)

    return position.(ABPE)
end

function simulate!(ABPE, matrices, Nt, δt)
    # PΘ = [ (position(abpe), orientation(abpe)) ]
    # pθ = PΘ[1]
    
    for nt in 1:Nt
        ABPE[nt+1] = update(ABPE[nt],matrices,δt)
        println("Step $nt")
    end
    return nothing
end

position(abpe::ABPE2) = [ abpe.x abpe.y ]
orientation(abpe::ABPE2) = abpe.θ

function update(abpe::ABPE, matrices::Tuple{Matrix{Float64}, BitMatrix, BitMatrix}, δt::Float64) where {ABPE <: ABPsEnsemble}
    pθ = ( position(abpe), orientation(abpe) ) .+ step(abpe,δt)

    periodic_BC_array!(pθ[1],abpe.L, abpe.R)

    hardsphere!(pθ[1], matrices[1], matrices[2], matrices[3], abpe.R)
    # @btime hardsphere!($p[:,1:2], $matrices[1], $matrices[2], $matrices[3], $params.R)
    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.v, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )

    return new_abpe
end


function step(abpe::ABPE, δt::Float64) where {ABPE <: ABPsEnsemble}
    
    if size(position(abpe),2) == 2
        δp = sqrt.(2*δt*abpe.DT)*randn(abpe.Np,2) .+ abpe.v*δt*[cos.(abpe.θ) sin.(abpe.θ)]
        δθ = sqrt(2*abpe.DR*δt)*randn(abpe.Np)
    else
        println("No step method available")
    end
    #if nt == 1 
        #println("lo step vero di questo giro è: $δp") 
    #end
    return (δp, δθ)
end

function periodic_BC_array!(xy::Array{Float64,2},L::Float64, R)   #quando una particella supera un bordo ricompare dalla parte opposta
	# Boundary conditions: horizontal edge
	idx = abs.(xy[:,1]) .> L/2 + R #creo vettore idx in cui ho 1 dove il valore assoluto della coordinata x delle particelle è fuori dall'area di osservazione 
	if any(idx)
		xy[idx,1] .-= sign.(xy[idx,1]).*L   #dove ho uni in idx faccio ricomparire la particella dalla parte opposta di x rispetto allo 0 
	end
	# Boundary conditions: vertical edge
	idy = abs.(xy[:,2]) .> L/2 + R
	if any(idy)
		xy[idy,2] .-= sign.(xy[idy,2]).*L
	end
	return nothing
end
# ---------------------------------------------------
# impongo la condizione di rimbalzo delle particelle

function wall_condition!(xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) #quando una particella tocca un bordo rimbalza
    #if nt == 1 
        #println("lo step che uso poi in questo giro è: $step_mem") 
    #end
	# Boundary conditions: horizontal edge
	idx = abs.(xy[:,1]) .> (L/2 - R)
	if any(idx)
		xy[idx,1] .-= 2*sign.(xy[idx,1]).*(abs.(xy[idx,1]) .- (L/2 - R)) 
	end
	# Boundary conditions: vertical edge
	idy = abs.(xy[:,2]) .> (L/2 - R)
	if any(idy)
        xy[idy,2] .-= 2*sign.(xy[idy,2]).*(abs.(xy[idy,2]) .- (L/2 - R))
	end
	return nothing
end

function update_wall(abpe::ABPE, matrices::Tuple{Matrix{Float64}, BitMatrix, BitMatrix}, δt::Float64) where {ABPE <: ABPsEnsemble}
    memory_step = step(abpe,δt)
    #if nt == 1 
        #println("lo step di questo giro è: $memory") 
    #end
    pθ = ( position(abpe), orientation(abpe) ) .+ memory_step

    wall_condition!(pθ[1],abpe.L, abpe.R, memory_step[1])

    hardsphere!(pθ[1], matrices[1], matrices[2], matrices[3], abpe.R)
    # @btime hardsphere!($p[:,1:2], $matrices[1], $matrices[2], $matrices[3], $params.R)
    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.v, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )

    return new_abpe
end

function simulate_wall!(ABPE, matrices, Nt, δt)
    # PΘ = [ (position(abpe), orientation(abpe)) ]
    # pθ = PΘ[1]
    
    for nt in 1:Nt
        ABPE[nt+1] = update_wall(ABPE[nt],matrices,δt)
        println("Step $nt")
    end
    return nothing
end

function multiparticleE_wall(Np::Integer, L::Float64, R::Float64, v::Float64, Nt::Int64=2, δt::Float64=1e-3)
    (Nt isa Int64) ? Nt : Nt=convert(Int64,Nt)
    
    ABPE = Vector{ABPE2}(undef,Nt+1)
    ABPE[1], matrices = initABPE( Np, L, R, v ) # including initial hardsphere correction
    
    simulate_wall!(ABPE, matrices, Nt, δt)
    
    return position.(ABPE)
end