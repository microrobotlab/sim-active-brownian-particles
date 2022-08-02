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

## Update position and orientation of the particle (create new)
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
function trajectory(abp::ABP, N; δt::Float64=1e-3) where {ABP <: ActiveBrownianParticle}
	p = [position(abp)]
	t = 0:δt:δt*N
    
	for i in 1:N
		abp = update(abp, step(abp,δt))
		push!(p, position(abp))
	end
	return p, t
end
## INTERVALLO DI TEMPO TOTALE DENTRO CUI CONSIDERO LA TRAIETTORIA:
# δt = 1e-3
# N = 100000    (scelto da me in fase di scrittura del programma)
# t_tot = δt * N = 100s = 1min e 40sec

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