# PURPOSE: Output of ABP main 
# all codes in repository are complied in this function
#VARIABLES: Destination folder path and filename
include("ABP main.jl")
include("ABP file.jl")
include("ABP analysis.jl")
include("ABP SD.jl")

using BenchmarkTools,Plots,Distances,NaNStatistics,CSV, DataFrames

gr()

N = 10000
Delta_t = 1e-2
t_tot= N*Delta_t

tauMax = t_tot/10               #is the actual maximum delta t over which I can calculate the MSD
N_Max = Int64(tauMax/Delta_t)   # is the maximum number of frames of the camera. For me that I am simulating and that's it it will be given by the delta_t_MAX / delta_t_MIN on which I can calculate the MSD

# -----------------------------------------------------------------------------------------------------------------------------------------------------
# THIS IS THE CODE TO CALL MAIN FUNCTION
# We plot the set of particles considering the correction of hard spheres

L = 100.0 	# μm box length
R = 2.0		# μm particle radius
v = 10.0 	# μm/s particle velocity
a=L/2
b=L/4
#pf_factor = (R^2)/(a*b)
pf_factor = (R^2)
DT, DR = diffusion_coeff(R).*[1e12, 1]
packing_fraction = 0.1
Np = round(Int,packing_fraction*L^2/(2R^2))  #Np is the number of particles in my set and I choose it random?
#π
Nt = 10000# Nt is the number of steps 
#println(" Number of particles: $Np") 
@timed graph = multiparticleE(Np,L,R,v,Nt);

println("multiparticleE complied\n")

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# THIS IS THE CODE TO CALL WALL FUNCTIONS IN THE MAIN FUNCTION

# Same with the wall condition (particles bounce off the edge)

@timed graph_wall = multiparticleE_wall(Np,L,R,v,Nt) # has values of x and y posiiton in each frame in graph_wall[1]

println("multiparticleE_wall complied\n")
#-------------------------------------------------------------------------------------------------------------------
# destination folder/filename selection
path="C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\Coding\\2023\\"  # destination folder path

filename="data_ellipse"   # filename base for all data

pathf= path*filename  
#---------------------------------------------------------------------------------------------------------------------
# file store
file_store(graph_wall,Nt,pathf)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# analysis
inside_Np=stat_analysis1(a,b,R,pathf)
# mean and standard deviation

analysis_SD= stat_analysis2(a,b,R,pathf)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# making animation
anim = @animate for i = 1:100:Nt
    scatter(graph_wall[1][i][:,1], graph_wall[1][i][:,2], aspect_ratio=:equal, lims=(-L/2, L/2),markersize=350R/L, background_color_inside=:Blues_6,marker =:circle,legend=false, title = "$(inside_Np) particles, steps $i, ellipse a=L/2, b= L/4")
    plot!(L/2*cos.(-π:0.01:π), L/4*sin.(-π:0.01:π))
    quiver!(graph_wall[1][i][:,1],graph_wall[1][i][:,2],quiver=(4*cos.(graph_wall[2][i,1]),4*sin.(graph_wall[2][i,1])), color=:red)
end
#marker_z=graph_wall[2][i,1], color=:rainbow, for 

f1= pathf*".gif"
gif(anim, f1)



