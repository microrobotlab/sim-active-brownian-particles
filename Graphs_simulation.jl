
include("MyActiveBrownianParticle.jl");
using BenchmarkTools
using Plots 
using Distances
using NaNStatistics

N = 10000
Delta_t = 1e-2
t_tot= N*Delta_t

tauMax = t_tot/10               # è l'effettivo delta t massimo su cui posso calcolare l'MSD
N_Max = Int64(tauMax/Delta_t)   # è il numero massimo di frame della videocamera. Per me che sto simulando e basta sarà dato dal delta_t_MAX/delta_t_MIN su cui posso calcolare l'MSD

# ---------------------------------------------------------------
# Analisi ActiveParticle1 con R = 0.5µm e v = 0µm/s 
plot1, matrMSD1 = traj_and_MSD(0.0, 0.0, 0.5, 0.0, 100, N, Delta_t, N_Max) 
plot(plot1)

graphMSD1 = plot();

xMSD1 = Array(0:Delta_t:tauMax)
yMSD1 = vec(nanmean(matrMSD1, dims=2))
std_MSD1 = vec(nanstd(matrMSD1, dims=2))

plot!(xMSD1, yMSD1, ribbon=std_MSD1, fillalpha=.1, marker=false, legend=false);
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")

# ---------------------------------------------------------------
# Analisi ActiveParticle1 con R = 0.5µm e v = 5µm/s 
plot2, matrMSD2 = traj_and_MSD(0.0, 0.0, 0.5, 5.0, 100, N, Delta_t, N_Max) 
plot(plot2)

graphMSD2 = plot();

xMSD2 = Array(0:Delta_t:tauMax)
yMSD2 = vec(nanmean(matrMSD2, dims=2))
std_MSD2 = vec(nanstd(matrMSD2, dims=2))

plot!(xMSD2, yMSD2, ribbon=std_MSD2, fillalpha=.1, marker=false, legend=false);
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")

# ---------------------------------------------------------------
# Realizzo subplot delle due tipologie di ActiveParticle1 delle traiettorie
plot(plot1, plot2; layout = 2, plot_title = "ActiveParticle1 (R = 0.5µm)", title = ["v = 0µm/s" "v = 5µm/s"], plot_titlevspan = 0.1)

# Plotto in un unico grafico l'MSD di 1 per le due velocità 
MSD = plot();
plot!(xMSD1, yMSD1, ribbon=std_MSD1, fillalpha=.1, marker=false, legend=false);
plot!(xMSD2, yMSD2, ribbon=std_MSD2, fillalpha=.1, marker=false, legend=false);
plot(MSD)

#PROBLEMA: l'ultimo elemento della matrMSD2 è molto maggiore di quello del matrMSD1 (circa 200 contro circa 17)
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Analisi ActiveParticle2 con R = 1.5µm e v = 0µm/s 
plot3, matrMSD3 = traj_and_MSD(0.0, 0.0, 1.5, 0.0, 100, N, Delta_t, N_Max) 
plot(plot3)

graphMSD3 = plot();

xMSD3 = Array(0:Delta_t:tauMax)
yMSD3 = vec(nanmean(matrMSD3, dims=2))
std_MSD3 = vec(nanstd(matrMSD3, dims=2))

plot!(xMSD3, yMSD3, ribbon=std_MSD3, fillalpha=.1, marker=false, legend=false);
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")

# ---------------------------------------------------------------
# Analisi ActiveParticle1 con R = 0.5µm e v = 5µm/s 
plot4, matrMSD4 = traj_and_MSD(0.0, 0.0, 1.5, 5.0, 100, N, Delta_t, N_Max) 
plot(plot4)

graphMSD4 = plot();

xMSD4 = Array(0:Delta_t:tauMax)
yMSD4 = vec(nanmean(matrMSD4, dims=2))
std_MSD4 = vec(nanstd(matrMSD4, dims=2))

plot!(xMSD4, yMSD4, ribbon=std_MSD4, fillalpha=.1, marker=false, legend=false);
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")

# ---------------------------------------------------------------
# Realizzo subplot delle due tipologie di ActiveParticle1 delle traiettorie
plot(plot3, plot4; layout = 2, plot_title = "ActiveParticle2 (R = 1.5µm)", title = ["v = 0µm/s" "v = 5µm/s"], plot_titlevspan = 0.1)

# Plotto in un unico grafico l'MSD di 1 per le due velocità 
MSD = plot();
plot!(xMSD3, yMSD3, ribbon=std_MSD3, fillalpha=.1, marker=false, legend=false);
plot!(xMSD4, yMSD4, ribbon=std_MSD4, fillalpha=.1, marker=false, legend=false);
plot(MSD)

#PROBLEMA: l'ultimo elemento della matrMSD2 è molto maggiore di quello del matrMSD1 (circa 2200 contro circa 6)
# ---------------------------------------------------------------