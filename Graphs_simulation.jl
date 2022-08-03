#Script Lapo Corti per simulazione movimento particella attiva
#OBIETTIVO: realizzare 4 grafici diversi (ognuno per una particella 
# di dimensione nota e velocità iniziale nota) e all'interno di ciascun 
# grafico vado a plottare 4-5 traiettorie diverse della stessa particella. 


include("MyActiveBrownianParticle.jl");

using BenchmarkTools
using Plots 
using Distances 

#Definizione parametri di interesse per la generazione delle traiettorie 
N = 10000
Delta_t = 1e-2


## ActiveParticle1 con R=0.5um e v=0um/s sarà chiamata abp1_1##

#inizializzazione particella (R=0.5e-6m , v=0m/s) 
orientazione1 = rand()*2*pi 
orientazione1 = round(orientazione1, digits=3)
abp1_1 = initABP( (0.0, 0.0, orientazione1), 0.5, 0.0 );


#creo il primo grafico grafico 
p, t = trajectory( abp1_1, N, Delta_t);   #NOTA: p è un vettore di tuple, ovvero è un vettore con N elementi e ciascun elemento è una tupla (x,y, teta)

#creo i vettori degli assi 
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]

#NOTA: questo comando mi permette di creare un vettore i cui elementi 
# sono dati dai primi elementi di pi dove pi 
# assume tutti i valori dentro p (cioè parte da p)

#plotto in x,y
plot1 = plot(x,y, range=[-100,75],  title = "ActiveParticle1 (R=0.5µm, v=0µm/s)", aspect_ratio= :equal, legend=false) 
#scatter!([x[end]],[y[end]], legend=false)
xlabel!("x [μm]")
ylabel!("y [μm]")

#realizzo più traiettorie sullo stesso plot
orientazione2 = rand()*2*pi 
orientazione2 = round(orientazione2, digits=3)
abp1_1 = initABP( (0.0, 0.0, orientazione2), 0.5, 0.0 );
p, t = trajectory( abp1_1, N, Delta_t);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)

orientazione3 = rand()*2*pi 
orientazione3 = round(orientazione3, digits=3)
abp1_1 = initABP( (0.0, 0.0, orientazione3), 0.5, 0.0 );
p, t = trajectory( abp1_1, N, Delta_t);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75], legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)

orientazione4 = rand()*2*pi 
orientazione4 = round(orientazione4, digits=3)
abp1_1 = initABP( (0.0, 0.0, orientazione4), 0.5, 0.0 );
p, t = trajectory( abp1_1, N, Delta_t);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)
    
#subplot e metto i 4 grafici in un unico grafico facendo subplot
#plot(plot1, plot2 , plot3, plot4; layout = 4, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4" ], plot_title = "ActiveParticle1 (R=0.5µm, v=0µm/s)", plot_titlevspan=0.1)


## ActiveParticle1 con R=0.5um e v=5um/s sarà chiamata abp1_2##

#inizializzazione particella (R=0.5e-6m , v=5µm/s)
orientazione5 = rand()*2*pi 
orientazione5 = round(orientazione5, digits=3)
abp1_2 = initABP( (0.0, 0.0, orientazione5), 0.5, 5.0 );

#creo il primo grafico grafico 
p, t = trajectory( abp1_2, N, Delta_t);

#NOTA: p è un vettore di tuple, ovvero è un vettore con N elementi 
# e ciascun elemento è una tupla (x,y, teta)

#creo i vettori degli assi 
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]

#NOTA: questo comando mi permette di creare un vettore i cui elementi 
# sono dati dai primi elementi di pi dove pi 
# assume tutti i valori dentro p (cioè parte da p)

#plotto in x y
plot2 = plot(x,y,range=[-100,75], title = "ActiveParticle2 (R=0.5µm, v=5µm/s)", aspect_ratio= :equal,legend=false) 
#scatter!([x[end]],[y[end]], legend=false)
xlabel!("x [μm]")
ylabel!("y [μm]")

#realizzo più grafici nello stesso plot (subplot) 
# (vedere se riesco a migliorare questa parte)
orientazione6 = rand()*2*pi 
orientazione6 = round(orientazione6, digits=3)
abp1_2 = initABP( (0.0, 0.0, orientazione6), 0.5, 5.0 );
p, t = trajectory( abp1_2, N, Delta_t);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],aspect_ratio= :equal,legend=false)
#scatter!([x[end]],[y[end]], legend=false)


orientazione7 = rand()*2*pi 
orientazione7 = round(orientazione7, digits=3)
abp1_2 = initABP( (0.0, 0.0, orientazione7), 0.5, 5.0 );
p, t = trajectory( abp1_2, N, Delta_t);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],aspect_ratio= :equal,legend=false)
#scatter!([x[end]],[y[end]], legend=false)


orientazione8 = rand()*2*pi 
orientazione8 = round(orientazione8, digits=3)
abp1_2 = initABP( (0.0, 0.0, orientazione8), 0.5, 5.0 );
p, t = trajectory( abp1_2, N, Delta_t);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],aspect_ratio= :equal,legend=false)
#scatter!([x[end]],[y[end]], legend=false)

    
#subplot e metto i 4 grafici in un unico grafico
#plot(plot1, plot2 , plot3, plot4; layout = 4, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4" ], plot_title = "ActiveParticle2 (R=0.5µm, v=5µm/s)", plot_titlevspan=0.1)

# subplot deei due grafici della stessa particella:
plot(plot1, plot2; layout = 2, plot_title = "ActiveParticle1 (R = 0.5µm)", title = ["v = 0µm/s" "v = 5µm/s"], plot_titlevspan = 0.1)










## ActiveParticle2 con R=1.5µm E v=0µm/s ##

#inizializzazione particella (R=1.5e-6m , v=0µm/s)
orientazione1 = rand()*2*pi 
orientazione1 = round(orientazione1, digits=3)
abp2_1 = initABP( (0.0, 0.0, orientazione1), 1.5, 0.0 );

#creo il primo grafico grafico 
p, t = trajectory( abp2_1, N, Delta_t); #NOTA: p è un vettore di tuple, ovvero è un vettore con N elementi e ciascun elemento è una tupla (x,y, teta)
 
#creo i vettori degli assi 
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]

#NOTA: questo comando mi permette di creare un vettore i cui elementi 
# sono dati dai primi elementi di pi dove pi 
# assume tutti i valori dentro p (cioè parte da p)

#plotto in x y
plot3 = plot(x,y, range=[-100,75], title = "ActiveParticle3 (R=1.5µm, v=0µm/s)", aspect_ratio= :equal) 
#scatter!([x[end]],[y[end]], legend=false)
xlabel!("x [μm]")
ylabel!("y [μm]")

#realizzo più grafici nello stesso plot (subplot) 
# (vedere se riesco a migliorare questa parte)
orientazione2 = rand()*2*pi 
orientazione2 = round(orientazione2, digits=3)
abp2_1 = initABP( (0.0, 0.0, orientazione2), 1.5, 0.0 );
p, t = trajectory( abp2_1, N, Delta_t);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)

orientazione3 = rand()*2*pi 
orientazione3 = round(orientazione3, digits=3)
abp2_1 = initABP( (0.0, 0.0, orientazione3), 1.5, 0.0 );
p, t = trajectory( abp2_1, N, Delta_t);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)

orientazione4 = rand()*2*pi 
orientazione4 = round(orientazione4, digits=3)
abp2_1 = initABP( (0.0, 0.0, orientazione4), 1.5, 0.0 );
p, t = trajectory( abp2_1, N, Delta_t);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)
    
#subplot e metto i 4 grafici in un unico grafico
#plot(plot1, plot2 , plot3, plot4; layout = 4, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4"], plot_title = "ActiveParticle3 (R=1.5µm, v=0µm/s)", plot_titlevspan=0.1)



## ActiveParticle2 con R=1.5µm E v=5µm/s ##

#inizializzazione particella (R=1.5e-6m , v=5µm/s)
orientazione5 = rand()*2*pi 
orientazione5 = round(orientazione5, digits=3)
abp2_2 = initABP( (0.0, 0.0, orientazione5), 1.5, 5.0 );

#creo il primo grafico grafico 
p, t = trajectory( abp2_2, N, Delta_t);

#creo i vettori degli assi 
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]

#NOTA: questo comando mi permette di creare un vettore i cui elementi 
# sono dati dai primi elementi di pi dove pi 
# assume tutti i valori dentro p (cioè parte da p)

#plotto in x y
plot4 = plot(x,y,range=[-100,75], title = "ActiveParticle4 (R=1.5µm, v=5µm/s)", aspect_ratio= :equal) 
#scatter!([x[end]],[y[end]], legend=false)
xlabel!("x [μm]")
ylabel!("y [μm]")

#realizzo più grafici nello stesso plot (subplot) 
# (vedere se riesco a migliorare questa parte)
orientazione6 = rand()*2*pi 
orientazione6 = round(orientazione6, digits=3)
abp2_2 = initABP( (0.0, 0.0, orientazione6), 1.5, 5.0 );
p, t = trajectory( abp2_2, N, Delta_t);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)

orientazione7 = rand()*2*pi 
orientazione7 = round(orientazione7, digits=3)
abp2_2 = initABP( (0.0, 0.0, orientazione7), 1.5, 5.0 );
p, t = trajectory( abp2_2, N, Delta_t);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)

orientazione8 = rand()*2*pi 
orientazione8 = round(orientazione8, digits=3)
abp2_2 = initABP( (0.0, 0.0, orientazione8), 1.5, 5.0 );
p, t = trajectory( abp2_2, N, Delta_t);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)
    
#subplot e metto i 4 grafici in un unico grafico
#plot(plot1, plot2 , plot3, plot4; layout = 4, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4"], plot_title = "ActiveParticle4 (R=1.5µm, v=5µm/s)", plot_titlevspan=0.1)

#subplot per visualizzare plot3 e plot4 vicini
plot(plot3, plot4; layout = 2, plot_title = "ActiveParticle2 (R = 1.5µm)", title = ["v = 0µm/s" "v = 5µm/s"], plot_titlevspan = 0.1)






## ESTRAZIONE FUNZIONE MSD DAL FILE DI GAIA ##
#------------------------------
tauMax=20000     ##da inserire. E' il numero massimo di frame
tr = 1/abp23.DR  ## da calcolare a seconda della particella analizzata. Mi serve per trovare il valore del tempo che discrimina il passaggio da andamento parabolico a lineare dell'MSD
#------------------------------
#Funzione calcolo MSD
function MSDcalculate(x,y,tauMax)
    ltrack= N
    msd=zeros(tauMax+1)
    for tau in 1:tauMax
        for i in tau+1:ltrack
            msd[tau+1]+=((x[i]-x[i-tau])^2+(y[i]-y[i-tau])^2)/(ltrack-tau)
        end
    end
    return msd
end
#--------------------------------
## plot MSD
graph=plot();

yMSD = MSDcalculate(x,y, tauMax)
xMSD = Array(0:1e-3:tauMax*1e-3) 

plot!(xMSD,yMSD, marker=true,legend=false);
xlabel!("Δt [s]");
ylabel!("MSD [μm²]")

display(graph)
