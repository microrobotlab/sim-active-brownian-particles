#Script Lapo Corti per simulazione movimento particella attiva
#OBIETTIVO: realizzare 4 grafici diversi (ognuno per una particella 
# di dimensione nota e velocità iniziale nota) e all'interno di ciascun 
# grafico vado a plottare 4-5 traiettorie diverse della stessa particella. 


include("MyActiveBrownianParticle.jl");

using BenchmarkTools
using Plots 
using Distances 

## PARTIAMO DAL PRIMO GRAFICO PER UNA PARTICELLA DI RAGGIO 0.5MICRON E V=0 ##

#inizializzazione particella (R=0.5e-6m , v=0m/s) 
orientazione = rand(float(0:pi/2))
abp1 = initABP( (0.0, 0.0, orientazione), 0.5, 0.0 );


#creo il primo grafico grafico 
N = 100000
p, t = trajectory( abp1, N);

#NOTA: p è un vettore di tuple, ovvero è un vettore con N elementi 
# e ciascun elemento è una tupla (x,y, teta)

#creo i vettori degli assi 
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]

#NOTA: questo comando mi permette di creare un vettore i cui elementi 
# sono dati dai primi elementi di pi dove pi 
# assume tutti i valori dentro p (cioè parte da p)

#plotto in x y

plot1 = plot(x,y, range=[-100,75], title = "ActiveParticle1 (R=0.5µm, v=0µm/s)", aspect_ratio= :equal,legend=false) 
#scatter!([x[end]],[y[end]], legend=false)
xlabel!("x [μm]")
ylabel!("y [μm]")


#realizzo più grafici nello stesso plot (subplot) 
# (vedere se riesco a migliorare questa parte)
p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75], legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)
    
#subplot e metto i 4 grafici in un unico grafico facendo subplot
#plot(plot1, plot2 , plot3, plot4; layout = 4, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4" ], plot_title = "ActiveParticle1 (R=0.5µm, v=0µm/s)", plot_titlevspan=0.1)




## RIPETIAMO PER UNA PARTICELLA DI RAGGIO 0.5MICRON E V=5µm/s ##

#inizializzazione particella (R=0.5e-6m , v=5µm/s)
abp1 = initABP( (0.0, 0.0, orientazione), 0.5, 5.0 );

#creo il primo grafico grafico 
N = 100000
p, t = trajectory( abp1, N);

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
p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],aspect_ratio= :equal,legend=false)
#scatter!([x[end]],[y[end]], legend=false)


p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],aspect_ratio= :equal,legend=false)
#scatter!([x[end]],[y[end]], legend=false)


p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,range=[-100,75],aspect_ratio= :equal,legend=false)
#scatter!([x[end]],[y[end]], legend=false)

    
#subplot e metto i 4 grafici in un unico grafico
#plot(plot1, plot2 , plot3, plot4; layout = 4, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4" ], plot_title = "ActiveParticle2 (R=0.5µm, v=5µm/s)", plot_titlevspan=0.1)

# subplot deei due grafici della stessa particella:
plot(plot1, plot2; layout = 2, plot_title = "ActiveParticle1 (R = 0.5µm, θ = $orientazione rad)", title = ["v = 0µm/s" "v = 5µm/s"], plot_titlevspan = 0.1)



## RIPETIAMO PER UNA PARTICELLA DI RAGGIO 2.5µm E v=0µm/s ##

#inizializzazione particella (R=2.5e-6m , v=0µm/s)
orientazione2 = rand(0.0:pi/2)
abp1 = initABP( (0.0, 0.0, orientazione2), 2.5, 0.0 );

#creo il primo grafico grafico 
N = 100000
p, t = trajectory( abp1, N);

#NOTA: p è un vettore di tuple, ovvero è un vettore con N elementi 
# e ciascun elemento è una tupla (x,y, teta)

#creo i vettori degli assi 
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]

#NOTA: questo comando mi permette di creare un vettore i cui elementi 
# sono dati dai primi elementi di pi dove pi 
# assume tutti i valori dentro p (cioè parte da p)

#plotto in x y
plot3 = plot(x,y, xlim=[-25,200],ylim=[-100,100], title = "ActiveParticle3 (R=2.5µm, v=0µm/s)", aspect_ratio= :equal) 
#scatter!([x[end]],[y[end]], legend=false)
xlabel!("x [μm]")
ylabel!("y [μm]")

#realizzo più grafici nello stesso plot (subplot) 
# (vedere se riesco a migliorare questa parte)
p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,xlim=[-25,200],ylim=[-100,100],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,xlim=[-25,200],ylim=[-100,100],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,xlim=[-25,200],ylim=[-100,100],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)
    
#subplot e metto i 4 grafici in un unico grafico
#plot(plot1, plot2 , plot3, plot4; layout = 4, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4"], plot_title = "ActiveParticle3 (R=2.5µm, v=0µm/s)", plot_titlevspan=0.1)






## RIPETIAMO PER UNA PARTICELLA DI RAGGIO 2.5MICRON E v=5µm/s ##

#inizializzazione particella (R=2.5e-6m , v=5µm/s)
abp1 = initABP( (0.0, 0.0, orientazione2), 2.5, 5.0 );

#creo il primo grafico grafico 
N = 100000
p, t = trajectory( abp1, N);

#NOTA: p è un vettore di tuple, ovvero è un vettore con N elementi 
# e ciascun elemento è una tupla (x,y, teta)

#creo i vettori degli assi 
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]

#NOTA: questo comando mi permette di creare un vettore i cui elementi 
# sono dati dai primi elementi di pi dove pi 
# assume tutti i valori dentro p (cioè parte da p)

#plotto in x y
plot4 = plot(x,y,xlim=[-25,200],ylim=[-100,100], title = "ActiveParticle4 (R=2.5µm, v=5µm/s)", aspect_ratio= :equal) 
#scatter!([x[end]],[y[end]], legend=false)
xlabel!("x [μm]")
ylabel!("y [μm]")

#realizzo più grafici nello stesso plot (subplot) 
# (vedere se riesco a migliorare questa parte)
p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,xlim=[-25,200],ylim=[-100,100],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,xlim=[-25,200],ylim=[-100,100],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot!(x,y,xlim=[-25,200],ylim=[-100,100],legend=false,aspect_ratio= :equal)
#scatter!([x[end]],[y[end]], legend=false)
    
#subplot e metto i 4 grafici in un unico grafico
#plot(plot1, plot2 , plot3, plot4; layout = 4, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4"], plot_title = "ActiveParticle4 (R=2.5µm, v=5µm/s)", plot_titlevspan=0.1)

#subplot per visualizzare plot3 e plot4 vicini
plot(plot3, plot4; layout = 2, plot_title = "ActiveParticle2 (R = 2.5µm, θ = $orientazione2 rad)", title = ["v = 0µm/s" "v = 5µm/s"], plot_titlevspan = 0.1)




