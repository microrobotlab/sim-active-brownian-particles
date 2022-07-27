#Script Lapo Corti per simulazione movimento particella attiva
#OBIETTIVO: realizzare 4 grafici diversi (ognuno per una particella 
# di dimensione nota e velocità iniziale nota) e all'interno di ciascun 
# grafico vado a plottare 4-5 traiettorie diverse della stessa particella. 


include("ActiveBrownianParticles2.jl");

using BenchmarkTools
using Plots 
using Distances 

## PARTIAMO DAL PRIMO GRAFICO PER UNA PARTICELLA DI RAGGIO 1.5MICRON E V=0 ##

#inizializzazione particella (R=1.5e-6m , v=0m/s)
abp1 = initABP( (0.0, 0.0, 0.0), 1.5, 0.0 );

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
plot1 = plot_trajectory(x,y)

#realizzo più grafici nello stesso plot (subplot) 
# (vedere se riesco a migliorare questa parte)
p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot2 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot3 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot4 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot5 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot6 = plot_trajectory(x,y)
    
#subplot e metto i 6 grafici in un unico grafico
plot(plot1, plot2 , plot3, plot4, plot5, plot6 ; layout = 6, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4"  "Trajectory 5" "Trajectory 6"], plot_title = "ActiveParticle1 (R=1.5µm, v=0µm/s)", plot_titlevspan=0.1)

 



## RIPETIAMO PER UNA PARTICELLA DI RAGGIO 1.5MICRON E V=10µm/s ##

#inizializzazione particella (R=1.5e-6m , v=10µm/s)
abp1 = initABP( (0.0, 0.0, 0.0), 1.5, 10.0 );

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
plot1 = plot_trajectory(x,y)

#realizzo più grafici nello stesso plot (subplot) 
# (vedere se riesco a migliorare questa parte)
p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot2 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot3 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot4 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot5 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot6 = plot_trajectory(x,y)
    
#subplot e metto i 6 grafici in un unico grafico
plot(plot1, plot2 , plot3, plot4, plot5, plot6 ; layout = 6, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4"  "Trajectory 5" "Trajectory 6"], plot_title = "ActiveParticle2 (R=1.5µm, v=10µm/s)", plot_titlevspan=0.1)





## RIPETIAMO PER UNA PARTICELLA DI RAGGIO 3µm E V=0µm/s ##

#inizializzazione particella (R=3e-6m , v=0µm/s)
abp1 = initABP( (0.0, 0.0, 0.0), 3.0, 0.0 );

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
plot1 = plot_trajectory(x,y)

#realizzo più grafici nello stesso plot (subplot) 
# (vedere se riesco a migliorare questa parte)
p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot2 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot3 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot4 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot5 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot6 = plot_trajectory(x,y)
    
#subplot e metto i 6 grafici in un unico grafico
plot(plot1, plot2 , plot3, plot4, plot5, plot6 ; layout = 6, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4"  "Trajectory 5" "Trajectory 6"], plot_title = "ActiveParticle3 (R=3µm, v=0µm/s)", plot_titlevspan=0.1)






## RIPETIAMO PER UNA PARTICELLA DI RAGGIO 1.5MICRON E V=10µm/s ##

#inizializzazione particella (R=1.5e-6m , v=10µm/s)
abp1 = initABP( (0.0, 0.0, 0.0), 1.5, 10.0 );

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
plot1 = plot_trajectory(x,y)

#realizzo più grafici nello stesso plot (subplot) 
# (vedere se riesco a migliorare questa parte)
p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot2 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot3 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot4 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot5 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot6 = plot_trajectory(x,y)
    
#subplot e metto i 6 grafici in un unico grafico
plot(plot1, plot2 , plot3, plot4, plot5, plot6 ; layout = 6, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4"  "Trajectory 5" "Trajectory 6"], plot_title = "ActiveParticle2 (R=1.5µm, v=10µm/s)", plot_titlevspan=0.1)





## RIPETIAMO PER UNA PARTICELLA DI RAGGIO 3µm E V=10µm/s ##

#inizializzazione particella (R=3e-6m , v=10µm/s)
abp1 = initABP( (0.0, 0.0, 0.0), 3.0, 10.0 );

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
plot1 = plot_trajectory(x,y)

#realizzo più grafici nello stesso plot (subplot) 
# (vedere se riesco a migliorare questa parte)
p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot2 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot3 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot4 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot5 = plot_trajectory(x,y)

p, t = trajectory( abp1, N);
x = [pi[1] for pi in p]
y = [pi[2] for pi in p]
plot6 = plot_trajectory(x,y)
    
#subplot e metto i 6 grafici in un unico grafico
plot(plot1, plot2 , plot3, plot4, plot5, plot6 ; layout = 6, title = ["Trajectory 1" "Trajectory 2" "Trajectory 3" "Trajectory 4"  "Trajectory 5" "Trajectory 6"], plot_title = "ActiveParticle4 (R=3µm, v=10µm/s)", plot_titlevspan=0.1)



