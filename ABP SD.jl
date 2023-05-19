
# PURPOSE: Calculate the mean and standard deviation of the packing fraction
# INPUT: Data file genereated by ABP analysis.jl

using BenchmarkTools, Plots, LaTeXStrings, Statistics, CSV, DataFrames

gr()

function stat_analysis1(a,b,R,pathf)
  f1= pathf*"_p.csv"
  df= CSV.read(f1, DataFrame)
  neq= df[!,:p1]
  npole= df[!,:p2]
  eθ = atan(b/a)   # angle at which area will be same
pf_factor = R*R
Aeq= a*b*(atan(a*tan(eθ)/b))  # equator area

Ap= a*b*(atan(b/(tan(eθ)*a)))
  meq= mean(neq)
  sdeq= stdm(neq,meq)    # standard deviation of equator
  mpole= mean(npole)
  sdpole= stdm(npole,mpole)
  mpfe = meq*(0.5*π/Aeq)*pf_factor 
  mpfp = mpole*(0.5*π/Ap)*pf_factor 
  println("i am in ABP SD")
return meq, sdeq, mpole,sdpole
#return mpfe, mpfp
end