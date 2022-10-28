#this code will check and calculate the gradients 
# JS
using ForwardDiff
#--------------------------------------------------------------------------------------------

# here one element is called 
#=
function  scale(c::Float64)
    x = zeros(length(2),2)

    x = c * rand(2,2) # it makes a 2x1 matric
    #x = c * ones(2,2) # it makes a 2x1 matric
    check = grad(x)
    
     return check
end

function grad(x::Array{Float64})
     
    f(x) = 4x[1]^2 + 3x[1] + 2x[2]
    
   
    df=ForwardDiff.gradient(f, [x[1],x[2]])
    
   
    return df
end


result1=scale(10.0) 

=#
#--------------------------------------------------------------------------------------------

# I am calling the functions and changing the x values randomly like in our code
# here mutiple element are called using the dot function

function  scale(c::Float64)
    
    x= [c*rand(2) for _ in 1:5] # random elements vector
    check = grad.(x) # check is a vector of the gradient w.r.t x--> x[1] and y --> x[2] resp
    
     return check
end

function grad(x::Array{Float64})
     
    f(x) = 4x[1]^2 + 3x[1] + 10*(sin.(x[2]))
    df=ForwardDiff.gradient(f, [x[1],x[2]]) # calculate gradient of function f at points x[1] amd x[2]
    
    return df
end

result=scale(10.0) # a scaling factor input

