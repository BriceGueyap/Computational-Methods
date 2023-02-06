#=-------------------------------------------------------------------------------
                    Problem Set 3

            Brice Romuald Gueyap Kounga
                    February 2023
      Note: This is a joint work with Fengfan Xiang
------------------------------------------------------------------------------=#

cd() # Go to root directory
cd("/Users/bgueyap/Library/CloudStorage/GoogleDrive-bricegueyap@gmail.com/My Drive/University - Studies/University of Western Ontario/Winter 2023/Computational-Methods/Problem Sets/PS3")
mkpath("Figures")
using Plots

# CRRA utility function with σ = 2
σ = 2
u(c) = c.^(1-σ)./(1-σ)

# Log utility function
u1(c) = log.(c)

# square root utility function
u2(c) = sqrt.(c)

# Fine Grid
caxis = range(0.05,2;length=10000) ;

# Exact evaluation of utility functions
crra = map(u,caxis)
logu = map(u1,caxis)
sqru = map(u2,caxis)

#-------------------------------------------------------------------------------
# Interpolation Routine 1: Monomial

function monomial_interpolation(x,x_grid::Array,y_grid::Array)
    N = length(x_grid) # Number of nodes

    A = ones(N,1)
    for i = 1:N-1
        A = hcat(A,x_grid.^i)
    end

    a = A\y_grid # Monomial interpolation coefficients

    V_hat = 0
    for i = 1:N
        V_hat = V_hat + a[i]*x^(i-1)
    end

    return V_hat
end

#-------------------------------------------------------------------------------
# Interpolation Routine 2: Spline Piecewise Linear

function spline_linear_interpolation(x,x_grid::Array,y_grid::Array)
    N = length(x_grid)
    ind = findmax(sign.(x_grid .- x))[2] - 1
    A_x = (x_grid[ind+1] - x)/(x_grid[ind+1]-x_grid[ind])

    V_hat = A_x*y_grid[ind] + (1-A_x)*y_grid[ind+1]

    return V_hat
end

#-------------------------------------------------------------------------------
# Interpolation Routine 3: Spline Cubic

function spline_cubic_interpolation(x,x_grid::Array,y_grid::Array)
    N = length(x_grid)
    h = zeros(N-1)
    q = zeros(N-1)
    for i = 1:N-1
        h[i] = x_grid[i+1] - x_grid[i]
        q[i] = y_grid[i+1] - y_grid[i]
    end

   M1 = hcat(1,zeros(1,N-1))
   M2 = vcat([h[1], 2*(h[1]+h[2]),h[2]], zeros(N-3,1))
   M = vcat(M1,M2')

   for i = 2:N-3
       temp = vcat(zeros(i-1,1),[h[i], 2*(h[i]+h[i+1]),h[i+1]])
       temp = vcat(temp,zeros(N-i-2,1))
       M = vcat(M,temp')
   end

   Mm = vcat(zeros(N-3,1),[h[N-2], 2*(h[N-2]+h[N-1]),h[N-1]])
   Mn = hcat(zeros(1,N-1),1)
   Mn = vcat(Mm',Mn)
   M = vcat(M,Mn)./6
   b = zeros(N)

   for i = 1:N-2
       b[i+1] = q[i+1]/h[i+1]-q[i]/h[i]
   end

   a = M\b

   ind = findmax(sign.(x_grid .- x))[2] - 1

   x1 = (x_grid[ind+1]-x)/(x_grid[ind+1]-x_grid[ind])
   x2 = (x-x_grid[ind])/(x_grid[ind+1]-x_grid[ind])

   V_hat  = x1*y_grid[ind] + x2*y_grid[ind+1] +
   (1/6)*(x_grid[ind+1]-x_grid[ind])^2*
   (a[ind]*(x1^3-x1)+a[ind+1]*(x2^3-x2))

   return V_hat
end

#-------------------------------------------------------------------------------
## CRRA
# Plot Monomial Interpolation
gr()
plt = plot(title="Monimial Interpolation")
plot!(caxis,crra,linewidth=3,label = "CRRA Utility Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = collect(range(0.05,2;length=n)); # Collect makes it an array instead of a collection
    yi = map(u,xi);                       # the corresponding y-coordinates

    # Monomial Interpolation
    interp = map(x->monomial_interpolation(x,xi,yi),caxis)
    plot!(caxis,interp,linewidth=3,label="n=$n")
end
savefig("./Figures/CRRA_Monimial_Polynomial.png")

# Plot Spline Piecewise Linear Interpolation
gr()
plt = plot(title="Spline Linear Interpolation")
plot!(caxis,crra,linewidth=3,label = "CRRA Utility Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = collect(range(0.05,2;length=n)); # Collect makes it an array instead of a collection
    yi = map(u,xi);                       # the corresponding y-coordinates

    # Spline Piecewise Linear Interpolation
    interp = map(x->spline_linear_interpolation(x,xi,yi),caxis)
    plot!(caxis,interp,linewidth=3,label="n=$n")
end
savefig("./Figures/CRRA_Spline_Linear.png")

# Plot Spline Cubic Interpolation
gr()
plt = plot(title="Spline Cubic Interpolation")
plot!(caxis,crra,linewidth=3,label = "CRRA Utility Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = collect(range(0.05,2;length=n)); # Collect makes it an array instead of a collection
    yi = map(u,xi);                       # the corresponding y-coordinates

    # Spline Cubic Interpolation
    interp = map(x->spline_cubic_interpolation(x,xi,yi),caxis)
    plot!(caxis,interp,linewidth=3,label="n=$n")
end
savefig("./Figures/CRRA_Spline_Cubic.png")

## Log
# Plot Monomial Interpolation
gr()
plt = plot(title="Monimial Interpolation")
plot!(caxis,logu,linewidth=3,label = "Log Utility Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = collect(range(0.05,2;length=n)); # Collect makes it an array instead of a collection
    yi = map(u1,xi);                       # the corresponding y-coordinates

    # Monomial Interpolation
    interp = map(x->monomial_interpolation(x,xi,yi),caxis)
    plot!(caxis,interp,linewidth=3,label="n=$n")
end
savefig("./Figures/Log_Monimial_Polynomial.png")

# Plot Spline Piecewise Linear Interpolation
gr()
plt = plot(title="Spline Linear Interpolation")
plot!(caxis,logu,linewidth=3,label = "Log Utility Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = collect(range(0.05,2;length=n)); # Collect makes it an array instead of a collection
    yi = map(u1,xi);                       # the corresponding y-coordinates

    # Spline Piecewise Linear Interpolation
    interp = map(x->spline_linear_interpolation(x,xi,yi),caxis)
    plot!(caxis,interp,linewidth=3,label="n=$n")
end
savefig("./Figures/Log_Spline_Linear.png")

# Plot Spline Cubic Interpolation
gr()
plt = plot(title="Spline Cubic Interpolation")
plot!(caxis,logu,linewidth=3,label = "Log Utility Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = collect(range(0.05,2;length=n)); # Collect makes it an array instead of a collection
    yi = map(u1,xi);                       # the corresponding y-coordinates

    # Spline Cubic Interpolation
    interp = map(x->spline_cubic_interpolation(x,xi,yi),caxis)
    plot!(caxis,interp,linewidth=3,label="n=$n")
end
savefig("./Figures/Log_Spline_Cubic.png")

## square root
# Plot Monomial Interpolation
gr()
plt = plot(title="Monimial Interpolation")
plot!(caxis,sqru,linewidth=3,label = "Sqrt Utility Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = collect(range(0.05,2;length=n)); # Collect makes it an array instead of a collection
    yi = map(u2,xi);                       # the corresponding y-coordinates

    # Monomial Interpolation
    interp = map(x->monomial_interpolation(x,xi,yi),caxis)
    plot!(caxis,interp,linewidth=3,label="n=$n")
end
savefig("./Figures/Sqrt_Monimial_Polynomial.png")

# Plot Spline Piecewise Linear Interpolation
gr()
plt = plot(title="Spline Linear Interpolation")
plot!(caxis,sqru,linewidth=3,label = "Sqrt Utility Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = collect(range(0.05,2;length=n)); # Collect makes it an array instead of a collection
    yi = map(u2,xi);                       # the corresponding y-coordinates

    # Spline Piecewise Linear Interpolation
    interp = map(x->spline_linear_interpolation(x,xi,yi),caxis)
    plot!(caxis,interp,linewidth=3,label="n=$n")
end
savefig("./Figures/Sqrt_Spline_Linear.png")

# Plot Spline Cubic Interpolation
gr()
plt = plot(title="Spline Cubic Interpolation")
plot!(caxis,sqru,linewidth=3,label = "Sqrt Utility Function",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
for n = (4,6,11,21)
    # Grid of nodes for interpolation
    xi = collect(range(0.05,2;length=n)); # Collect makes it an array instead of a collection
    yi = map(u2,xi);                       # the corresponding y-coordinates

    # Spline Cubic Interpolation
    interp = map(x->spline_cubic_interpolation(x,xi,yi),caxis)
    plot!(caxis,interp,linewidth=3,label="n=$n")
end
savefig("./Figures/Sqrt_Spline_Cubic.png")
