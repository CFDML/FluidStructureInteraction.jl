using Peridynamics
include("thermstep2.0.jl")
include("therm2.0.jl")
include("coth.jl")

l = 0.04
w = 0.01
Δx = l/50
δ = 3.015Δx 
nth = 500
b = 1.0
p_const = 3e6/Δx

pos, vol = uniform_box(l, w, Δx, Δx; center_x=-0Δx, center_y=0, center_z=0)

body = Body(BBTMaterial(), pos, vol)
material!(body; horizon=3.015Δx, E=72.4e9, rho=1850.0, epsilon_c=b, 
          kc=35.0, aph=8.59e-5, cv=1106.0, rft=273.0, dx = Δx, h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, D=2) 
          
       
#flux_index = findall(x -> x <= Δx - l/2, pos[1, :])
hsource_matrix = zeros(1, size(pos, 2))
for i in  1:size(pos, 2)
    hsource_matrix[i] = p_const
end

#hsource_bc!(hsource_matrix, body, :all_points)
point_set!(p -> (p[1] < Δx-l/2), body, :flux_area)
hsource_bc!(hsource_matrix, body, :flux_area)

vv = Thermstep(steps=nth)
job = Job(body, vv; freq=10, path="results/test", fields=(:temperature, :hsource))
A = cou_submit(job, "therm")
