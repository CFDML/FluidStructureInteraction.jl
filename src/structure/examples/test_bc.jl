using Peridynamics
include("thermstep2.0.jl")
include("therm2.0.jl")
include("coth.jl")

l = 0.02
w = 0.02
h = 0.002
Δx = l/100
δ = 3.015Δx 
dt = 1e-5
nth = 45000
b = 1.0

pos, vol = uniform_box(l, w+6Δx, h, Δx; center_x=-0Δx, center_y=0, center_z=0)

body = Body(BBTMaterial(), pos, vol)
material!(body; horizon=3.015Δx, E=2.0e11, rho=1e6, epsilon_c=b, 
          kc=114.0, aph=1.0e-6, cv=1.0, rft=273, dx = Δx, h=1e5, hσ=1.0, hϵ=1.0, tem∞=0.0) 
point_set!(p -> (p[2] < -w/2), body, :tem_area1)
point_set!(p -> (p[2] >  w/2), body, :tem_area2)

temperature_bc!(t -> -100.0, body, :tem_area1)
temperature_bc!(t -> 100.0, body, :tem_area2)

point_set!(p -> (abs(p[1]) > (l/2-Δx)), body, :conv_area)

for i in 1:2
    let 
        t = 2π*i/6
        spx = 0.0033 * cos(t)
        spy = -sqrt(3)*0.0033 + 0.0033 * sin(t)
        sp = [spx, spy]
        lc = 0.006
        ex = cos(t) 
        ey = sin(t)  
        nx = -ey
        ny = ex

        u_n = Symbol("cu_$(i)_area")
        d_n = Symbol("cd_$(i)_area")

        point_set!(p -> (0 < (p[1]-sp[1])*ex + (p[2]-sp[2])*ey < lc) && 
        (0.0 < (p[1]-sp[1])*nx + (p[2]-sp[2])*ny < 3Δx) , body, u_n)
    
        point_set!(p -> (0 < (p[1]-sp[1])*ex + (p[2]-sp[2])*ey < lc) && 
        (0.0 > (p[1]-sp[1])*nx + (p[2]-sp[2])*ny > -3Δx) , body, d_n)
    
        precrack!(body, u_n, d_n)       
    end
end

vv = Thermstep(steps=nth, stepsize = dt)
job = Job(body, vv; freq=100, path="results/therm_crack_conv", fields=(:temperature))
A = cou_submit(job, "therm")
