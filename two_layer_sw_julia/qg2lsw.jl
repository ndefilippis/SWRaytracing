function qg2layersw(nx, T_Fr_days, U_g, f, Cg)
#= 
Input:
nx: resolution of QG flow
T_Fr_days: number of days to simulate
U_g: amplitude of background geostrophic flow
f: Coriolis parameter
Cg: Group velocity of the waves (equal to sqrt(gH))
=#

L = 20;
dx = L/nx;
x = 
end