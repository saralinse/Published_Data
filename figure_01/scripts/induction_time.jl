using CairoMakie, Roots

function diffusion_limited_induction_time(
	c::Float64, # [1/um^3]
	gamma0::Float64; # [1/um^2]
	diffCoeff::Float64 = 45.0, # [um^2/s] (Rh = 4.9 nm)
	R::Float64 = 2100.0, # [um]
	N::Int = 100) # []
	
	function diffusion_limited_gamma(
		c::Float64, # [1/um^3]
		t::Float64) # [s]
		
		M = Minf*( 1 - 6/pi^2 * sum([exp(-diffCoeff*n^2*pi^2*t / R^2)/n^2 for n in 1:N]) ) # [1/um^2]
		
		return M::Float64
		
	end
	
	Minf = ( c*(4*pi*R^3 / 3) ) / ( 4*pi*R^2 ) # [1/um^2]
	
	if gamma0 >= Minf
		tind = Inf # [s]
	else
		tind = max(0.0, find_zero(t -> diffusion_limited_gamma(c, t) - gamma0, 0)) # [s]
	end
	
	return tind::Float64
	
end

#####################################

include("data.jl")

plotTheme = Theme(
	size = (400, 400),
	figure_padding = (30, 20, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 2,
	markersize = 15,
	Axis = (
		xlabel = rich("[JB6] / μM"),
		ylabel = rich("Induction time / s"),
		limits = ((0.05, 10), (1, 100000)),
		xscale = log10,
		xticks = LogTicks(-2:2),
		xminorticks = IntervalsBetween(9),
		yscale = log10,
		yticks = LogTicks(-1:5),
		yminorticks = IntervalsBetween(9),
		xticksmirrored = true,
		yticksmirrored = true,
		xminorticksvisible = true,
		yminorticksvisible = true,
		spinewidth = 2,
		xgridvisible = false,
		ygridvisible = false)
	)

set_theme!(plotTheme)

figure = Figure()
axis = Axis(figure[1, 1])

color = RGBAf(0.0, 0.0, 0.0, 0.75)
scatter!(axis, [(tinds[i][1], tinds[i][2]) for i in eachindex(tinds)],
	color = color,
	marker = :cross,
	label = nothing)
cs = [10 .^(-1.5:0.0001:-1); 10 .^(-1:0.01:2)] # [uM]
gamma0 = 24E3 # [1/um^2]
lines!(axis, cs, diffusion_limited_induction_time.(cs*1E-6*6.022E23*1E-15, gamma0);
	linestyle = :dot,
	color = color,
	label = rich("Γ", subscript("0"), "= ", string(round(Int, gamma0/1E3)), "·10", superscript("3"), " μm", superscript("-2")))
weights = [1; 1/2; 1/2; 1; 1; 1; 1/4; 1/4; 1/4; 1/4; 1/4; 1/4; 1/4; 1/4]
deviation = sum([weights[i]*( log(diffusion_limited_induction_time.(tinds[i][1]*1E-6*6.022E23*1E-15, gamma0)) - log(tinds[i][2]) )^2 for i in eachindex(tinds)])
println("Γ0 = "*string(gamma0)*" [1/μm^2]. Diffusion model deviation: "*string(deviation))

axislegend(axis;
	titlefont = :regular,
	titlesize = 20,
	labeljustification = :left,
	position = :rt,
	backgroundcolor = :transparent,
	framevisible = false,
	unique = true,
	labelsize = 20)

save("output/induction_time.pdf", figure)

#####################################
