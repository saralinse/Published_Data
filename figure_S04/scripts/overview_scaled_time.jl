using CairoMakie, Printf, Roots

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

function diffusion_limited_gamma(
	c::Float64, # [1/um^3]
	t::Float64; # [s]
	diffCoeff::Float64 = 45.0, # [um^2/s] (Rh = 4.9 nm)
	R::Float64 = 2100.0, # [um]
	N::Int = 100) # []
	
	Minf = ( c*(4*pi*R^3 / 3) ) / ( 4*pi*R^2 ) # [1/um^2]
	M = Minf*( 1 - 6/pi^2 * sum([exp(-diffCoeff*n^2*pi^2*t / R^2)/n^2 for n in 1:N]) ) # [1/um^2]
	
	return M::Float64
	
end

function surface_tension(
	gamma::Float64; # [1/um^2]
	y0::Float64 = 72.0E-3, # [J/m^2]
	gamma0::Float64 = 24.0E3) # [1/um^2]
	
	p1 = 1.029
	p2 = -32.4
	p3 = 77.4
	p4 = 49.8
	
	f(gamma) = y0 - 1E-3*log(p1^(p2*gamma/gamma0 + p3) + p1^p4)/log(p1)
	
	pressure = max(0.0, f(gamma)) # [J/m^2]
	
	surfaceTension = y0 - pressure
	
	return surfaceTension::Float64 # [J/m^2]
	
end

#####################################

include("data.jl")

plotTheme = Theme(
	size = (600, 400),
	figure_padding = (10, 150, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 2,
	markersize = 15,
	Axis = (
		xlabel = rich("t/t", subscript("ind")),
		ylabel = rich("Surface tension / mJ⋅m", superscript("-2")),
		limits = ((0, 16), (52.5, 75)),
		xticks = 0:2:16,
		# limits = ((1E-1, 1E3), (52.5, 75)),
		# xscale = log10,
		# xticks = LogTicks(-3:3),
		# xminorticks = IntervalsBetween(9),
		yticks = 50:5:75,
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

colors = cgrad([RGBf(1/4, 0, 1/3), RGBf(1/4, 1/4, 2/3), RGBf(0, 3/4, 3/4), RGBf(0, 3/4, 0), RGBf(1, 3/4, 1/8)],
	length(uniqueConcentrations);
	alpha = 0.7,
	categorical = true,
	rev = false)

for i in reverse(eachindex(fileNames))
	
	(concentration, timePoint, surfaceTension) = read_data_from_file("data/"*fileNames[i])
	colorIndex = findfirst(uniqueConcentrations .== concentration)
	
	#if concentration == uniqueConcentrations[1]
	if round(concentration; sigdigits = 2) == 0.99 # HACK
		concString = (@sprintf "%#.2g" concentration)*" µM"
	else
		concString = @sprintf "%#.2g" concentration
	end
	
	# gamma0 = 24E3 # [1/um^2]
	# tind = diffusion_limited_induction_time(concentration*1E-6*6.022E23*1E-15, gamma0) # [s]
	# lines!(axis, timePoint[10 .<= timePoint .<= 5E3]/tind, surfaceTension[10 .<= timePoint .<= 5E3];
		# color = colors[colorIndex],
		# label = concString)
	
	k = findfirst([tinds[j][4] for j in eachindex(tinds)] .== fileNames[i])
	if !isnothing(k)
		
		lines!(axis, timePoint[timePoint .<= 10E3]/tinds[k][2], surfaceTension[timePoint .<= 10E3];
			color = colors[colorIndex],
			label = concString)
		
	end
	
	vlines!(axis, [1];
		linestyle = :dot,
		color = RGBAf(0, 0, 0, 1),
		label = nothing)
	
end

axislegend(axis;
	titlefont = :regular,
	titlesize = 15,
	labeljustification = :left,
	position = (1.55, 0.5),
	backgroundcolor = :transparent,
	framevisible = false,
	unique = true,
	labelsize = 20)

save("output/overview_scaled_time.pdf", figure)

#####################################

include("data_model.jl")

plotTheme = Theme(
	size = (600, 400),
	figure_padding = (10, 150, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 2,
	markersize = 15,
	Axis = (
		xlabel = rich("t/t", subscript("ind")),
		ylabel = rich("Surface tension / mJ⋅m", superscript("-2")),
		limits = ((0, 16), (52.5, 75)),
		xticks = 0:2:16,
		# limits = ((1E-1, 1E3), (52.5, 75)),
		# xscale = log10,
		# xticks = LogTicks(-3:3),
		# xminorticks = IntervalsBetween(9),
		yticks = 50:5:75,
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

colors = cgrad([RGBf(1/4, 0, 1/3), RGBf(1/4, 1/4, 2/3), RGBf(0, 3/4, 3/4), RGBf(0, 3/4, 0), RGBf(1, 3/4, 1/8)],
	length(uniqueConcentrations);
	alpha = 0.7,
	categorical = true,
	rev = false)

for i in reverse(eachindex(fileNames))
	
	(concentration, timePoint, surfaceTension) = read_data_from_file("data_model/"*fileNames[i])
	colorIndex = findfirst(uniqueConcentrations .== concentration)
	
	#if concentration == uniqueConcentrations[1]
	if round(concentration; sigdigits = 2) == 0.99 # HACK
		concString = (@sprintf "%#.2g" concentration)*" µM"
	else
		concString = @sprintf "%#.2g" concentration
	end
	
	# gamma0 = 24E3 # [1/um^2]
	# tind = diffusion_limited_induction_time(concentration*1E-6*6.022E23*1E-15, gamma0) # [s]
	# lines!(axis, timePoint[10 .<= timePoint .<= 5E3]/tind, surfaceTension[10 .<= timePoint .<= 5E3];
		# color = colors[colorIndex],
		# label = concString)
	
	(_, timePointDifflim, surfaceTensionDifflim) = read_data_from_file("data_model_difflim/"*fileNames[i])
	
	k = findfirst([tinds[j][4] for j in eachindex(tinds)] .== fileNames[i])
	#if !isnothing(k)
	if concentration < 1.0 # Model
		
		# lines!(axis, timePointDifflim[timePointDifflim .<= 10E3]/tinds[k][2], surfaceTensionDifflim[timePointDifflim .<= 10E3];
			# linestyle = :dash,
			# color = colors[colorIndex],
			# label = nothing)
		
		gamma0 = 24E3 # [1/um^2]
		ts = 1.0:1.0:10000.0
		lines!(axis, ts/diffusion_limited_induction_time(concentration*1E-6*6.022E23*1E-15, gamma0), 1E3*surface_tension.(diffusion_limited_gamma.(concentration*1E-6*6.022E23*1E-15, ts));
			linestyle = :dash,
			color = colors[colorIndex],
			label = nothing)
		
		lines!(axis, timePoint[timePoint .<= 10E3]/tinds[k][2], surfaceTension[timePoint .<= 10E3];
			linestyle = :solid,
			color = colors[colorIndex],
			label = concString)
		
	end
	
	vlines!(axis, [1];
		linestyle = :dot,
		color = RGBAf(0, 0, 0, 1),
		label = nothing)
	
end

axislegend(axis;
	titlefont = :regular,
	titlesize = 15,
	labeljustification = :left,
	position = (1.55, 0.5),
	backgroundcolor = :transparent,
	framevisible = false,
	unique = true,
	labelsize = 20)

save("output/overview_scaled_time_model.pdf", figure)
