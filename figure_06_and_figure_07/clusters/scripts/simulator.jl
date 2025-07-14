using CairoMakie

include("adsorption_model.jl")

###############################

cbulks = [ # [1/um^3]
round(5.000*1E-6*6.022E23*1E-15; sigdigits = 4);
round(3.330*1E-6*6.022E23*1E-15; sigdigits = 4);
round(2.220*1E-6*6.022E23*1E-15; sigdigits = 4);
round(1.480*1E-6*6.022E23*1E-15; sigdigits = 4);
round(0.988*1E-6*6.022E23*1E-15; sigdigits = 4);
round(0.658*1E-6*6.022E23*1E-15; sigdigits = 4);
round(0.439*1E-6*6.022E23*1E-15; sigdigits = 4);
round(0.293*1E-6*6.022E23*1E-15; sigdigits = 4);
round(0.195*1E-6*6.022E23*1E-15; sigdigits = 4);
round(0.130*1E-6*6.022E23*1E-15; sigdigits = 4);
round(0.087*1E-6*6.022E23*1E-15; sigdigits = 4)]

# cbulks = [round(0.988*1E-6*6.022E23*1E-15; sigdigits = 4)] # [1/um^3]

tmax = 5000.0 # [s]
xmax = 2100.0 # [um] (drop)

dataset = Vector{Tuple{Float64, Vector{Tuple{Float64, Float64}}, Vector{Tuple{Float64, Float64}}}}()

for cbulk in cbulks
	println("Concentration = "*string(cbulk)*" [1/um^3]")
	gammas, surfaceTensions = simulate(cbulk, tmax, xmax)
	push!(dataset, (cbulk, gammas, surfaceTensions))
end

################################

plotTheme = Theme(
	size = (600, 400),
	figure_padding = (10, 150, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 3,
	markersize = 15,
	Axis = (
		xlabel = "Time, t [s]",
		ylabel = "Surface tension, γ [mJ/m²]",
		limits = ((0, tmax), (52.5, 75)),
		xticks = 0:1000:tmax,
		yticks = 0:5:100,
		xticksmirrored = true,
		yticksmirrored = true,
		xminorticksvisible = true,
		yminorticksvisible = true,
		spinewidth = 2,
		xgridvisible = false,
		ygridvisible = false,
		xminorgridvisible = false,
		yminorgridvisible = false)
	)

set_theme!(plotTheme)

figure = Figure()
axis = Axis(figure[1, 1])

colors = cgrad(:viridis, length(cbulks); alpha = 0.5, categorical = true)

for data in dataset
	
	cbulk, _, surfaceTensions = data
	
	concentrationString = lpad(split(string(round((cbulk/6.022E23)*1E15*1E6; sigdigits = 3)), ".")[1], 2, "0")*"."*rpad(split(string(round((cbulk/6.022E23)*1E15*1E6; sigdigits = 3)), ".")[2], 3, "0")
	io = open("output/data/model_"*replace(concentrationString, "." => "_")*"_uM.csv", "w")
	println(io, "IDs:;model_"*replace(concentrationString, "." => "_")*"_uM")
	println(io, "Concentration [uM]:;"*concentrationString)
	println(io, "Time [s];Surface tension [mJ/m^2];Area [mm^2];Volume [mm^3]")
	for surfaceTension in surfaceTensions
		println(io, string(surfaceTension[1])*";"*string(1E3*surfaceTension[2])*";Inf;Inf")
	end
	close(io)
	
	colorIndex = findfirst(cbulks .== cbulk)
	lines!(axis, [(surfaceTension[1], 1E3*surfaceTension[2]) for surfaceTension in surfaceTensions];
		color = colors[colorIndex],
		label = string(round((cbulk/6.022E23)*1E15*1E6; sigdigits = 3))* " μM")
	
end

axislegend(axis;
	labeljustification = :left,
	position = (1.6, 0.5),
	backgroundcolor = :transparent,
	framevisible = false,
	unique = true,
	labelsize = 20)

save("output/model_surface_tension.pdf", figure)
