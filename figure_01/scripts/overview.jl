using CairoMakie, Printf

include("data.jl")

plotTheme = Theme(
	size = (600, 400),
	figure_padding = (10, 150, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 2,
	markersize = 15,
	Axis = (
		xlabel = rich("Time / s"),
		ylabel = rich("Surface tension / mJ⋅m", superscript("-2")),
		limits = ((0, 5000), (52.5, 75)),
		xticks = 0:1000:5000,
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
	
	if concentration == uniqueConcentrations[1]
		concString = (@sprintf "%#.2g" concentration)*" µM"
	else
		concString = @sprintf "%#.2g" concentration
	end
	
	lines!(axis, timePoint, surfaceTension;
		color = colors[colorIndex],
		label = concString)
	
	k = findfirst([tinds[j][4] for j in eachindex(tinds)] .== fileNames[i])
	if !isnothing(k)
		scatter!(axis, tinds[k][2], tinds[k][3];
			color = colors[colorIndex],
			strokewidth = 1.5,
			strokecolor = RGBAf(0.0, 0.0, 0.0, 1.0),
			marker = :cross,
			label = nothing)
	end
	
	k = findfirst([yeqs[j][4] for j in eachindex(yeqs)] .== fileNames[i])
	if !isnothing(k)
		scatter!(axis, yeqs[k][2], yeqs[k][3];
			color = colors[colorIndex],
			strokewidth = 2,
			strokecolor = RGBAf(0.0, 0.0, 0.0, 1.0),
			marker = :circle,
			label = nothing)
	end
	
end

axislegend(axis;
	titlefont = :regular,
	titlesize = 15,
	labeljustification = :left,
	position = (1.5, 0.5),
	backgroundcolor = :transparent,
	framevisible = false,
	unique = true,
	labelsize = 20)

let i = 11
	
	(concentration, timePoint, surfaceTension) = read_data_from_file("data/"*fileNames[i])
	colorIndex = findfirst(uniqueConcentrations .== concentration)
	
	text!(400, 72.0,
		text = "1",
		#font = "Arial bold",
		fontsize = 20,
		#color = RGBAf(0.0, 0.0, 0.0, 1.0))
		color = colors[colorIndex])
	
	text!(1800, 67.0,
		text = "2",
		#font = "Arial bold",
		fontsize = 20,
		#color = RGBAf(0.0, 0.0, 0.0, 1.0))
		color = colors[colorIndex])
	
	text!(3700, 61.75,
		text = "3",
		#font = "Arial bold",
		fontsize = 20,
		#color = RGBAf(0.0, 0.0, 0.0, 1.0))
		color = colors[colorIndex])
	
end

save("output/overview.pdf", figure)

###

plotTheme = Theme(
	size = (600, 400),
	figure_padding = (10, 150, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 2,
	markersize = 15,
	Axis = (
		xlabel = rich("Time / s"),
		ylabel = rich("Surface tension / mJ⋅m", superscript("-2")),
		limits = ((0, 50000), (45, 75)),
		xticks = 0:20000:50000,
		yticks = 45:5:75,
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
	
	if concentration == uniqueConcentrations[1]
		concString = (@sprintf "%#.2g" concentration)*" µM"
	else
		concString = @sprintf "%#.2g" concentration
	end
	
	lines!(axis, timePoint, surfaceTension;
		color = colors[colorIndex],
		label = concString)
	
	k = findfirst([tinds[j][4] for j in eachindex(tinds)] .== fileNames[i])
	if !isnothing(k)
		scatter!(axis, tinds[k][2], tinds[k][3];
			color = colors[colorIndex],
			strokewidth = 1.5,
			strokecolor = RGBAf(0.0, 0.0, 0.0, 1.0),
			marker = :cross,
			label = nothing)
	end
	
	k = findfirst([yeqs[j][4] for j in eachindex(yeqs)] .== fileNames[i])
	if !isnothing(k)
		scatter!(axis, yeqs[k][2], yeqs[k][3];
			color = colors[colorIndex],
			strokewidth = 2,
			strokecolor = RGBAf(0.0, 0.0, 0.0, 1.0),
			marker = :circle,
			label = nothing)
	end
	
end

axislegend(axis;
	titlefont = :regular,
	titlesize = 15,
	labeljustification = :left,
	position = (1.6, 0.5),
	backgroundcolor = :transparent,
	framevisible = false,
	unique = true,
	labelsize = 20)

save("output/overview_long.pdf", figure)

###

plotTheme = Theme(
	size = (200, 150),
	figure_padding = (10, 20, 10, 10), # left, right, bottom, top
	fontsize = 20,
	linewidth = 2,
	markersize = 15,
	Axis = (
		xlabel = "",
		ylabel = "",
		limits = ((0, 400), (52.5, 75)),
		xticks = 0:200:400,
		yticks = 50:10:75,
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
	
	if concentration == uniqueConcentrations[1]
		concString = (@sprintf "%#.2g" concentration)*" µM"
	else
		concString = @sprintf "%#.2g" concentration
	end
	
	lines!(axis, timePoint, surfaceTension;
		color = colors[colorIndex],
		label = concString)
	
	k = findfirst([tinds[j][4] for j in eachindex(tinds)] .== fileNames[i])
	if !isnothing(k)
		scatter!(axis, tinds[k][2], tinds[k][3];
			color = colors[colorIndex],
			strokewidth = 1.5,
			strokecolor = RGBAf(0.0, 0.0, 0.0, 1.0),
			marker = :cross,
			label = nothing)
	end
	
	k = findfirst([yeqs[j][4] for j in eachindex(yeqs)] .== fileNames[i])
	if !isnothing(k)
		scatter!(axis, yeqs[k][2], yeqs[k][3];
			color = colors[colorIndex],
			strokewidth = 2,
			strokecolor = RGBAf(0.0, 0.0, 0.0, 1.0),
			marker = :circle,
			label = nothing)
	end
	
end

save("output/overview_inset.pdf", figure)
