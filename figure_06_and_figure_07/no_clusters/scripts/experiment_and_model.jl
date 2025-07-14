using CairoMakie, KissSmoothing, Printf

function read_data_from_file(FileName::String)
	
	dataString = read(FileName, String)
	
	dataString = replace(dataString, "\r\n" => "\n")
	dataString = replace(dataString, "\r" => "\n")
	dataString = replace(dataString, "\t" => "")
	dataString = replace(dataString, " " => "")
	dataString = replace(dataString, "," => ".")
	dataString = String.(split(dataString, "\n"; keepempty = false))
	
	concentration = parse(Float64, split(dataString[2], ";")[2])
	
	dataString = dataString[4:end]
	
	timePoint = zeros(length(dataString))
	surfaceTension = zeros(length(dataString))
	
	for i in eachindex(dataString)
		timePoint[i] = parse(Float64, split(dataString[i], ";")[1])
		surfaceTension[i] = parse(Float64, split(dataString[i], ";")[2])
	end
	
	return concentration::Float64, timePoint::Vector{Float64}, surfaceTension::Vector{Float64}
	
end

function induction_time(timePoints::Vector{Float64}, surfaceTension::Vector{Float64}, y0::Float64, drop::Float64)
	
	curve = denoise(surfaceTension; factor = 1.0)[1]
	
	index = findfirst(surfaceTension .< y0-drop)
	if index == 1
		return NaN, NaN
	end
	
	t = timePoints[index]
	v = y0 - drop
	
	return t::Float64, v::Float64
	
end

function meso_equilibrium(timePoints::Vector{Float64}, surfaceTension::Vector{Float64}, slope::Float64)
	
	curve = denoise(surfaceTension; factor = 4.0)[1]
	diffCurve = (curve[2:end] .- curve[1:end-1]) ./ (timePoints[2:end] .- timePoints[1:end-1])
	
	index = findfirst(reverse(diffCurve) .< slope)
	if index == nothing || index == 1
		return NaN, NaN
	end
	index = length(diffCurve) - findfirst(reverse(diffCurve) .< slope) + 1
	
	t = timePoints[index]
	v = surfaceTension[index]
	
	return t::Float64, v::Float64
	
end

#####################################

fileNames = [
"01_DNAJB6_05_000_uM.csv";
"02_DNAJB6_05_000_uM.csv";
"03_DNAJB6_03_330_uM.csv";
"04_DNAJB6_02_220_uM.csv";
"05_DNAJB6_01_480_uM.csv";
"06_DNAJB6_00_988_uM.csv";
"07_DNAJB6_00_658_uM.csv";
"08_DNAJB6_00_658_uM.csv";
"09_DNAJB6_00_439_uM.csv";
"10_DNAJB6_00_293_uM.csv";
"11_DNAJB6_00_195_uM.csv";
"12_DNAJB6_00_130_uM.csv";
"13_DNAJB6_00_130_uM.csv";
"13_b01_DNAJB6_00_130_uM.csv";
"13_b02_DNAJB6_00_130_uM.csv";
"14_DNAJB6_00_087_uM.csv";
"15_DNAJB6_00_087_uM.csv";
"15_b01_DNAJB6_00_087_uM.csv";
"15_b02_DNAJB6_00_087_uM.csv"]

uniqueConcentrations = unique([read_data_from_file("experimental_data/"*fileName)[1] for fileName in fileNames])

modelFileNames = [
"model_05_000_uM.csv";
"model_03_330_uM.csv";
"model_02_220_uM.csv";
"model_01_480_uM.csv";
"model_00_988_uM.csv";
"model_00_658_uM.csv";
"model_00_439_uM.csv";
"model_00_293_uM.csv";
"model_00_195_uM.csv";
"model_00_130_uM.csv";
"model_00_087_uM.csv"]

modelUniqueConcentrations = unique([read_data_from_file("output/data/"*fileName)[1] for fileName in modelFileNames])

#####################################

y0 = 72.0
drop = 1.5
slope = -2.5/3600

yeqs = Vector{Tuple{Float64, Float64}}()
tinds = Vector{Tuple{Float64, Float64}}()

plotTheme = Theme(
	size = (600, 400),
	figure_padding = (10, 150, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 5,
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

# colors = cgrad(:viridis, length(uniqueConcentrations); alpha = 0.15, categorical = true)
# modelColors = cgrad(:viridis, length(uniqueConcentrations); alpha = 0.8, categorical = true)

colors = cgrad([RGBf(1/4, 0, 1/3), RGBf(1/4, 1/4, 2/3), RGBf(0, 3/4, 3/4), RGBf(0, 3/4, 0), RGBf(1, 3/4, 1/8)],
	length(uniqueConcentrations);
	alpha = 0.25,
	categorical = true,
	rev = false)

modelColors = cgrad([RGBf(1/4, 0, 1/3), RGBf(1/4, 1/4, 2/3), RGBf(0, 3/4, 3/4), RGBf(0, 3/4, 0), RGBf(1, 3/4, 1/8)],
	length(uniqueConcentrations);
	alpha = 0.8,
	categorical = true,
	rev = false)

for i in reverse(eachindex(fileNames))
	
	(concentration, timePoint, surfaceTension) = read_data_from_file("experimental_data/"*fileNames[i])
	colorIndex = findfirst(uniqueConcentrations .== concentration)
	
	lines!(axis, timePoint, surfaceTension;
		color = colors[colorIndex],
		label = nothing)
	
end

for i in reverse(eachindex(modelFileNames))
	
	(concentration, timePoint, surfaceTension) = read_data_from_file("output/data/"*modelFileNames[i])
	colorIndex = findfirst(uniqueConcentrations .== concentration)
	
	lines!(axis, timePoint, surfaceTension;
		linestyle = :dash,
		color = modelColors[colorIndex],
		linewidth = 3,
		label = nothing)
	
	if concentration == uniqueConcentrations[1]
		concString = (@sprintf "%#.2g" concentration)*" µM"
	else
		concString = @sprintf "%#.2g" concentration
	end
	
	lines!(axis, [NaN], [NaN];
		linestyle = :solid,
		color = modelColors[colorIndex],
		linewidth = 3,
		label = concString)
	
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

save("output/experiment_and_model/overview.pdf", figure)

###

plotTheme = Theme(
	size = (200, 150),
	figure_padding = (10, 20, 10, 10), # left, right, bottom, top
	fontsize = 20,
	linewidth = 5,
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

# colors = cgrad(:viridis, length(uniqueConcentrations); alpha = 0.15, categorical = true)
# modelColors = cgrad(:viridis, length(uniqueConcentrations); alpha = 0.8, categorical = true)

colors = cgrad([RGBf(1/4, 0, 1/3), RGBf(1/4, 1/4, 2/3), RGBf(0, 3/4, 3/4), RGBf(0, 3/4, 0), RGBf(1, 3/4, 1/8)],
	length(uniqueConcentrations);
	alpha = 0.25,
	categorical = true,
	rev = false)

modelColors = cgrad([RGBf(1/4, 0, 1/3), RGBf(1/4, 1/4, 2/3), RGBf(0, 3/4, 3/4), RGBf(0, 3/4, 0), RGBf(1, 3/4, 1/8)],
	length(uniqueConcentrations);
	alpha = 0.8,
	categorical = true,
	rev = false)

for i in reverse(eachindex(fileNames))
	
	(concentration, timePoint, surfaceTension) = read_data_from_file("experimental_data/"*fileNames[i])
	colorIndex = findfirst(uniqueConcentrations .== concentration)
	
	lines!(axis, timePoint, surfaceTension;
		color = colors[colorIndex],
		label = nothing)
	
	for i in reverse(eachindex(modelFileNames))
	
	(concentration, timePoint, surfaceTension) = read_data_from_file("output/data/"*modelFileNames[i])
	colorIndex = findfirst(uniqueConcentrations .== concentration)
	
	lines!(axis, timePoint, surfaceTension;
		linestyle = :dash,
		color = modelColors[colorIndex],
		linewidth = 3,
		label = nothing)
	
	end
	
end

save("output/experiment_and_model/overview_inset.pdf", figure)
