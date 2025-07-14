using CairoMakie

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

#####################################

fileNames = [
"01_DNAJB6_00_439_uM_compression_corrected.csv";
"02_DNAJB6_05_000_uM_compression_corrected.csv"]

#####################################

y0 = 72.0

plotTheme = Theme(
	size = (500, 500),
	figure_padding = (10, 30, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 2,
	markersize = 15,
	Axis = (
		xlabel = rich("Time / s"),
		ylabel = rich("Surface tension / mJ⋅m", superscript("-2")),
		limits = ((0, 10000), (37.5, 75)),
		xticks = 0:2000:10000,
		yticks = 40:5:75,
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

(concentration, timePoint, surfaceTension) = read_data_from_file("data/"*fileNames[1])
lines!(axis, timePoint, surfaceTension;
	color = RGBAf(0.0, 0.0, 1.0, 1.0),
	label = "0.44 µM")

axislegend(axis;
	titlefont = :regular,
	titlesize = 15,
	labeljustification = :left,
	position = :rt,
	backgroundcolor = :transparent,
	framevisible = false,
	unique = true,
	labelsize = 20)

text!(1200, 63.0,
	text = rich("A", subscript("r")," = A/A", subscript("0")," ="),
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(1200, 60.0,
	text = "1.0",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(2200, 56.5,
	text = "0.89",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(2700, 54.75,
	text = "0.79",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(3200, 53.0,
	text = "0.68",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(4700, 53.0,
	text = "0.56",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

save("output/surface_tension_time_0_439uM.pdf", figure)

#####################################

y0 = 72.0

plotTheme = Theme(
	size = (500, 500),
	figure_padding = (10, 30, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 2,
	markersize = 15,
	Axis = (
		xlabel = rich("Time / s"),
		ylabel = rich("Surface tension / mJ⋅m", superscript("-2")),
		limits = ((0, 10000), (37.5, 75)),
		xticks = 0:2000:10000,
		yticks = 40:5:75,
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

(concentration, timePoint, surfaceTension) = read_data_from_file("data/"*fileNames[2])
lines!(axis, timePoint, surfaceTension;
	color = RGBAf(0.0, 1.0, 0.0, 1.0),
	label = "5.0 µM")

axislegend(axis;
	titlefont = :regular,
	titlesize = 15,
	labeljustification = :left,
	position = :rt,
	backgroundcolor = :transparent,
	framevisible = false,
	unique = true,
	labelsize = 20)

text!(1200, 60.0,
	text = rich("A", subscript("r")," = A/A", subscript("0")," ="),
	fontsize = 20,
	color = RGBAf(0.0, 1.0, 0.0, 1.0))

text!(1200, 57.0,
	text = "1.0",
	fontsize = 20,
	color = RGBAf(0.0, 1.0, 0.0, 1.0))

text!(2200, 56.0,
	text = "0.89",
	fontsize = 20,
	color = RGBAf(0.0, 1.0, 0.0, 1.0))

text!(2700, 54.25,
	text = "0.78",
	fontsize = 20,
	color = RGBAf(0.0, 1.0, 0.0, 1.0))

text!(3200, 52.5,
	text = "0.68",
	fontsize = 20,
	color = RGBAf(0.0, 1.0, 0.0, 1.0))

text!(4700, 52.5,
	text = "0.57",
	fontsize = 20,
	color = RGBAf(0.0, 1.0, 0.0, 1.0))

save("output/surface_tension_time_5uM.pdf", figure)
