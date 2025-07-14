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
"04_DNAJB6_00_439_uM_expansion_corrected.csv"]

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
		limits = ((0, 10000), (47.5, 85)),
		xticks = 0:2000:10000,
		yticks = 50:5:85,
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
	#position = (1.6, 0.5),
	position = :rt,
	backgroundcolor = :transparent,
	framevisible = false,
	unique = true,
	labelsize = 20)

text!(1000, 70.0,
	text = rich("A", subscript("r")," = A/A", subscript("0")," ="),
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(1000, 63.0,
	text = "1.0",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(2100, 67.5,
	text = "1.2",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(2600, 65.75,
	text = "1.4",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(3100, 64.0,
	text = "1.6",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(4700, 60.0,
	text = "1.8",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

save("output/surface_tension_time.pdf", figure)
