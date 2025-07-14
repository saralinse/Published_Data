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
"02_DNAJB6_05_000_uM_compression_corrected.csv";
"03_DNAJB6_00_439_uM_exchange_expansion_compression_corrected.csv"]

uniqueConcentrations = unique([read_data_from_file("data/"*fileName)[1] for fileName in fileNames])

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

(concentration, timePoint, surfaceTension) = read_data_from_file("data/"*fileNames[3])
colorIndex = findfirst(uniqueConcentrations .== concentration)
lines!(axis, timePoint, surfaceTension;
	color = RGBAf(0.0, 0.0, 1.0, 1.0),
	label = nothing)

axislegend(axis;
	titlefont = :regular,
	titlesize = 15,
	labeljustification = :left,
	position = :rt,
	backgroundcolor = :transparent,
	framevisible = false,
	unique = true,
	labelsize = 20)

text!(1000, 75.0,
	text = rich("A", subscript("r")," = A/A", subscript("0")," ="),
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(1000, 64.75,
	text = "1.0",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(1600, 69.75,
	text = "1.2",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(2400, 71.5,
	text = "1.5",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(3100, 66.0,
	text = "1.2",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(3700, 62.75,
	text = "1.0",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

text!(5000, 60.0,
	text = "0.75",
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 1.0, 1.0))

save("output/surface_tension_time_exchange.pdf", figure)
