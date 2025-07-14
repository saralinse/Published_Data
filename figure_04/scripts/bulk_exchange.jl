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
"01_DNAJB6_00_439_uM_to_buffer.csv";
"02_DNAJB6_00_439_uM_to_DNAJB6_05_000_uM.csv"]

#####################################

y0 = 72.0

plotTheme = Theme(
	size = (400, 400),
	figure_padding = (10, 30, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 2,
	markersize = 15,
	Axis = (
		xlabel = rich("Time / s"),
		ylabel = rich("Surface tension / mJ⋅m", superscript("-2")),
		limits = ((0, 5000), (45, 75)),
		xticks = 0:1000:5000,
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

(concentration, timePoint, surfaceTension) = read_data_from_file("data/"*fileNames[1])
lines!(axis, timePoint, surfaceTension;
	color = RGBAf(0.0, 0.0, 1.0, 1.0),
	label = "0.44 to 0.0 µM")

(concentration, timePoint, surfaceTension) = read_data_from_file("data/"*fileNames[2])
lines!(axis, timePoint, surfaceTension;
	color = RGBAf(0.0, 1.0, 0.0, 1.0),
	label = "0.44 to 5.0 µM")

axislegend(axis;
	titlefont = :regular,
	titlesize = 15,
	labeljustification = :left,
	position = :rt,
	backgroundcolor = :transparent,
	framevisible = false,
	unique = true,
	labelsize = 20)

save("output/bulk_exchange.pdf", figure)
