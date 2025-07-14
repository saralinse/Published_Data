using KissSmoothing

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

uniqueConcentrations = unique([read_data_from_file("data/"*fileName)[1] for fileName in fileNames])

#####################################

yeqs::Vector{Tuple{Float64, Float64, Float64, String}} = []
tinds::Vector{Tuple{Float64, Float64, Float64, String}} = []

for i in eachindex(fileNames)
	
	y0 = 72.0
	drop = 1.5
	slope = -2.5/3600
	
	(concentration, timePoint, surfaceTension) = read_data_from_file("data/"*fileNames[i])
	
	t, v = induction_time(timePoint, surfaceTension, y0, drop)
	if t > 10
		push!(tinds, (concentration, t, v, fileNames[i]))
	end
	
	t, v = meso_equilibrium(timePoint, surfaceTension, slope)
	if t < 5000
		push!(yeqs, (concentration, t, v, fileNames[i]))
	end
	
end
