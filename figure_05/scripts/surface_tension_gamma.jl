using CairoMakie, LsqFit

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
	
	area = zeros(length(dataString))
	surfaceTension = zeros(length(dataString))
	
	for i in eachindex(dataString)
		area[i] = parse(Float64, split(dataString[i], ";")[1])
		surfaceTension[i] = parse(Float64, split(dataString[i], ";")[2])
	end
	
	return concentration::Float64, area::Vector{Float64}, surfaceTension::Vector{Float64}
	
end

#####################################

fileNames = [
"DNAJB6_00_439_uM_area_tension.csv";
"DNAJB6_05_000_uM_area_tension.csv";
"DNAJB6_00_439_uM_exchange_area_tension.csv"]

uniqueConcentrations = unique([read_data_from_file("data/"*fileName)[1] for fileName in fileNames])

#####################################

y0 = 72.0

p4(p1, p2, p3) = log(-exp(p1*(p2+p3)) + exp(p1*y0))/p1
model(gamma, p) = log.(exp.(p[1] .* (p[2] .* gamma .+ p[3])) .+ exp.(p[1] .* p4(p[1], p[2], p[3]))) ./ p[1]
p0 = [log(1.1); -10.0; 80.0]

lower_p = [0.0; -Inf; 0.0]
upper_p = [Inf; 0.0; Inf]

(_, area1, surfaceTension1) = read_data_from_file("data/"*fileNames[1]) # 0.439 µM
(_, area2, surfaceTension2) = read_data_from_file("data/"*fileNames[2]) # 5 µM
(_, area3, surfaceTension3) = read_data_from_file("data/"*fileNames[3]) # 0.439 µM to 0 µM

gamma1 = 120.0*1.00 ./ area1
gamma2 = 120.0*1.33 ./ area2
gamma3 = 120.0*0.40 ./ area3

(gamma, surfaceTension) = ([gamma1; gamma2; gamma3], [surfaceTension1; surfaceTension2; surfaceTension3])

fit = curve_fit(model, gamma, surfaceTension, p0; lower = lower_p,  upper = upper_p)
a1 = round(exp(fit.param[1]); sigdigits = 4)
(a2, a3, a4) = round.((fit.param[2], fit.param[3], p4(fit.param[1], fit.param[2], fit.param[3])); sigdigits = 3)
curve(gamma) = log(a1^(a2*gamma + a3) + a1^a4)/log(a1)

dcurve(gamma) = a1^(a2*gamma+a3)*a2 / (a1^(a2*gamma+a3) + a1^a4)
ddcurve(gamma) = log(a1)*a1^(a2*gamma+a3+a4)*a2^2 / (a1^a4 + a1^(a2*gamma+a3))^2
curvature(gamma) = abs(ddcurve(gamma)) / (1 + dcurve(gamma)^2)^(3/2)

@show (a1, a2, a3, a4)

plotTheme = Theme(
	size = (500, 500),
	figure_padding = (10, 30, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 2,
	markersize = 15,
	Axis = (
		xlabel = rich("Relative surface excess, Γ", subscript("r"), " = Γ/Γ", subscript("0")),
		ylabel = rich("Surface tension / mJ⋅m", superscript("-2")),
		limits = ((0, 7), (45, 75)),
		xticks = 0:1.0:7,
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

lines!(axis, [0.0; 1.5], [y0; y0];
	color = RGBAf(0.0, 0.0, 0.0, 1.0),
	label = rich("γ", subscript("0")))

scatter!(axis, gamma1, surfaceTension1;
	color = RGBAf(0.0, 0.0, 1.0, 1.0),
	label = "0.44 µM")

scatter!(axis, gamma3, surfaceTension3;
	color = RGBAf(0.0, 0.0, 1.0, 0.0),
	strokecolor = RGBAf(0.0, 0.0, 1.0, 1.0),
	strokewidth = 2.0,
	label = "0.44 to 0.0 µM")

scatter!(axis, gamma2, surfaceTension2;
	color = RGBAf(0.0, 1.0, 0.0, 1.0),
	label = "5.0 µM")

gs = collect(0:0.01:7)

lines!(axis, gs, 0 .* gs .+ a4;
	linestyle = :dot,
	color = RGBAf(0.0, 0.0, 0.0, 0.4),
	label = nothing)

lines!(axis, gs, curve.(gs);
	linestyle = :dot,
	color = RGBAf(0.0, 0.0, 0.0, 1.0),
	label = rich("ln(", string(a1), superscript(string(a2), "·Γ", subscript("r"), "+", string(a3)), "+\n+", string(a1), superscript(string(a4)), ")/ln(", string(a1), ")"; fontsize = 20) )

scatter!(axis, (3.5, curve(3.5));
	marker = :xcross,
	markersize = 30,
	color = RGBAf(0.0, 0.0, 0.0, 0.4),
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

text!(3.7, 54.0,
	text = rich("Γ", subscript("collapse"), " ≈ 3.5·Γ", subscript("0")),
	fontsize = 20,
	color = RGBAf(0.0, 0.0, 0.0, 1.0))

save("output/surface_tension_gamma.pdf", figure)
