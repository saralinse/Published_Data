using CairoMakie

include("data.jl")

plotTheme = Theme(
	size = (400, 400),
	figure_padding = (30, 20, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 2,
	markersize = 15,
	Axis = (
		xlabel = rich("[JB6] / μM"),
		ylabel = rich("Surface tension / mJ⋅m", superscript("-2")),
		limits = ((0.05, 10), (52.5, 75)),
		xscale = log10,
		xticks = LogTicks(-1:2),
		yticks = 50:5:75,
		xticksmirrored = true,
		yticksmirrored = true,
		xminorticksvisible = true,
		xminorticks = IntervalsBetween(9),
		yminorticksvisible = true,
		spinewidth = 2,
		xgridvisible = false,
		ygridvisible = false)
	)

set_theme!(plotTheme)

figure = Figure()
axis = Axis(figure[1, 1])

color = RGBAf(0.0, 0.0, 0.0, 0.75)
scatter!(axis, [(yeqs[i][1], yeqs[i][3]) for i in eachindex(yeqs)];
	color = color,
	label = nothing)

axislegend(axis;
	titlefont = :regular,
	titlesize = 20,
	labeljustification = :left,
	position = :rt,
	backgroundcolor = :transparent,
	framevisible = false,
	unique = true,
	labelsize = 20)

save("output/semi_stability.pdf", figure)
