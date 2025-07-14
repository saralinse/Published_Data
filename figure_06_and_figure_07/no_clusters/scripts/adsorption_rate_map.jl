using CairoMakie

include("adsorption_model.jl")

################################

cs = 10 .^ (-2:0.01:1) # [uM]
gs = collect(0:500.0:200000) # [1/um^2]

################################

plotTheme = Theme(
	size = (600, 400),
	figure_padding = (10, 20, 20, 20), # left, right, bottom, top
	fontsize = 20,
	linewidth = 2,
	markersize = 15,
	Axis = (
		xlabel = rich("Concentration, c / μM"),
		ylabel = rich("Surface excess, Γ / 10", superscript("3"), "⋅μm", superscript("-2")),
		limits = ((0.05, 10), (0, 200)),
		xscale = log10,
		xticks = LogTicks(-2:2),
		xminorticks = IntervalsBetween(9),
		yticks = 0:50:200,
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

p = heatmap!(axis, cs, 1E-3*gs, adsorption_rate_constant.(1, gs') .* (1E-6*6.022E23*1E-15*cs);
	colormap = Reverse(:blues),
	colorrange = (1E-2, 1E1),
	colorscale = log10)

Colorbar(figure[1, 2], p;
	ticks = LogTicks(-6:1:6),
	minorticksvisible = true,
	minorticks = IntervalsBetween(9),
	label = rich("Adsorption rate, dΓ/dt / s", superscript("-1"), "·μm", superscript("-2")))

text!(0.06, 163,
	text = rich("B = 70·10", superscript("3"), " s", superscript("-1"), "·μm", superscript("-2"), "\nq = 0.15"),
	fontsize = 20,
	color = RGBAf(1.0, 1.0, 1.0, 1.0))

save("output/model_adsorption_rate_map.pdf", figure)
