using DelimitedFiles
using PGFPlotsX

loadcsv(n) = (
	open(n) do f
	    return readdlm(f, skipstart=1)
	end) 
table(data) = Table([:x => vec(7:length(data)), :y => vec(data[7:end])])
addtoitems!(items, n, color, mark, style) = begin
	t = table(loadcsv(n))
	options = [color, "only marks", style, "mark=$mark"]
	_p = Plot(@pgf({options...}), t)
    push!(items, _p)
    options = [color, style]
    _p = Plot(@pgf({options...}), t)
    push!(items, _p)
    return items
end
styles = vec(["solid"
"densely dashed"
"dashed"
"densely dotted"
"dotted"
"loosely dotted"])
colors = vec(["red" "green" "blue" "magenta" "black" "yellow"  "brown"  "teal"  "orange"  "violet"  "cyan"  "gray"])
marks = vec(["*" "square*" "triangle*" "star" "pentagon*" "diamond*" "otimes*" "p" "a"])

items = []


# ["bend_hex_spectrum_full", "bend_hex_spectrum_underintegrated", "bend_hex_spectrum_ms", "bend_hex_spectrum_im"]
for (j, basen) in enumerate(["bend_hex_spectrum_full", "bend_hex_spectrum_im"])
	for (i, aspect) in enumerate([1.0 10.0 100.0 1000.0])
		tag = "aspect=$(aspect)"
		n = "$(basen)-$(tag).csv"
		addtoitems!(items, n, colors[j], marks[j], styles[i])
	end
end

@pgf _a = SemiLogYAxis({
	xlabel = "Mode [ND]",
	ylabel = "Generalized stiffness [N/m]",
	grid="major",
	legend_pos  = "north west"
	},
	items...)
display(_a)
