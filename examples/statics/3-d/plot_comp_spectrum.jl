using DelimitedFiles
using PGFPlotsX

loadcsv(n) = (
	open(n) do f
	    return readdlm(f, skipstart=1)
	end) 
table(data) = Table([:x => vec(7:length(data)), :y => vec(data[7:end])])
addtoitems!(items, n, color, mark, style) = begin
@show n, color, mark, style
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
colors = vec(["red" "green" "blue" "black" "yellow"  "brown"  "teal"  "orange"  "violet"  "cyan" "magenta"  "gray"])
marks = vec(["*" "square*" "triangle*" "star" "diamond*" "otimes*" "pentagon*" "p" "a"])

items = []


for (j, basen) in enumerate(["comp_hex_spectrum_full", "comp_hex_spectrum_underintegrated", "comp_hex_spectrum_ms"])
	for (i, nu) in enumerate([0.3 0.4999 0.499999])
		tag = "nu=$(nu)"
		n = "$(basen)-$(tag).csv"
		addtoitems!(items, n, colors[j], marks[j], styles[i])
	end
end
@show items

@pgf _a = SemiLogYAxis({
	xlabel = "Mode [ND]",
	ylabel = "Generalized stiffness [N/m]",
	grid="major",
	legend_pos  = "north west"
	},
	items...)
display(_a)
