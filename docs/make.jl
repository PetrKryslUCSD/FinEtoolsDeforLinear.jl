using Documenter, FinEtools, FinEtoolsDeforLinear

makedocs(
	modules = [FinEtoolsDeforLinear],
	doctest = false, clean = true,
	format = Documenter.HTML(prettyurls = false),
	authors = "Petr Krysl",
	sitename = "FinEtoolsDeforLinear.jl",
	pages = Any[
		"Home" => "index.md",
		"Tutorials" => "tutorials/tutorials.md",
		"How to guide" => "guide/guide.md",
		"Types and Functions" => Any[
			"man/types.md",
			"man/functions.md"]
			]
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl.git",
)
