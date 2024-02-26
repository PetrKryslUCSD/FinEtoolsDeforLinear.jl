using Documenter, FinEtools, FinEtoolsDeforLinear

makedocs(
	modules = [FinEtoolsDeforLinear],
	doctest = false, clean = true,
	warnonly = Documenter.except(:linkcheck, :footnote),
	format = Documenter.HTML(prettyurls = false),
	authors = "Petr Krysl",
	sitename = "FinEtoolsDeforLinear.jl",
	pages = Any[
		"Home" => "index.md",
		"How to guide" => "guide/guide.md",
		"Reference" => Any[
			"man/man.md"]
			]
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/FinEtoolsDeforLinear.jl.git",
)
