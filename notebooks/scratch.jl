### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 50509764-ac74-43e5-871d-52cd513321bc
import Pkg

# ╔═╡ 3a435419-5e38-4615-8bc2-f3b60c826c1a
Pkg.activate(joinpath(@__DIR__, ".."))

# ╔═╡ e9377074-ad32-454c-bba4-620ed9dca50f
import IceFloeTracker

# ╔═╡ 6eabbe9e-7610-4a74-8ca0-b3d3eebc6db4
IceFloeTracker.fetchdata(; output = "data")

# ╔═╡ Cell order:
# ╠═50509764-ac74-43e5-871d-52cd513321bc
# ╠═3a435419-5e38-4615-8bc2-f3b60c826c1a
# ╠═e9377074-ad32-454c-bba4-620ed9dca50f
# ╠═6eabbe9e-7610-4a74-8ca0-b3d3eebc6db4
