### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ b27148f2-0395-478c-bbff-ffa8deb15c51
import Pkg

# ╔═╡ dbc0f1b2-3525-4cc2-bb13-ff33e29cef53
Pkg.activate(".")

# ╔═╡ e9377074-ad32-454c-bba4-620ed9dca50f
import IceFloeTracker

# ╔═╡ 6eabbe9e-7610-4a74-8ca0-b3d3eebc6db4
IceFloeTracker.fetchdata(;output="data")

# ╔═╡ Cell order:
# ╠═b27148f2-0395-478c-bbff-ffa8deb15c51
# ╠═dbc0f1b2-3525-4cc2-bb13-ff33e29cef53
# ╠═e9377074-ad32-454c-bba4-620ed9dca50f
# ╠═6eabbe9e-7610-4a74-8ca0-b3d3eebc6db4
