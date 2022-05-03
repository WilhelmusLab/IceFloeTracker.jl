### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 5339db3d-ce14-4cda-b7e4-5ff76cf717b3
using Pkg: Pkg

# ╔═╡ 0b396bdc-0eb1-4da5-9320-e9ce8194bc18
Pkg.activate(".")

# ╔═╡ e9377074-ad32-454c-bba4-620ed9dca50f
using IceFloeTracker: IceFloeTracker

# ╔═╡ 6eabbe9e-7610-4a74-8ca0-b3d3eebc6db4
IceFloeTracker.fetchdata(; output="data")

# ╔═╡ Cell order:
# ╠═5339db3d-ce14-4cda-b7e4-5ff76cf717b3
# ╠═0b396bdc-0eb1-4da5-9320-e9ce8194bc18
# ╠═e9377074-ad32-454c-bba4-620ed9dca50f
# ╠═6eabbe9e-7610-4a74-8ca0-b3d3eebc6db4
