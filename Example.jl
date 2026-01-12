### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ e53f75f8-22ce-4720-a6fe-6f344b42e979
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 3270ef27-7447-47d4-9b5e-4fdc0ff97efd
using LinearAlgebra, Plots, FileIO, PlutoUI

# ╔═╡ 9bf54eaa-adea-4422-a6a9-b1e81fe8f720
using Base.Threads # Pour la performance

# ╔═╡ 4e855cd5-1b07-4b3d-8ad9-f0212260171c
begin
    include("./OrbitalConstellations/OrbitalConstellations.jl") # Module dev pour le projet
    using .OrbitalConstellations
end

# ╔═╡ 6e81dda0-ec0d-11f0-30db-a77c261949e4
html"""
 <! -- this adapts the width of the cells to display its being used on -->
<style>
	main {
		margin: 0 auto;
		max-width: 2000px;
    	padding-left: max(160px, 10%);
    	padding-right: max(160px, 10%);
	}
</style>
"""

# ╔═╡ d8440cd3-173e-4851-9142-909b9a18a470
FITCACHE_pdop

# ╔═╡ b856ac29-3431-47cf-8d5f-8672e68dc0bb
empty!(FITCACHE_pdop)

# ╔═╡ 6427fe21-9de3-47ea-aaea-c6da4a67ff20
begin
	iter = 10
	a_galileo = 23222*1e3+Re
	i_deg = 65
	eps_deg = 10.0
	Pmax = 8
	Ninit = 12
	grid_ga = GroundGrid(-i_deg, i_deg; dlat=3, dlon=3) # Initialisation de GroundGrid
	configs = Dict()
	for _ in 1:iter
		vec_pdop, mpdop, cov, N = evolve_vec_pdop(Pmax, Ninit, 1, i_deg, Re+24000e3, eps_deg)
		if !haskey(configs,(vec_pdop, mpdop, cov, N))
			configs[(vec_pdop, mpdop, cov, N)] = 1/iter * 100
		else 
			configs[(vec_pdop, mpdop, cov, N)] += 1/iter * 100
		end
	end
	configs = sort!(collect(configs), by = x -> x[1][2], rev=true)
end

# ╔═╡ 6643807e-c239-4713-b13e-8b7472b30a73
begin
    Ns = [c[1][4] for c in configs]
    w2  = [c[2] for c in configs]

    minN, maxN = extrema(Ns)
    ticksN = collect(minN:maxN)

    histogram(
        Ns;
        weights = w2,
        bins = (minN - 0.5):1:(maxN + 0.5),
        xticks = ticksN,
        xlims = (minN - 0.5, maxN + 0.5),
        xlabel = "Nombre de satellites N",
        ylabel = "Pourcentage [%]",
        title = "Histogramme de N ($iter itérations)",
        legend = false
    )
end

# ╔═╡ 929b36c0-0b12-4c5b-bbd1-2ec49d774bb8
begin
    mpdops_all = [c[1][2] for c in configs]
    w_all      = [c[2]    for c in configs]

    mask = isfinite.(mpdops_all)
    mpdops = mpdops_all[mask]
    w_mpdop = w_all[mask]
    w_mpdop .= 100 .* w_mpdop ./ sum(w_mpdop)

    histogram(
        mpdops;
        weights = w_mpdop,
        bins = 30,
        xlabel = "mpdop",
        ylabel = "Pourcentage [%]",
        title = "Histogramme de PDOP ($iter itérations)",
        legend = false
    )
end

# ╔═╡ c9caef05-beaa-4d48-a579-d706856c142f
begin
    covs = [c[1][3] for c in configs]
    w    = [c[2]    for c in configs]

    histogram(
        covs;
        weights = w,
        bins = 95:0.25:100,
        xticks = (95:1:100, string.(95:1:100)),
        xlims = (95, 100),
        xlabel = "Couverture cov [%]",
        ylabel = "Pourcentage [%]",
        title = "Histogramme de cov ($iter itérations)",
        legend = false,
        xminorgrid = true,
        minorgrid = true,
        xminorgridalpha = 0.2,
        minorgridalpha = 0.2,
        xminorgridlinewidth = 0.4,
        minorgridlinewidth = 0.4,
        xtickfontsize = 12,
        xguidefontsize = 12,
        yguidefontsize = 12,
        titlefontsize = 13,
        tickfontsize = 11
    )
end

# ╔═╡ b1369aca-0b51-4c68-a2d4-9d93f2249003
vec_pdop, mpdop, cov, N = evolve_vec_pdop(8, 12, 1, 65, Re+24000e3, 10.0)

# ╔═╡ 64591d68-3623-4a29-8d76-dd24b8852120
begin
	sats = myconstellation([7 0 1 0 0 6 6 0], 1 , 55, Re+24000e3)
	show_pdop_coverage_heatmap(sats, 0.0, 10.0; min_sats=5, pdop_max=7)
end

# ╔═╡ 8b352526-29cc-4f4d-a7c5-ff67c222fe33
plot_constellation(myconstellation(vec_pdop, 1, 55, Re+24000e3), 0.0)

# ╔═╡ d885b095-fd41-4694-8b6d-b07ecf49eb02
# Pour faire des GIFs sur CDN :)
if !isfile("./Figures/sats24000.gif")
    folder = mktempdir()
	vec0 = [7 0 1 0 0 6 6 0]
	satsgif = myconstellation(vec0, 1, 55, Re+24000e3)
	Tmax = 10000
	step0 = 100
    for t in 0:step0:Tmax
		plot_constellation(satsgif,t)
        savefig(joinpath(folder, "frame_$(Int(t/step0)).png"))
    end
    frames = [load(joinpath(folder, "frame_$i.png")) for i in 0:Int(Tmax/step0)]
    gr()
    save("./Figures/sats24000.gif", cat(frames..., dims=3))
end

# ╔═╡ f5f3d0c6-36c4-43ff-835e-2a7b3ac56395
md"""
Click here to reload the GIF : $(@bind reload Button("Reload"))
"""

# ╔═╡ 9d902edc-39b8-4244-ab13-49e56ab79390
begin
	reload
	PlutoUI.LocalResource("./Figures/sats24000.gif")
end

# ╔═╡ Cell order:
# ╟─6e81dda0-ec0d-11f0-30db-a77c261949e4
# ╟─e53f75f8-22ce-4720-a6fe-6f344b42e979
# ╠═3270ef27-7447-47d4-9b5e-4fdc0ff97efd
# ╠═9bf54eaa-adea-4422-a6a9-b1e81fe8f720
# ╠═4e855cd5-1b07-4b3d-8ad9-f0212260171c
# ╠═d8440cd3-173e-4851-9142-909b9a18a470
# ╠═b856ac29-3431-47cf-8d5f-8672e68dc0bb
# ╠═6427fe21-9de3-47ea-aaea-c6da4a67ff20
# ╟─6643807e-c239-4713-b13e-8b7472b30a73
# ╠═929b36c0-0b12-4c5b-bbd1-2ec49d774bb8
# ╟─c9caef05-beaa-4d48-a579-d706856c142f
# ╠═b1369aca-0b51-4c68-a2d4-9d93f2249003
# ╠═64591d68-3623-4a29-8d76-dd24b8852120
# ╠═8b352526-29cc-4f4d-a7c5-ff67c222fe33
# ╠═d885b095-fd41-4694-8b6d-b07ecf49eb02
# ╟─f5f3d0c6-36c4-43ff-835e-2a7b3ac56395
# ╟─9d902edc-39b8-4244-ab13-49e56ab79390
