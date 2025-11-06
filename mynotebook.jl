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

# ╔═╡ 49799ef0-ba51-11f0-28de-b37bf41e41e8
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 8969a0b2-50f7-4573-9c63-40dcf7ef773e
using LinearAlgebra, Plots

# ╔═╡ daca6429-e148-42d3-9499-bfd1dbc04531
begin
	using Logging
	Logging.disable_logging(LogLevel(1000));
end

# ╔═╡ 287f599c-0932-4d5c-ae28-16c6488d585a
using FileIO, PlutoUI

# ╔═╡ ffc4fe26-29ab-405c-abb4-441945e251f0
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

# ╔═╡ 488991a5-2fda-4bc2-8702-74ad4459f1d4
## Constantes

begin
	const μ = 3.986004418e14 	# Paramètre gravitationnel terrestre
	const Re = 6371e3 			# Rayon moyen de la Terre
	const ωe = 7.2921150e-5 	# Vitesse de rotation de la Terre (considérée constante)

	nothing
end

# ╔═╡ e40e03da-6dad-40b0-b739-ae42d8ed98fb
## Structure pour stocker les satellites (définis par a, i, Ω et M0).

struct Sat
    a::Float64   # Demi-grand axe (m) → distance moyenne au centre de la Terre
    i::Float64   # Inclinaison orbitale (rad) → angle entre le plan orbital et l’équateur
    Ω::Float64   # Longitude du nœud ascendant (rad) → orientation du plan orbital autour de la Terre
    M0::Float64  # Anomalie moyenne initiale (rad) → position du satellite sur son orbite
end


# ╔═╡ 2e3b8833-e011-41b1-be61-fdc9905e3115
## deg2rad et rad2deg functions

begin
	deg2rad(x)=x*pi/180
	rad2deg(x)=180*x/pi

	nothing
end

# ╔═╡ 03b47088-7ff3-4817-9012-7d156c638855
## Matrices de rotation

begin
	R1(θ)=[1 0 0; 0 cos(θ) -sin(θ); 0 sin(θ) cos(θ)] # Autour de l'axe X → Sert pour incliner le plan orbital d’un angle i (l'inclinaison).
	R3(θ)=[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] # Autour de l'axe Z → Sert pour tourner le plan orbital d’un angle Ω (ascension du nœud) ou de l’anomalie vraie ν.

	nothing
end

# ╔═╡ 326e41c9-3240-415a-a50a-135906a2e5c7
"""
Convertit un vecteur de position d’un satellite du repère ECI (Earth-Centered Inertial) vers le repère ECEF (Earth-Centered Earth-Fixed).
"""
function ecef_from_eci(r::Vector{},t::Int64)
	return R3(-ωe*t)*r # page 129 - SE216 (v.Aout 2024)
end

# ╔═╡ d60e29d5-eeb4-4d19-87a9-56a0e62a8a75
"""
Convertit un vecteur position exprimé dans le repère ECEF (Earth-Centered Earth-Fixed) en latitude et longitude (degrés).
"""
function latlon_from_ecef(r::Vector{})
    x, y, z = r # Vecteur position
    ρ = norm(r) # Distance entre le satellite et le centre de la Terre
	
    ϕ = asin(z / ρ) 	# Latitude -> angle entre le vecteur position et le plan équatorial.  
    λ = atan(y, x) 		# Longitude -> angle dans le plan équatorial entre l’axe X (Greenwich) et la projection du vecteur sur ce plan.

    return rad2deg(ϕ), rad2deg(λ)
end


# ╔═╡ 735be1be-0989-4e0d-8c36-910d80870534
"""
Calcule la position d’un satellite dans le repère inertiel (ECI) à un instant t, en supposant une orbite circulaire.
"""
function eci_pos(s,t::Int64)
    n = sqrt(μ / s.a^3) 			# Vitesse angulaire du satellite d'après la 3e loi de Kepler
    u = s.M0 + n * t 				# Position angulaire du satellite sur son orbite (Position initiale (M₀) + n * t)
    R = R3(s.Ω) * R1(s.i) * R3(u) 	# Matrice de transformation (Formule 7.35 - SE216 (v.Aout 2024))
    return R * [s.a, 0, 0] 			# Position finale du satellite
end

# ╔═╡ 48c2baeb-e62a-44a7-87d0-148705574470
"""
Initialise une constellation de type Walker-Delta.\n
Utilisée pour répartir uniformément des satellites sur plusieurs plans orbitaux inclinés.
"""
function walker_delta(P::Int64,S::Int64,F::Int64,i_deg::Int64,a::Float64)
    sats = Sat[]
    for p in 0:P-1, s in 0:S-1 # On itère sur chaque plan orbital et sur chaque satellite par plan
        
		Ω = 2π * p / P # Ascension du nœud (Ω) → Chaque plan orbital est séparé de 360°/P autour de l’axe z.
		
        # Répartition uniforme des satellites sur chaque orbite (s/S),
        M0 = 2π * (s / S + (F * p) / (S * P)) # Anomalie moyenne initiale (M₀) + déphasage (F*p)/(S*P) pour éviter que les plans soient alignés.

        push!(sats, Sat(a, deg2rad(i_deg), Ω, M0)) # Création du satellite en fonction de a, i, Ω, M₀
    end
    return sats
end

# ╔═╡ 82695929-6414-4e51-ae34-c882211530c0
"""
Vérifie si un satellite est visible depuis un point au sol donné,
en tenant compte d’un angle d’élévation minimal (ε).
"""
function visible(r_ecef::Vector{},lat_deg::Int64,lon_deg::Int64,eps_deg::Int64)
    ϕ = deg2rad(lat_deg)
    λ = deg2rad(lon_deg)

    g = Re * [cos(ϕ) * cos(λ), cos(ϕ) * sin(λ), sin(ϕ)] / Re # Vecteur qui pointe du centre de la Terre vers le point au sol défini par sa latitude ϕ et sa longitude λ. (repère ECEF)

    r = r_ecef 	# Vecteur position (ECEF)
	@show r
    ρ = norm(r) # Distance entre le satellite et le centre de la Terre

    cosγ = dot(r, g) / (ρ * 1.0) 	# Produit scalaire normalisé
    γ = acos(clamp(cosγ, -1, 1))  	# Angle entre la direction du satellite et celle du point au sol
	
    ψmax = acos((Re / ρ) * cos(deg2rad(eps_deg))) # Angle au centre de la Terre jusqu’où le satellite reste visible au-dessus d’un certain angle d’élévation ε.

    return γ ≤ ψmax # Satellite visible avec au moins l'élévation minimale spécifiée
end


# ╔═╡ b04c08c9-62ad-4039-a90c-86f9e5478c55
"""
Calcule la fraction (%) de points visibles au moins par un satellite à un instant donné.\n
Tient compte de l'élévation minimale nécessaire pour voir les satellites depuis le sol.
"""
function coverage_fraction(sats,t::Int64,latmin::Int64,latmax::Int64,eps_deg::Int64)
    pts = 0
    covered = 0

    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats] # Position des satellites à l'instant t

    for lat in latmin:2:latmax, lon in -180:2:180
        pts += 1 
        for r in r_ecef 
            if visible(r, lat, lon, eps_deg)
                covered += 1 # +1 si un satellite voit le point 
                break
            end
        end
    end
    return covered / pts # Fraction couverte
end

# ╔═╡ 4899cb2f-04a9-40a6-9f56-0a63a8a37db3
"""
Affiche les zones couvertes par les satellites.\n
Tient compte de l'élévation minimale nécessaire pour voir les satellites depuis le sol.
"""
function show_coverage_heatmap(sats,t::Int64,eps_deg::Int64)
    lats = collect(-90:1:90)
    lons = collect(-180:1:180)
    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats]

    M = falses(length(lats), length(lons))

    for (i, lat) in enumerate(lats), (j, lon) in enumerate(lons)
        M[i, j] = any(r -> visible(r, lat, lon, eps_deg), r_ecef)
    end
	p = heatmap(lons, lats, Int.(M))
	p = plot!(p, xlabel = "Longitude [°]", ylabel = "Latitude [°]", title = "Zones couvertes par $(length(sats)) satellites")
	
    return plot!(p, aspect_ratio=1, colorbar=false, framestyle=:none)
end


# ╔═╡ 08104ac0-7a6f-4f2b-ba38-43ba48b8eaae
"""
Plot une sphère (la Terre) avec des lignes pour les lattitudes et l'équateur en bleu.
"""
function plot_earth(;Rearth=Re::Float64)
	# Hémisphère arrière de la Terre (longitudes ~ [-π, 0])
	u = range(-pi/2, pi/2, 60); v1 = range(-pi, 0, 60)
	xs = [Re*cos(ui)*cos(vi) for ui in u, vi in v1]
	ys = [Re*cos(ui)*sin(vi) for ui in u, vi in v1]
	zs = [Re*sin(ui)         for ui in u, vi in v1]
	p = surface(xs, ys, zs, color=:lightblue, opacity=0.03, linecolor=:transparent)
	
	# Hémisphère avant de la Terre (longitudes ~ [0, π])
	v2 = range(0, pi, 60)
	xs = [Re*cos(ui)*cos(vi) for ui in u, vi in v2]
	ys = [Re*cos(ui)*sin(vi) for ui in u, vi in v2]
	zs = [Re*sin(ui)         for ui in u, vi in v2]
	p = surface!(p, xs, ys, zs, color=:lightblue, opacity=0.05, linecolor=:transparent)

	# Ajout de lignes pour les lattitudes
	for lat_deg in -80:10:-10
        lat_rad = deg2rad(lat_deg)
        r = Rearth * cos(lat_rad)
        z = Rearth * sin(lat_rad)
        θs = range(0, 2π, 100)
        x = [r * cos(θ) for θ in θs]
        y = [r * sin(θ) for θ in θs]
        z_line = fill(z, length(θs))
        p = plot3d!(p, x, y, z_line, lw=0.3, c=:black, alpha=0.3)
    end
	for lat_deg = 0
        lat_rad = deg2rad(lat_deg)
        r = Rearth * cos(lat_rad)
        z = Rearth * sin(lat_rad)
        θs = range(0, 2π, 100)
        x = [r * cos(θ) for θ in θs]
        y = [r * sin(θ) for θ in θs]
        z_line = fill(z, length(θs))
        p = plot3d!(p, x, y, z_line, lw=0.3, c=:blue)
    end
	for lat_deg in 10:10:80
        lat_rad = deg2rad(lat_deg)
        r = Rearth * cos(lat_rad)
        z = Rearth * sin(lat_rad)
        θs = range(0, 2π, 100)
        x = [r * cos(θ) for θ in θs]
        y = [r * sin(θ) for θ in θs]
        z_line = fill(z, length(θs))
        p = plot3d!(p, x, y, z_line, lw=0.3, c=:black, alpha=0.3)
    end
	
	return plot!(p, aspect_ratio=:equal, legend=false, colorbar=false, size=(600,600), framestyle=:none, ticks=:none, camera=(0,5))
end 

# ╔═╡ 4d1173b5-3158-4535-abdd-14da09801f09
"""
Plot les positions instantanées des satellites en 3D autour de la Terre.
"""
function plot_constellation!(sats,t::Int64; Rearth=Re::Float64)
	p = plot_earth()
	
	θs = range(0, 2π, 100)
    for s in sats
        R = R3(s.Ω) * R1(s.i)  # orientation du plan orbital
        orb = [R * [s.a * cos(θ), s.a * sin(θ), 0.0] for θ in θs]
        Xorb = [r[1] for r in orb]; Yorb = [r[2] for r in orb]; Zorb = [r[3] for r in orb]
    p = plot3d!(p, Xorb, Yorb, Zorb, lw=0.5, c=:gray, alpha=0.5)  # Tracé des plans orbitaux
    end
	
    r_ecef = [eci_pos(s, t) for s in sats] # Position des satellites
    X = [r[1] for r in r_ecef]
    Y = [r[2] for r in r_ecef]
    Z = [r[3] for r in r_ecef]
    p = scatter3d!(p, X, Y, Z, marker=:circle, ms=3.5, color=:orange, title="Visualisation des satellites en t = $(Int(t))s") # Plot des satellites

	return p
end

# ╔═╡ 6730fbd6-2cdd-4f88-a00c-182a601e6d97
## Paramètres

begin
	h_km=800 			# Altitude des satellites
	eps_deg=10 			# Elevation minimale nécéssaire pour voir le satellite depuis le sol
	i_deg=30 			# Inclinaison orbitale (angle avec l'équateur) pour savoir les lattitudes couvertes
	P=2 				# Nombre de plans orbitaux
	S=7 				# Nombre de satellites par plan
	F=0 				# Facteur de déphasage (pour décaler la position des satellites entre les plans afin d'éviter qu'ils soient alignés)
	a=Re + h_km*1e3 	# Demi-grand axe pour calculer la période orbitale
	t=0 				# Temps en secondes

	sats=walker_delta(P,S,F,i_deg,a)
end

# ╔═╡ 1332fe9d-cd3d-4c04-9009-42d2b80ddfb1
println("Couverture ±$(i_deg)° instantanée (ε=$(eps_deg)°) : $(round(coverage_fraction(sats, t, -i_deg, i_deg, eps_deg)*100,digits=2))% avec $P fois $S satellites.")

# ╔═╡ 4985fe03-5e41-4ad6-899e-d864639107f8
round(coverage_fraction(sats, t, -i_deg, i_deg, eps_deg)*100,digits=2)

# ╔═╡ 27790501-c853-4b7d-9fdc-ab78585745fa
begin
	p1 = show_coverage_heatmap(sats, t, eps_deg)
	p2 = plot_constellation!(sats, t)
	plot(p1, p2; layout=(1,2), size=(1300,600))
end

# ╔═╡ e67f0e51-c074-4ed6-b38a-fd4afbdaa9be
# r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats]

# ╔═╡ d97c7a1e-5929-4d50-aa6e-c1b029be1861
# ╠═╡ disabled = true
#=╠═╡
begin
	Ts = 0:100:10000
	covs = [coverage_fraction(sats, T, -90, 90, eps_deg)*100 for T in Ts]
	
	plot(Ts, round.(covs; digits=2),
	    xlabel = "Temps (s)",
	    ylabel = "Couverture (%)",
	    legend = false,
	    title = "Évolution de la couverture",
	    markersize = 3,
	    color = :blue)
end
  ╠═╡ =#

# ╔═╡ 84ff8c5f-5c92-4336-acb4-cd643d3e56e6
# Pour faire des GIFs sur CDN :)
if !isfile("./cov.gif") #!isfile("./sats.gif")
    folder = mktempdir()
	Tmax = 10000
	step = 100
    for t in 0:step:Tmax
		show_coverage_heatmap(sats, t, eps_deg)
		#plot_constellation!(sats,t)
        savefig(joinpath(folder, "frame_$(Int(t/step)).png"))
    end
	
    frames = [load(joinpath(folder, "frame_$i.png")) for i in 0:Int(Tmax/step)]
    gr()
    save("./cov.gif", cat(frames..., dims=3))
	#save("./sats.gif", cat(frames_..., dims=3))
end

# ╔═╡ 89b95e5c-684c-44ca-9455-469e3bb97129
md"""
Click here to reload the GIF : $(@bind reload Button("Reload"))
"""

# ╔═╡ 8926723b-d835-4818-9073-89eea4b0dea4
begin
	reload
	#PlutoUI.LocalResource("./sats.gif")
end

# ╔═╡ b2388c97-997d-4afd-a681-2b86b7c1458a
begin
	reload
	#PlutoUI.LocalResource("./cov.gif")
end

# ╔═╡ Cell order:
# ╟─ffc4fe26-29ab-405c-abb4-441945e251f0
# ╟─49799ef0-ba51-11f0-28de-b37bf41e41e8
# ╠═8969a0b2-50f7-4573-9c63-40dcf7ef773e
# ╠═daca6429-e148-42d3-9499-bfd1dbc04531
# ╠═488991a5-2fda-4bc2-8702-74ad4459f1d4
# ╠═e40e03da-6dad-40b0-b739-ae42d8ed98fb
# ╠═2e3b8833-e011-41b1-be61-fdc9905e3115
# ╠═03b47088-7ff3-4817-9012-7d156c638855
# ╟─326e41c9-3240-415a-a50a-135906a2e5c7
# ╟─d60e29d5-eeb4-4d19-87a9-56a0e62a8a75
# ╟─735be1be-0989-4e0d-8c36-910d80870534
# ╟─48c2baeb-e62a-44a7-87d0-148705574470
# ╟─82695929-6414-4e51-ae34-c882211530c0
# ╟─b04c08c9-62ad-4039-a90c-86f9e5478c55
# ╟─4899cb2f-04a9-40a6-9f56-0a63a8a37db3
# ╟─08104ac0-7a6f-4f2b-ba38-43ba48b8eaae
# ╟─4d1173b5-3158-4535-abdd-14da09801f09
# ╠═6730fbd6-2cdd-4f88-a00c-182a601e6d97
# ╠═1332fe9d-cd3d-4c04-9009-42d2b80ddfb1
# ╠═4985fe03-5e41-4ad6-899e-d864639107f8
# ╠═27790501-c853-4b7d-9fdc-ab78585745fa
# ╠═e67f0e51-c074-4ed6-b38a-fd4afbdaa9be
# ╠═d97c7a1e-5929-4d50-aa6e-c1b029be1861
# ╠═287f599c-0932-4d5c-ae28-16c6488d585a
# ╠═84ff8c5f-5c92-4336-acb4-cd643d3e56e6
# ╟─89b95e5c-684c-44ca-9455-469e3bb97129
# ╠═8926723b-d835-4818-9073-89eea4b0dea4
# ╠═b2388c97-997d-4afd-a681-2b86b7c1458a
