## Script utilisé pour faire les figures illustrant les différents phénomènes

# Pour le fait que la couverture change de façon périodique:
let
    h_km = 800
    eps_deg = 10
    i_deg = 30
    P = 2
    S = 178
    F = 0
    a = Re + h_km*1e3

	sats=walker_delta(P,S,F,i_deg,a)
	
	Ts = 0:100:86800
	covs = [coverage_fraction(sats, T, -90, 90, eps_deg) for T in Ts]
	
	plot = plot(Ts, round.(covs; digits=2),
	    xlabel = "Temps (s)",
	    ylabel = "Couverture (%)",
	    legend = false,
	    title = "Évolution de la couverture en fonction du temps\nWalker-Delta P=$(P), S=$(S), i=$(i_deg)°",
	    markersize = 3,
	    color = :blue)
	savefig(plot, "./coverage_fraction_période")
end

# Pour montrer que la couverture moyenne sur une période converge avec un certain nombre d'échantillons:
let
    h_km = 800
    eps_deg = 10
    i_deg = 30
    P = 2
    S = 5
    F = 0
    a = Re + h_km*1e3

    sats = walker_delta(P, S, F, i_deg, a)

    N = 200
    ns = 2 .* (1:N)
    vals = zeros(N)

    for k in 1:N
        vals[k] = mean_coverage_fraction(sats, -i_deg, i_deg, eps_deg; n=ns[k])
    end

    plot = plot(
        ns, vals;
        xlabel = "Nombre d'échantillons temporels",
        ylabel = "Couverture moyenne (%)",
        title = "Convergence de la couverture moyenne\nWalker-Delta P=$(P), S=$(S), i=$(i_deg)°",
        lw = 2,
        markershape = :circle,
        markerstrokewidth = 0,
        legend = false,
        grid = true,
    )
	savefig(plot,"./convergence_mean_coverage")
end

# Convergence de la couverture en fonction de l'altitude (Pmax = 6, N = 19)
let
    eps_deg = 10
    i_deg = 30
    F = 0
    htest = 500:20:1000
    cov_a = Vector{Float64}(undef, length(htest))
    Pmax = 6 
    Ntest = 19 
    for idx in eachindex(htest)
        a = htest[idx]*1e3 + Re
        empty!(FITCACHE)
        results = Vector{Tuple}(undef, 15)
        Threads.@threads for t in 1:15
            results[t] = evolve_vec(Pmax, Ntest, F, i_deg, a, eps_deg; popsize=20, generations=200, Pbonus=true)
        end
        covs = [r[2] for r in results]
        cov_a[idx] = maximum(covs)
    end
	boss = plot(
	    htest, cov_a;
	    marker = (:circle, 6, :blue, stroke(0)),
	    line = (:solid, 2, :blue),
	    legend = false,
	    grid = true,
	    title = "Convergence de la couverture en fonction de l'altitude\nPmax = 6, N = 19",
	    xlabel = "Altitude des satellites (km)",
	    ylabel = "Couverture (%)",
	    titlefont = (12, "Arial")
	)
	savefig(boss, "./convergence_cov_altitude.png")
end