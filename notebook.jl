### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ e6386ab2-386e-46b1-ba7e-f24f7cbd81ed
using Pkg

# ╔═╡ 4a0552ff-24a8-4567-aab3-e818f370a3cb
Pkg.activate("C://Users//habets.a/Documents//Projet ES313")

# ╔═╡ eb2c4800-ba16-11f0-1507-879cb1a882eb
# ╠═╡ show_logs = false
using OrbitalTrajectories

# ╔═╡ 721c315d-b2a6-43f7-8baf-1738a221a8b6
begin
	using DifferentialEquations
	using Plots
	using Unitful
	using DiffEqUncertainty
	using Distributions
end

# ╔═╡ 374c6557-378e-4a6e-b252-d0cc839840c8
using ModelingToolkit

# ╔═╡ 4df1a8bc-3d71-4c92-99ca-33180d34484d
state = State(
	ER3BP(:Earth, :Moon), # Pre-compiled ER3BP model
	SynodicFrame(), # Reference frame (can be omitted)
	[0.76710535, 0., 0., 0., 0.47262724, 0.], # Initial state, (x, y, z,x˙, y˙, z˙)
	(1.05π, 3π) # Propagation timespan, (t0 → tf)
)
# Propagate trajectory (default: Vern7 solver, 1.0×10−10 tolerance)

# ╔═╡ 7491947b-88a1-44e5-bacd-96a9d4dd091e
trajectory = solve(state)

# ╔═╡ 92d3c95a-765e-4f68-af7b-c3982790a443
plot(trajectory, SynodicFrame())

# ╔═╡ 5fd9711a-77c9-40d8-a4d4-376f61a8c3c0
let
	system = CR3BP(:jupiter, :europa)
	 u0 = [1.0486505808029702, 0., 0., 0., -0.09354853663949217, 0.]
	state = State(system, u0, (0., 37.))
	trajectory = solve(state)([33., 37.])

	# Define some Normal distributions for the initial values
	p_SD = uconvert(NoUnits, 100u"km" / system.props.L)
	v_SD = uconvert(NoUnits, 30u"m/s" / system.props.V)
	N_pos = [truncated(Normal(u, p_SD), u-3*p_SD, u+3*p_SD) for u in trajectory.u[1][1:2]]
	N_vel = [truncated(Normal(u, v_SD), u-3*v_SD, u+3*v_SD) for u in trajectory.u[1][4:5]]

	 prob_func(prob, _, _) = remake(prob, u0=[rand.(N_pos)..., prob.u0[3], rand.(N_vel)..., prob.u0[6]])
	
	callback = collision(system, secondary_body(system))
	ensemble = EnsembleProblem(State(trajectory); prob_func)
	trajectories = solve(ensemble, OrbitalTrajectories.Dynamics.DEFAULT_ALG,
	EnsembleThreads(); trajectories=200, callback)

	p1 = plot(trajectories; alpha=0.1, arrow=true,
		color_palette=:linear_blue_5_95_c73_n256, nomodel=true)
	p1 = plot!(p1, trajectory, size=(500,500); c=:black, lw=3)

	collision_prob_MC = cumsum(map(crashed, trajectories)) ./ (1:length(trajectories))

	u0_uncertain = [N_pos..., trajectory.u[1][3], N_vel..., trajectory.u[1][6]]

	collision_prob_Koopman = expectation(crashed, State(trajectory).prob, u0_uncertain, State(trajectory).p, Koopman(), OrbitalTrajectories.Dynamics.DEFAULT_ALG; callback)

	p2 = plot(collision_prob_MC; label="Monte Carlo expectation", legend=:bottomright, xlabel="Monte Carlo simulations", ylabel="Collision probability (%)")
	p2 = hline!(p2, [collision_prob_Koopman.u]; label="Koopman expectation")

	plot(p1,p2)
end

# ╔═╡ 93263eae-7edc-4126-9dde-3cb82024b5b5
let
    # Constantes physiques de la Terre
    μ = 3.986004418e14           # m^3/s^2, paramètre gravitationnel de la Terre
    R_earth = 6.371e6            # m, rayon moyen de la Terre

    # Paramètres orbitaux pour une orbite circulaire basse (LEO)
    altitude = 400e3             # m
    a = R_earth + altitude       # demi-grand axe
    e = 0.0                      # excentricité
    i = 0.0                      # inclinaison (orbite équatoriale)
    Ω = 0.0                      # longitude du nœud ascendant
    ω = 0.0                      # argument du périgée
    ν0 = 0.0                     # anomalie vraie initiale

    # Création de l'orbite
    orb = Orbit(a, e, i, Ω, ω, ν0, μ)

    # Échantillonnage de la trajectoire
    ν = range(0, 2π; length=500)
    traj = [position(orb, νi) for νi in ν]

    # Tracé de l'orbite
    x = [r[1] for r in traj]
    y = [r[2] for r in traj]
    plot(x, y,
        aspect_ratio = :equal,
        xlabel = "x (m)",
        ylabel = "y (m)",
        title = "Orbite circulaire autour de la Terre",
        legend = false)
end

# ╔═╡ cf4aa9e3-48cc-482c-b4f7-6509a006fda1
let
    # Paramètres du système Terre-Lune
    μ = 0.012150585609624  # masse réduite Terre-Lune

    # Variables symboliques
    @parameters t
    @variables x(t) y(t) vx(t) vy(t)
    D = Differential(t)

    # Équations du CR3BP
    r1 = sqrt((x + μ)^2 + y^2)
    r2 = sqrt((x - (1 - μ))^2 + y^2)

    ax = x - (1 - μ)*(x + μ)/r1^3 - μ*(x - (1 - μ))/r2^3 + 2*vy
    ay = y - (1 - μ)*y/r1^3 - μ*y/r2^3 - 2*vx

    eqs = [
        D(x) ~ vx,
        D(y) ~ vy,
        D(vx) ~ ax,
        D(vy) ~ ay
    ]

    # Système symbolique
    sys = ODESystem(eqs)

    # Conditions initiales (proche de Lagrange L4)
    u0 = [
        x => 0.48785,
        y => 0.86603,
        vx => 0.0,
        vy => 0.0
    ]

    # Intervalle de temps
    tspan = (0.0, 10.0)

    # Construction du problème
    prob = ODEProblem(sys, u0, tspan)
    sol = solve(prob, Vern9(), abstol=1e-10, reltol=1e-10)

    # Tracé de la trajectoire
    plot(sol, vars=(x, y), xlabel="x", ylabel="y", title="Trajectoire CR3BP Terre-Lune", aspect_ratio=:equal)
end


# ╔═╡ Cell order:
# ╠═e6386ab2-386e-46b1-ba7e-f24f7cbd81ed
# ╠═4a0552ff-24a8-4567-aab3-e818f370a3cb
# ╠═eb2c4800-ba16-11f0-1507-879cb1a882eb
# ╠═721c315d-b2a6-43f7-8baf-1738a221a8b6
# ╠═374c6557-378e-4a6e-b252-d0cc839840c8
# ╠═4df1a8bc-3d71-4c92-99ca-33180d34484d
# ╠═7491947b-88a1-44e5-bacd-96a9d4dd091e
# ╠═92d3c95a-765e-4f68-af7b-c3982790a443
# ╠═5fd9711a-77c9-40d8-a4d4-376f61a8c3c0
# ╠═93263eae-7edc-4126-9dde-3cb82024b5b5
# ╠═cf4aa9e3-48cc-482c-b4f7-6509a006fda1
