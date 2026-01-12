"""
    mutate_vec_pdop(vec, best_cov, Ctarget; p_move=0.4, p_add_max=0.3, p_rem_max=0.2, Nmax=30)

Applique une mutation à une configuration orbitale.
La mutation peut déplacer, ajouter ou retirer des satellites.

# Arguments
- vec       : Vecteur de taille P indiquant le nombre de satellites dans chaque plan orbital.
- best_cov  : Meilleur coverage courant (en %).
- Ctarget   : Coverage cible (en %).

# Paramètres optionnels
- p_move    : Probabilité de déplacer un satellite d'un plan vers un autre (répété P fois).
- p_add_max : Probabilité maximale d'ajouter un satellite.
- p_rem_max : Probabilité maximale d'enlever un satellite.
- Nmax      : Nombre maximal de satellites souhaités.

# Valeur retournée
- v         : Nouveau vecteur muté, basé sur `vec`.
"""
function mutate_vec_pdop(vec, best_cov, Ctarget; p_move=0.4, p_add_max=0.3, p_rem_max=0.2, Nmax=30)
    P = length(vec)
    v = copy(vec)
    N = sum(v)

    gap = (best_cov - Ctarget) / Ctarget
    p_add = p_add_max * clamp(-gap, 0.0, 1.0)  # si cov < Ctarget, on favorise légèrement l'ajout
    p_rem = p_rem_max * clamp(gap, 0.0, 1.0)   # si cov > Ctarget, on favorise légèrement le retrait

    for _ in 1:P
        if rand() < p_move && N > 1
            i = rand(1:P)
            tries = 0
            while v[i] == 0 && tries < P # On prend que les plans vides
                i = rand(1:P)
                tries += 1
            end
            v[i] == 0 && continue # On continue si on a pas trouvé de plan non-vide
            j = rand(1:P)
            i == j && continue
            v[i] -= 1 # Retrait
            v[j] += 1 # Ajout
        end
    end

    if rand() < p_add && N < Nmax  # Ajout d'un satellite sur un plan au hasard
        k = rand(1:P)
        v[k] += 1
        N += 1
    end

    if rand() < p_rem && N > 1 # Retrait d'un satellite sur un plan au hasard
        k = rand(1:P)
        tries = 0
        while v[k] == 0 && tries < P
            k = rand(1:P)
            tries += 1
        end
        v[k] == 0 && return v # On ne retire pas si on a pas trouvé de plan non-vide
        v[k] -= 1
    end

    return v
end

##########################################
##########################################
struct FitCachePDOP{K,V}
    threads::Vector{Dict{K,V}}
end

FitCachePDOP{K,V}() where {K,V} = FitCachePDOP{K,V}([Dict{K,V}() for _ in 1:Threads.nthreads()])

@inline _cache_threads_pdop(fc::FitCachePDOP, tid::Int) = fc.threads[tid == 0 ? Threads.threadid() : clamp(tid, 1, length(fc.threads))]

function Base.empty!(fc::FitCachePDOP)
    for c in fc.threads
        empty!(c)
    end
    fc
end

const FITCACHE_pdop = FitCachePDOP{Tuple{Vararg{Int}}, Tuple{Float64,Float64,Float64}}()
##########################################
##########################################

"""
    fitness_pdop(vec, F, i_deg, a, eps_deg, grid_ga; tid=0, n=10, min_sats=5,
                 pdop_weight=100.0, cov_pow=2.0, cov_bonus=5.0,
                 Ncoef=0.05, Pcoef=0.3, Ctarget=99.0, Kcov=2.0)

Évalue la qualité d'une configuration orbitale en priorisant la performance PDOP.

# Arguments
- vec      : Vecteur de taille P indiquant le nombre de satellites dans chaque plan orbital.
- F        : Paramètre de phasage (Walker-Delta) utilisé par `myconstellation`.
- i_deg    : Inclinaison orbitale en degrés.
- a        : Demi-grand axe de l'orbite (mètres).
- eps_deg  : Angle d'élévation minimal requis pour la visibilité.
- grid_ga  : Stucture `GroundGrid` contenant la grille lat/lon et des valeurs précalculées.

# Paramètres optionnels
- tid         : ID du thread si `fitness_pdop` est utilisé dans une boucle `Threads.@threads`.
- n           : Nombre de pas de temps (évaluation grossière).
- min_sats    : Nombre minimal de satellites visibles pour définir le PDOP (coverage = >= min_sats).
- pdop_weight : Poids principal associé au PDOP.
- cov_pow     : Puissance appliquée à (cov/100) dans le terme PDOP.
- cov_bonus   : Bonus additionnel lié au coverage.
- Ncoef       : Pénalité légère associée au nombre total de satellites.
- Pcoef       : Bonus léger associé aux plans orbitaux vides.
- Ctarget     : Coverage cible (en %).
- Kcov        : Pénalité appliquée si cov < Ctarget.

# Valeurs retournées
- mpdop    : PDOP moyen (spatio-temporel) sur les points où il est défini.
- cov      : Coverage moyen (%) associé au PDOP (>= min_sats visibles).
- fit      : Score global de qualité (à maximiser).
"""
function fitness_pdop(vec, F, i_deg, a, eps_deg, grid_ga; tid=0, n=10, min_sats=5,
                      pdop_weight=100.0, cov_pow=2.0, cov_bonus=5.0,
                      Ncoef=0.05, Pcoef=0.3, Ctarget=99.0, Kcov=2.0)

    key = Tuple(vec)
    cache = _cache_threads_pdop(FITCACHE_pdop, tid)

    if haskey(cache, key)
        return cache[key]
    end
    # Evaluation grossière de la configuration (ok pour l'algorithme)
    mpdop, cov, N = eval_constellation_pdop(vec, F, i_deg, a, eps_deg, grid_ga; n=n, min_sats=min_sats)
    P_empty = length(vec) - count(!iszero, vec) # Nombre de plans vides

    covf = cov / 100.0
    pdop_term = isfinite(mpdop) ? pdop_weight * (covf^cov_pow) / (mpdop + 1e-9) : 0.0
    penalty_cov = max(0.0, Ctarget - cov) * Kcov

    fit = pdop_term + cov_bonus * covf - Ncoef * N + Pcoef * P_empty - penalty_cov

    cache[key] = mpdop, cov, fit
    return mpdop, cov, fit
end

"""
    evolve_vec_pdop(P, N_init, F, i_deg, a, eps_deg; grid_ga=0, popsize=30, generations=40, Nmax=30,
                    n=10, min_sats=5, pdop_weight=100.0, cov_pow=2.0, cov_bonus=5.0,
                    Ncoef=0.05, Pcoef=0.3, Ctarget=99.0, Kcov=2.0,
                    p_move=0.4, p_add_max=0.3, p_rem_max=0.2)

Algorithme génétique qui optimise une constellation en priorisant le PDOP.

# Arguments
- P          : Nombre total de plans orbitaux disponibles.
- N_init     : Nombre total de satellites utilisés pour initialiser la population.
- F          : Paramètre de phasage (Walker-Delta).
- i_deg      : Inclinaison orbitale en degrés.
- a          : Demi-grand axe de l'orbite (en mètres), typiquement Re + altitude.
- eps_deg    : Angle d'élévation minimal requis pour la visibilité.

# Paramètres optionnels
- grid_ga      : Stucture `GroundGrid` contenant la grille lat/lon et des valeurs précalculées.
- popsize      : Taille de la population.
- generations  : Nombre de générations.
- Nmax         : Nombre maximal de satellites souhaités.
- n            : Nombre de pas de temps (évaluation grossière).
- min_sats     : Nombre minimal de satellites visibles pour définir le PDOP (coverage = >= min_sats).
- pdop_weight  : Poids principal associé au PDOP.
- cov_pow      : Puissance appliquée à (cov/100) dans le terme PDOP.
- cov_bonus    : Bonus additionnel lié au coverage.
- Ncoef        : Pénalité légère associée au nombre total de satellites.
- Pcoef        : Bonus léger associé aux plans orbitaux vides.
- Ctarget      : Coverage cible (en %).
- Kcov         : Pénalité appliquée si cov < Ctarget.
- p_move       : Probabilité de déplacement (mutation).
- p_add_max    : Probabilité maximale d'ajout (mutation).
- p_rem_max    : Probabilité maximale de retrait (mutation).

# Valeurs retournées
- best_vec     : Meilleure répartition trouvée.
- mpdop_final  : PDOP moyen final (évaluation plus fine).
- cov_final    : Coverage moyen final (évaluation plus fine).
- N_final      : Nombre total de satellites.
"""
function evolve_vec_pdop(P, N_init, F, i_deg, a, eps_deg; grid_ga=0, popsize=30, generations=40, Nmax=30,
                         n=10, min_sats=5, pdop_weight=100.0, cov_pow=2.0, cov_bonus=5.0,
                         Ncoef=0.05, Pcoef=0.3, Ctarget=99.0, Kcov=2.0,
                         p_move=0.4, p_add_max=0.3, p_rem_max=0.2)

    if grid_ga == 0  # Initialisation de GroundGrid si il n'est pas donné
        grid_ga = GroundGrid(-i_deg, i_deg; dlat=6, dlon=6)
    end

    population = Vector{Vector{Int}}(undef, popsize)
    Threads.@threads for i in 1:popsize
        population[i] = random_vec(P, N_init) # Population initiale
    end

    best_vec = population[1]
    best_pdop, best_cov, best_fit = fitness_pdop(best_vec, F, i_deg, a, eps_deg, grid_ga; n=n, min_sats=min_sats, pdop_weight=pdop_weight, 
                                                cov_pow=cov_pow, cov_bonus=cov_bonus, Ncoef=Ncoef, Pcoef=Pcoef, Ctarget=Ctarget, Kcov=Kcov)

    fits = Vector{Float64}(undef, popsize) # Vecteur avec les valeurs de fit de la population actuelle
    covs = Vector{Float64}(undef, popsize) # Vecteur avec les valeurs de cov de la population actuelle

    for _ in 1:generations
        Threads.@threads for i in 1:popsize
            tid = Threads.threadid()
            mp, cov, fit = fitness_pdop(population[i], F, i_deg, a, eps_deg, grid_ga; n=n, min_sats=min_sats, pdop_weight=pdop_weight, 
                                                cov_pow=cov_pow, cov_bonus=cov_bonus, Ncoef=Ncoef, Pcoef=Pcoef, Ctarget=Ctarget, Kcov=Kcov)
            covs[i] = cov
            fits[i] = fit
        end

        order = sortperm(fits, rev=true)  # On ordonne la population en fonction des fitness
        elite_count = clamp(popsize ÷ 4, 1, popsize)
        elite = population[order[1:elite_count]] # On prend que les meilleurs

        if fits[order[1]] > best_fit
            best_fit = fits[order[1]]
            best_cov = covs[order[1]]
            best_vec = elite[1]
            best_pdop, _, _ = fitness_pdop(best_vec, F, i_deg, a, eps_deg, grid_ga; n=n, min_sats=min_sats, pdop_weight=pdop_weight, 
                                                cov_pow=cov_pow, cov_bonus=cov_bonus, Ncoef=Ncoef, Pcoef=Pcoef, Ctarget=Ctarget, Kcov=Kcov)
        end

        newpop = Vector{typeof(best_vec)}()
        append!(newpop, elite) # Ajout des elites dans la nouvelle population

        while length(newpop) < popsize - 1
            parent = elite[rand(1:elite_count)] # Mutation d'un parent aléatoire
            child = mutate_vec_pdop(parent, best_cov, Ctarget; p_move=p_move, p_add_max=p_add_max, p_rem_max=p_rem_max, Nmax=Nmax)
            push!(newpop, child) # Ajout de la nouvelle configuration
        end

        push!(newpop, random_vec(P, N_init)) # Ajout d'un vecteur aléatoire (avec un peu de chance, on découvre une nouvelle partie de l'espace des configurations)
        population = newpop
    end

    mpdop_final, cov_final, N_final = eval_constellation_pdop(best_vec, F, i_deg, a, eps_deg, grid_ga; n=75, min_sats=min_sats) # évaluation finale plus fine
    return best_vec, mpdop_final, cov_final, N_final
end