"""
    pdop_calcul(r_ecef, ρs, cosψmax_s, gx, gy, gz; min_sats=4)

Calcule le PDOP 3D à un point sol défini par le vecteur unitaire (gx, gy, gz).

# Arguments
- r_ecef      : Positions satellites en ECEF (vecteurs 3D).
- ρs          : Normes des positions satellites (‖r_sat‖).
- cosψmax_s   : Seuils de visibilité (cosinus) associés à `eps_deg`, pré-calculés par satellite.
- gx, gy, gz  : Composantes du vecteur unitaire du point sol (ECEF unitaire).

# Paramètres optionnels
- min_sats    : Nombre minimal de satellites visibles requis pour calculer le PDOP (défaut = 4).

# Valeur retournée
- PDOP 3D. Retourne `Inf` si le nombre de satellites visibles est insuffisant.
"""
function pdop_calcul(r_ecef, ρs, cosψmax_s, gx, gy, gz; min_sats=4)
    rx = Re * gx # Distance en x entre le point au sol et le centre de la terre
    ry = Re * gy # Distance en y entre le point au sol et le centre de la terre
    rz = Re * gz # Distance en z entre le point au sol et le centre de la terre

    # Initialisation des composantes de G
    g11 = 0.0; g12 = 0.0; g13 = 0.0; g14 = 0.0
    g22 = 0.0; g23 = 0.0; g24 = 0.0
    g33 = 0.0; g34 = 0.0
    g44 = 0.0

    m = 0  # nombre de satellites visibles utilisés dans H
    
    @inbounds for k in eachindex(r_ecef)
        sx, sy, sz = r_ecef[k]
        cosγ = (sx*gx + sy*gy + sz*gz) / ρs[k]
        cosγ >= cosψmax_s[k] || continue  # Il faut que l'observateur au sol puisse voir le satellite

        # Vecteur entre le satellite et le point au sol
        dx = sx - rx
        dy = sy - ry
        dz = sz - rz
        ρ = sqrt(dx*dx + dy*dy + dz*dz) # Distance entre le satellite et le point au sol
        invρ = 1.0 / ρ

        # Vecteur unitaire qui pointe vers le satellite
        ux = dx * invρ
        uy = dy * invρ
        uz = dz * invρ

        # Construction de la matrice G
        g11 += ux*ux; g12 += ux*uy; g13 += ux*uz; g14 += ux
        g22 += uy*uy; g23 += uy*uz; g24 += uy
        g33 += uz*uz; g34 += uz
        g44 += 1.0

        m += 1
    end

    m >= min_sats || return Inf  # PDOP non défini si il y a moins de 4 satellites au dessus du point au sol

    # On met toutes les composantes de G dans une matrice (symétrique)
    G = Symmetric([
        g11 g12 g13 g14;
        g12 g22 g23 g24;
        g13 g23 g33 g34;
        g14 g24 g34 g44
    ])

    try
        C = cholesky(G) # Calcul de la décomposition de Cholesky (G = produit d'une matrice triangulaire et de sa transposée)
        Q = C \ I(4) # Inversion
        return sqrt(Q[1,1] + Q[2,2] + Q[3,3]) # Le PDOP est la racine de la trace de la sous-matrice position (juste en 3D et pas 4D)
    catch
        return Inf # Infini si on ne peut pas calculer la décomposition de Cholesky car la matrice est mal conditionnée?
    end
end

"""
    pdop_point(sats, t, lat_deg, lon_deg, eps_deg; min_sats=4)

Calcule le Position Dilution of Precision (PDOP 3D) pour un point sol donné à un instant donné, en fonction de la géométrie de la constellation.

# Arguments
- sats      : Constellation de satellites.
- t         : Temps considéré (en secondes).
- lat_deg   : Latitude du point sol (en degrés).
- lon_deg   : Longitude du point sol (en degrés).
- eps_deg   : Angle d'élévation minimal (en degrés).

# Paramètres optionnels
- min_sats  : Nombre minimal de satellites visibles requis pour le calcul (défaut = 4).

# Valeur retournée
- PDOP 3D. Retourne `Inf` si le nombre de satellites visibles est insuffisant.
"""
function pdop_point(sats, t, lat_deg, lon_deg, eps_deg; min_sats=4)
    ϕ = deg2rad(lat_deg)
    λ = deg2rad(lon_deg)

    # Vecteur qui pointe du centre de la Terre vers le point au sol défini par sa latitude et sa longitude. (repère ECEF)
    gx = cos(ϕ) * cos(λ)
    gy = cos(ϕ) * sin(λ)
    gz = sin(ϕ)

    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats]
    ρs = map(norm, r_ecef)
    cosψmax_s = (Re ./ ρs) .* cos(deg2rad(eps_deg)) # Angle max pour l'observateur

    return pdop_calcul(r_ecef, ρs, cosψmax_s, gx, gy, gz; min_sats=min_sats)
end

"""
    pdop(sats, t, grid::GroundGrid, eps_deg; min_sats=4)

Calcule le PDOP ainsi que la couverture sur une grille de points sol pour une constellation donnée à un instant donné.

# Arguments
- sats      : Constellation de satellites.
- t         : Temps considéré (en secondes).
- grid      : Grille de points sol (`GroundGrid`).
- eps_deg   : Angle d'élévation minimal (en degrés).

# Paramètres optionnels
- min_sats  : Nombre minimal de satellites visibles requis (défaut = 4).

# Valeur retournée
- mean_pdop : PDOP moyen calculé sur les points où il est défini.
- cov       : Coverage (%) = % de points de la grille avec PDOP défini (>= min_sats visibles).
"""
function pdop(sats, t, grid::GroundGrid, eps_deg; min_sats=4)
    r_ecef = [ecef_from_eci(eci_pos(s, t), t) for s in sats]
    ρs = map(norm, r_ecef)
    cosψmax_s = (Re ./ ρs) .* cos(deg2rad(eps_deg))  # Angle max pour l'observateur

    # Préallocation des différents vecteurs pour le multi-threading
    nt = Threads.nthreads()
    sum_pdop = zeros(Float64, nt)
    cnt_ok = zeros(Int, nt)
    cnt_tot = zeros(Int, nt)

    nlon = length(grid.lons)

    Threads.@threads for i in eachindex(grid.lats)
        tid = Threads.threadid()
        cosϕ = grid.cosϕ[i] # On précalcule cosϕ pour ne pas le calculer plusieurs fois
        sinϕ = grid.sinϕ[i] # On précalcule sinϕ pour ne pas le calculer plusieurs fois

        local_sum = 0.0
        local_ok = 0
        local_tot = 0

        @inbounds for j in 1:nlon
            gx = cosϕ * grid.cosλ[j]
            gy = cosϕ * grid.sinλ[j]
            gz = sinϕ

            p = pdop_calcul(r_ecef, ρs, cosψmax_s, gx, gy, gz; min_sats=min_sats)
            local_tot += 1
            if isfinite(p)
                local_sum += p
                local_ok += 1
            end
        end

        sum_pdop[tid] += local_sum
        cnt_ok[tid] += local_ok
        cnt_tot[tid] += local_tot
    end

    S = sum(sum_pdop)
    ok = sum(cnt_ok)
    tot = sum(cnt_tot)

    mean_pdop = ok == 0 ? Inf : S / ok
    cov = tot == 0 ? 0.0 : 100.0 * ok / tot  # coverage = % points avec PDOP défini (>= min_sats visibles)

    return mean_pdop, cov
end

"""
    mean_pdop(sats, grid::GroundGrid, eps_deg; n=50, min_sats=4)

Calcule des statistiques moyennes de PDOP sur une période orbitale complète,
en moyennant les résultats obtenus à plusieurs instants.

Le coverage `cov` est défini comme le pourcentage de points de la grille
pour lesquels le PDOP est défini (>= min_sats satellites visibles).

# Arguments
- sats      : Constellation de satellites.
- grid      : Grille de points sol (`GroundGrid`).
- eps_deg   : Angle d'élévation minimal (en degrés).

# Paramètres optionnels
- n         : Nombre d'instants d'échantillonnage sur une période orbitale (défaut = 50).
- min_sats  : Nombre minimal de satellites visibles requis (défaut = 4).

# Valeur retournée
- mean_pdop : PDOP moyen spatio-temporel (sur les points où il est défini).
- cov       : Coverage moyen (%) sur la période orbitale.
"""
function mean_pdop(sats, grid::GroundGrid, eps_deg; n=50, min_sats=4)
    a = sats[1].a
    T = 2π * sqrt(a^3 / μ)
    ts = range(0, T; length=n)

    s_pdop = 0.0
    s_cov = 0.0

    for t in ts
        mp, cov = pdop(sats, t, grid, eps_deg; min_sats=min_sats)
        s_pdop += mp
        s_cov += cov
    end

    return s_pdop / length(ts), s_cov / length(ts)
end

"""
    eval_constellation_pdop(vec, F, i_deg, a, eps_deg, grid::GroundGrid; n=50, min_sats=4)

Évalue les performances PDOP d'une constellation définie par un vecteur de paramètres,
indépendamment du processus d'optimisation.

Le PDOP est moyenné sur une grille de points sol et sur une période orbitale complète.
Le coverage `cov` est défini comme le pourcentage de points où le PDOP est défini
(>= min_sats satellites visibles).

# Arguments
- vec       : Vecteur de paramètres définissant la constellation.
- F         : Paramètre de structure/phasage utilisé par `myconstellation`.
- i_deg     : Inclinaison orbitale (en degrés).
- a         : Demi-grand axe orbital (en mètres).
- eps_deg   : Angle d'élévation minimal (en degrés).
- grid      : Grille de points sol (`GroundGrid`).

# Paramètres optionnels
- n         : Nombre d'instants d'échantillonnage sur une période orbitale (défaut = 50).
- min_sats  : Nombre minimal de satellites visibles requis (défaut = 4).

# Valeur retournée
- mean_pdop : PDOP moyen spatio-temporel de la constellation.
- cov       : Coverage moyen (%) associé au PDOP (>= min_sats visibles).
- Ns        : Nombre total de satellites.
"""
function eval_constellation_pdop(vec, F, i_deg, a, eps_deg, grid::GroundGrid; n=50, min_sats=4)
    sats = myconstellation(vec, F, i_deg, a)
    mpdop, cov = mean_pdop(sats, grid, eps_deg; n=n, min_sats=min_sats)
    return mpdop, cov, length(sats)
end
