### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 048533b8-7fce-4f31-9e5f-e8ee2046192b
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 8c97cfeb-8ba7-49e8-a351-175edc2b5e93
using SatelliteToolbox, ReferenceFrameRotations, Geodesy, Graphs, Unitful, DataFrames, Orbits

# ╔═╡ 0aed51ba-d469-4888-bfd9-ffb164611dcf
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

# ╔═╡ 844571a7-840f-4e86-b9bc-7d20140bbdb5
const μ_earth = 3.986004418e14  # m³/s²

# ╔═╡ d1861088-71a5-4af3-ab72-230662ae226a
function rv_from_elements(a,e,i,Ω,ω,M,μ)
    E = M
    ν = 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
    r_norm = a*(1 - e*cos(E))
    r_perif = [r_norm*cos(ν), r_norm*sin(ν), 0.0]
    R3Ω = [cos(Ω) -sin(Ω) 0; sin(Ω) cos(Ω) 0; 0 0 1]
    R1i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)]
    R3ω = [cos(ω) -sin(ω) 0; sin(ω) cos(ω) 0; 0 0 1]
    Q = R3Ω*R1i*R3ω
    r = Q*r_perif
    return r
end

# ╔═╡ 33d65ab0-ba3e-11f0-0967-19042683053c
begin
	c=299_792_458.0
	Re=6371e3
	alts=collect(500e3:100e3:1500e3)
	Nsats=collect(48:48:288)
	latgrid=-60.0:2.0:60.0; longrid=-180.0:5.0:180.0
	elevmin=10.0
	function make_walker(N,P,F; i=60.0, h=550e3)
	    a=Re+h
	    sats=[]
	    for p in 0:P-1, s in 0:(N÷P)-1
	        Ω=rad2deg(360p/P); M=rad2deg(360s/(N÷P) + 360F*p/N)
	        push!(sats, (a=a,e=0.0,i=i,Ω=Ω,ω=0.0,M=M))
	    end
	    sats
	end
	function latency_stats(h,N)
	    P=round(Int, sqrt(N)); sats=make_walker(N,P,1; i=60.0, h=h)
	    tspan=0:60:3600
	    lats=[]
	    for t in tspan
	        xyz = [first(rv_from_elements(s.a, s.e, deg2rad(s.i), deg2rad(s.Ω), deg2rad(s.ω), deg2rad(s.M), μ_earth)) for s in sats]
	        R = rot_ECI_to_ECEF(gmst(date_to_jd(DateTime(2000,1,1) + Millisecond(t*1000))))
	        xyz_ecef=[R*x for x in xyz]
	        for ϕ in latgrid, λ in longrid
	            p=llh2ecef(LatLon(ϕ,λ), h0=0.0)
	            vis=false; dmin=Inf
	            for x in xyz_ecef
	                v=x.-p; up=p./norm(p); elev=asin( dot(v./norm(v), up) )
	                if elev>deg2rad(elevmin)
	                    d=norm(v); dmin=min(dmin,2d)
	                    vis=true
	                end
	            end
	            if vis; push!(lats,dmin/c); end
	        end
	    end
	    (;cover=length(lats)/(length(tspan)*length(latgrid)*length(longrid)), Lmean=mean(lats), Lp95=quantile(lats,0.95))
	end
	res=[(h=h,N=N, latency_stats(h,N)...) for h in alts for N in Nsats]
	df=DataFrame(res)
	
end

# ╔═╡ Cell order:
# ╟─0aed51ba-d469-4888-bfd9-ffb164611dcf
# ╟─048533b8-7fce-4f31-9e5f-e8ee2046192b
# ╠═8c97cfeb-8ba7-49e8-a351-175edc2b5e93
# ╠═844571a7-840f-4e86-b9bc-7d20140bbdb5
# ╟─d1861088-71a5-4af3-ab72-230662ae226a
# ╠═33d65ab0-ba3e-11f0-0967-19042683053c
