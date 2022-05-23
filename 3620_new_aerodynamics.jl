### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ dce5c906-6425-471a-9b89-114feaf02a38
begin
	using AeroMDAO      # Main package
	using Plots# Plotting library
	using PlutoUI
	using JuMP # Wrapper for Optimization solver
	using Ipopt # Optimization model
	using Optim# Optimization Package
	using LinearAlgebra
	gr(dpi = 300)       # Plotting backend
	using LaTeXStrings  # For LaTeX printing in plots
	theme(:default)
	PlutoUI.TableOfContents(title = "AeroMDAO", depth = 4)
end

# ╔═╡ 8a679de5-f08f-44c1-bd14-7bac0e58da90
md"## Javenlin XS
Aerodynamic Analysis -DAEYEONG KWUN
"

# ╔═╡ a550479a-df78-4b09-b213-d0c237d052b7
md"Load Dependencies.."

# ╔═╡ 8b6d0ff5-c7c2-465f-8de3-779441fc13cc
# Tasks
# 1. Optimizer Variables
# Changes
# 1. More detail on Wing Airfoil Optimizatino mid chord - decreased drag
# 2. Htail design
# 3. Profile drag
# 4. Nearfield for Spanwise Analysis, farfield for other things

# ╔═╡ 6e4db00a-e601-4519-aad2-0f456be3ba6a
md"### Visualization"

# ╔═╡ 64b849b4-e761-4e50-95c3-996aea337e27
begin 
	φ = @bind φ Slider(0:1e-2:90, default = 35)
	ψ = @bind ψ Slider(0:1e-2:90, default = 47)
	strlines = @bind strmlines CheckBox()
	alpha_n = @bind alpha_n NumberField(0:25; default = 0)
	beta_n = @bind beta_n Slider(-45:1e-2:45, default = 0)
	airspeed = @bind airspeed Slider(0:1e-2:229.805, default = 1)
	# z_s = @bind z_limit Slider(0:1e-2:10, default = 8.)
	# *z*-scale: $(z_s) # Put this in md if necessary
	md"""
	Horizontal: $(φ)
	Vertical: $(ψ)
	
	AoA: $(alpha_n)
	Beta: $(beta_n)
	Air Speed: $(airspeed)
	
	Streamlines: $strlines
	"""
end

# ╔═╡ 06d14669-17dd-42e6-9287-02cc518a796c
md"Max Air Speed is set as the cruising speed."

# ╔═╡ 72fd45c4-ce91-4f49-9a4f-b6e695fc405b
md"### Codes"

# ╔═╡ c6a480e4-a43d-4a1e-a0a3-01269a414272
md"Load airfoil files"

# ╔═╡ 9b5900f3-de26-48e0-887b-e40d31c9ba77
begin
	# Airfoil coordinates file path
	foilpath_18 = string(@__DIR__, "/naca643418.dat")
	foilpath_15 = string(@__DIR__, "/naca642415.dat")
	
	# Read coordinates file
	airfoil_18 = read_foil(foilpath_18;
                    name   = "NACA 643418" # To overwrite name in header
                   )
	airfoil_15 = read_foil(foilpath_15;
                    name   = "NACA 642415" # To overwrite name in header
                   )
end;

# ╔═╡ 59da49fb-d572-4a29-94f7-1e9532cf2211
md"Define Airfoils"

# ╔═╡ 2ba9a5d6-9b30-4b62-a8ad-e0271bf519bb
begin
	# Define vector of airfoils
	airfoils  = [ airfoil_18, airfoil_18, airfoil_15 ]
	# Convert Airfoil Coordinates into Camber Thickness
	xcamthick_18 = camber_thickness(airfoil_18, 60)
	xcamthick_15 = camber_thickness(airfoil_15, 60)
	xs_18, cambers_18, thiccs_18 = xcamthick_18[:,1], xcamthick_18[:,2], xcamthick_18[:,3];
	xs_15, cambers_15, thiccs_15 = xcamthick_15[:,1], xcamthick_15[:,2], xcamthick_15[:,3];
	# Maximum thickness-to-chord ratio 
	# (t/c)_max
	tc_max_18 = maximum(thiccs_18)
	tc_max_15 = maximum(thiccs_15)
	# Location of (t/c)_max
	# (x/c)_max
	max_tc_arg_18 = argmax(thiccs_18)
	max_tc_arg_15 = argmax(thiccs_15)
	x_tc_max_18 = xs_18[max_tc_arg_18]
	x_tc_max_15 = xs_15[max_tc_arg_15]
end;

# ╔═╡ 81bb09b6-ebbb-4e90-a119-a62ba10d7a68
begin 
	n = 10
	plot(xlabel = "x", ylabel = "y", aspect_ratio = 1)
	plot!(airfoil_18.x, airfoil_18.y, label = "$(airfoil_18.name) Coordinates")
	plot!(xs_18, thiccs_18, ls = :dash, label = "Thickness")
	plot!(xs_18, cambers_18, ls = :dash, label = "Camber")
	plot!(airfoil_18.x, zero(airfoil_18.x), ls = :dot, label = :none)
	plot!(fill(x_tc_max_18, n), LinRange(0, tc_max_18, n), ls = :dot, label = :none)
	scatter!([x_tc_max_18], [tc_max_18], label = "Max (t/c) at (x/c) = $(round(x_tc_max_18; digits = 4))")
end

# ╔═╡ 21394943-8cf5-4ee0-9027-7800aff65649
begin 
	plot(xlabel = "x", ylabel = "y", aspect_ratio = 1)
	plot!(airfoil_15.x, airfoil_15.y, label = "$(airfoil_15.name) Coordinates")
	plot!(xs_15, thiccs_15, ls = :dash, label = "Thickness")
	plot!(xs_15, cambers_15, ls = :dash, label = "Camber")
	plot!(airfoil_15.x, zero(airfoil_15.x), ls = :dot, label = :none)
	plot!(fill(x_tc_max_15, n), LinRange(0, tc_max_15, n), ls = :dot, label = :none)
	scatter!([x_tc_max_15], [tc_max_15], label = "Max (t/c) at (x/c) = $(round(x_tc_max_15; digits = 4))")
end

# ╔═╡ 3bdb49c0-de68-4b83-ba2f-02688c91a899
md"Define Wing"

# ╔═╡ 176509dc-b4f0-11ec-08a6-317c4af234c4
begin
	wing = Wing( # 2 spanwise sections ⇒ 3 chord specifications
		 foils = airfoils,
	     chords = [2.8, 1.4, 0.92], #1.73
	     twists   = [ 3.17, 3.17, -1.37 ], # Root to tip (degrees)
	     spans   = fill(6.9, 2),
	     dihedrals = fill(4., 2),
	     sweeps = fill(30, 2)
	     )
	# Horizontal tail
	htail = Wing(foils     = fill(naca4(0,0,1,2), 2),
             chords    = [1., 0.6],
             twists    = [0.0, 0.0],
             spans     = [4],
             dihedrals = [0.],
             sweeps    = [20],
             w_sweep   = 0.5,
             position  = [5.8, 0, 0.5],
             angle     = -2.,
             axis      = [0., 1., 0.])
	# Vertical tail
	vtail = HalfWing(foils     = fill(naca4(0,0,0,9), 2),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.],
                 dihedrals = [0.],
                 sweeps    = [0],
                 w_sweep   = 0.65,
                 position  = [5.8, 0, 0.5],
                 angle     = 90.,
                 axis      = [1., 0., 0.])
	# Generate the coordinates of the wing's outline
	wing_outline = plot_planform(wing) 
	htail_outline = plot_planform(htail)
	vtail_outline = plot_planform(vtail)
end;

# ╔═╡ 542efc07-818e-42ca-a667-3d13b7d40756
md"Meshing"

# ╔═╡ 2a1e745d-dce0-464a-ae23-cf94beb1767b
begin
	wing_mesh = WingMesh(wing, [20, 5], 10)
	htail_mesh = WingMesh(htail, [6], 4)
	vtail_mesh = WingMesh(vtail, [4], 3)
	# Compute camber panel distribution
	wing_cam_panels = camber_panels(wing_mesh)
	htail_cam_panls = camber_panels(htail_mesh)
	vtail_cam_panls = camber_panels(vtail_mesh)
	# Generate plotting points
	plt_wing_pans = plot_panels(wing_cam_panels)
	plt_htail_pans = plot_panels(htail_cam_panls)
	plt_vtail_pans = plot_panels(vtail_cam_panls)
end;

# ╔═╡ 18e0076a-c70e-452b-a1f6-d0dad0303a39
md"Assemble horseshoes on each surface"

# ╔═╡ 65a3845f-c6c5-4d86-8be4-d83c34800d1f
begin
	# Assemble horseshoes on each surface into a component vector
	aircraft = ComponentVector(
	                           wing  = make_horseshoes(wing_mesh),
	                           htail = make_horseshoes(htail_mesh),
	                           vtail = make_horseshoes(vtail_mesh)
	)
end;

# ╔═╡ 5644b1b4-12de-4628-8a07-a03a27244025
md"Visualization of streamines around the wing"

# ╔═╡ 24a8d3fc-61fe-42f8-8173-5e8557ae28ad
md"Freestream Condition"

# ╔═╡ fd92271f-cdfc-49c6-8487-ab7d57b04d28
begin
	fs  = Freestream(
	                 alpha = alpha_n, # degrees
	                 beta  = beta_n, # degrees
	                 omega = [0., 0., 0.]
	                )
	refs = References(
                  speed    = airspeed,
                  area     = projected_area(wing),
                  span     = span(wing),
                  chord    = mean_aerodynamic_chord(wing),
                  density  = 1.225,
				  viscosity = 1.5e-5,
                  location = mean_aerodynamic_center(wing)
                 )
	system = solve_case(
            aircraft, fs, refs;
            print            = true, # Prints the results for only the aircraft
            print_components = true, # Prints the results for all components
	)
	nf = nearfield(system)
	ff = farfield(system)
	nf_coeffs = nearfield_coefficients(system)
	ff_coeffs = farfield_coefficients(system)
end;

# ╔═╡ 20b2fccf-6c9e-477a-b795-2dacd7d1d497
begin
	CFs, CMs = surface_coefficients(system)

	# Compute spanwise loads
	span_loads  = spanwise_loading(wing_mesh, CFs.wing, projected_area(wing))
	CL_loads    = vec(sum(system.circulations.wing, dims = 1)) / (0.5 * refs.speed * refs.chord);
	
	# Plot spanwise loadings
	plot_CD = plot(span_loads[:,1], span_loads[:,2], label = :none, ylabel = L"C_{D_i}")
	plot_CY = plot(span_loads[:,1], span_loads[:,3], label = :none, ylabel = L"C_Y")
	plot_CL = begin
	            plot(span_loads[:,1], span_loads[:,4], label = :none, xlabel = L"y", ylabel = L"C_L")
	            plot!(span_loads[:,1], CL_loads, label = "Normalized", xlabel = L"y")
	          end
	plot(plot_CD, plot_CY, plot_CL, size = (800, 700), layout = (3,1))

end

# ╔═╡ 136e6b80-0f1b-43d5-b35e-2a82b8777cf8
begin
	# Streamlines generation
	span_points = 20
	if fs.beta == 0
		surf = wing.right
	else
		surf = wing
	end
	init        = chop_leading_edge(surf, span_points) 
	dx, dy, dz  = 0, 0, 1e-3
	seed_2      = [ init .+ Ref([dx, dy, dz])  ;
					init .+ Ref([dx, dy,-dz])  ]
	seed = seed_2;
	streams = ifelse(strmlines, plot_streamlines(system, seed, 2, 100), nothing);
end;

# ╔═╡ 0fc8cd9d-d268-49a3-bef5-81bc205c9dac
begin
	plt = plot(
	           wing_outline[:,1], wing_outline[:,2], wing_outline[:,3],
	           label = "Wing",
	           xaxis = "x", yaxis = "y", zaxis = "z",
	           aspect_ratio = 0.6,
	           camera = (φ, ψ),
	           zlim = (-0.1, span(wing) / 2),
	)
	plot!(plt,
    htail_outline[:,1], htail_outline[:,2], htail_outline[:,3],
    label = "Horizontal Tail"
   )
	plot!(plt,
   vtail_outline[:,1], vtail_outline[:,2], vtail_outline[:,3],
   label = "Vertical Tail"
  )
	# Mesh
	[ plot!(plt, panel, label = "", color = :lightblue)
	    for panel in plt_wing_pans ]
	
	[ plot!(plt, panel, label = "", color = :orange)
	    for panel in plt_htail_pans ]
	[ plot!(plt, panel, label = "", color = :lightgreen)
	    for panel in plt_vtail_pans ]
	# Streamline
	if strmlines
		[ plot!(Tuple.(stream), color = :green, label = :none) 
			for stream in eachcol(streams) ]
	end
	plot!(plt)
end

# ╔═╡ 5ab3b285-47db-49a5-854f-c45d8023f8c4
md"Drag Polar with varying AoA"

# ╔═╡ a2e12b5b-91d1-496b-a327-c3bbf2ab897f
md"Skin Friction Drag Estimation"

# ╔═╡ 68ae084d-a483-491b-8f84-fa8eb445fa6a
begin
	# Fully Laminar
	cf_lam(Re_ℓ) = 1.328 / √Re_ℓ
	# Fully Turbulent with compressible correctinos
	cf_turb(Re_ℓ, M) = 0.455 / log10(Re_ℓ)^2.58 / (1 + 0.144M^2)^0.65
	#Fitting from laminar to turbulent with transition
	cf_trans(Re_ℓ, Re_xtr, M) = max(cf_lam(Re_ℓ), cf_turb(Re_ℓ, M) - (Re_xtr / 320 - 39) / Re_ℓ)
	# Empirical Correction for curved shapes
	form_factor(x_c, t_c, Λ_m, M) = (1 + 0.6t_c / x_c + 100t_c^4) * (1.34M^0.18 * cos(Λ_m)^0.28)
	# I don't know what this is
	xcs, tcs, lams = max_thickness_to_chord_ratio_sweeps(wing.right, 60)
	xcs_h, tcs_h, lams_h = max_thickness_to_chord_ratio_sweeps(htail.right, 60)
	xcs_v, tcs_v, lams_v = max_thickness_to_chord_ratio_sweeps(vtail, 60)
	Kfs_w = form_factor.(xcs, tcs, lams, (airspeed/335))
	Kfs_h = form_factor.(xcs_h, tcs_h, lams_h, (airspeed/335))
	Kfs_v = form_factor.(xcs_v, tcs_v, lams_v, (airspeed/335))
	# Average form factor over the entire planform
	Kf_w = sum(Kfs_w) / length(Kfs_w)
	Kf_h = sum(Kfs_h) / length(Kfs_h)
	Kf_v = sum(Kfs_v) / length(Kfs_v)
	# Calculate the wetted area ratio
	S_wet_Sw = wetted_area_ratio(wing_mesh)
	S_wet_Sw_h = wetted_area_ratio(htail_mesh)
	S_wet_Sw_v = wetted_area_ratio(vtail_mesh)
	# Density (kg/m³)
	ρ    = 1.225    
	# Dynamic viscosity (kg m⁻¹ s⁻¹)
	μ    = 1e-6 	
	# Transition location as ratio of chord length
	x_tr_w = 0.7
	x_tr_h = 0.3
	# Chord length (m)
	c_w  = mean_aerodynamic_chord(wing)
	c_h  = mean_aerodynamic_chord(htail)
	c_v  = mean_aerodynamic_chord(vtail)
	# Reynolds Number Function
	reynolds_number(ρ, V, c, μ) = ρ * V * c / μ
	# Reynolds number Re
	Re_w = reynolds_number(ρ, airspeed, c_w, μ)
	Re_h = reynolds_number(ρ, airspeed, c_h, μ)
	Re_v = reynolds_number(ρ, airspeed, c_v, μ)
	# Reynolds number at transition
	Re_xtr = reynolds_number(ρ, airspeed, x_tr_w * c_w, μ)
	Re_xtr_h = reynolds_number(ρ, airspeed, x_tr_h * c_h, μ)
	Re_xtr_v = reynolds_number(ρ, airspeed, x_tr_h * c_v, μ)
	cf = cf_trans(Re_w, Re_xtr, (airspeed/335))
	cf_h = cf_trans(Re_h, Re_xtr_h, (airspeed/335))
	cf_v = cf_trans(Re_v, Re_xtr_v, (airspeed/335))
	# Profile Drag Coefficients
	CDp_w = 2 * S_wet_Sw * cf * Kf_w # 2 for the entire wing
	CDp_h = 2 * S_wet_Sw_h * cf_h * Kf_h 
	CDp_v = S_wet_Sw_v * cf_v * Kf_v
	# Viscous Analysis (Empirical)
	CDv_wing = profile_drag_coefficient(wing_mesh, x_tr_w, refs)
	CDv_htail = profile_drag_coefficient(htail, x_tr_h, refs)
	CDv_vtail = profile_drag_coefficient(vtail, x_tr_h, refs)
	CD_ff_wing = CDv_wing + ff_coeffs.wing[1]
	CD_ff_htail = CDv_htail + ff_coeffs.htail[1]
	CD_ff_vtail = CDv_vtail + ff_coeffs.vtail[1]
	# Total Profile Drag Coefficient
	CDv = CDv_wing + CDv_htail + CDv_vtail
	# Total Drag Coefficient
	CD_nf = CDv + nf.CDi
	CD_ff = CDv + ff.CDi_ff
	est_l_d = nf.CL / CD_nf
end;

# ╔═╡ beaf63ea-08cd-4ad9-819e-2fb18de27bc0
begin
	# Define function to compute system varying with angle of attack.
	vary_alpha(aircraft, α, refs) = solve_case(aircraft, Freestream(alpha = α), refs)
	
	# Run loop
	αs      = -5:0.5:25
	systems = [ vary_alpha(aircraft, α, refs) for α in αs ]
	# Cleaner: map(α -> vary_alpha(...), αs)
	
	# Get coefficients
	coeffs =  farfield.(systems) # Should i use nearfield or farfield?
	CDis   = [ c[1] for c in coeffs ]
	CDs   = [ c[1] + CDv for c in coeffs ] # CDi + CDp
	CLs    = [ c[3] for c in coeffs ];
end;

# ╔═╡ 30456a2d-904c-459a-a48b-92b5e2a5d030
begin
	# Concatenate results into one array
	data = permutedims(reduce(hcat, [α; c...] for (α, c) in zip(αs, coeffs)))
	
	# Plot
	plot(data[:,1], round.(data[:,2:end], digits = 4),
	     layout = (3,2),
	     xlabel = L"\alpha",
	     ylabel = [L"C_{D_i}" L"C_Y" L"C_L" L"C_\ell" L"C_m" L"C_n" ],
	     labels = "",
	     size   = (800, 600)
	    )
end

# ╔═╡ d40c579e-7bc2-4d3d-9473-80a966676d74
md"V-n Diagram"

# ╔═╡ 15e7100e-8b71-40b8-9535-dd03f6c57ede
begin
	# Number of points for generating lines
	num = 150
	# Wing Area
	wing_area = projected_area(wing)
	# Mean Aerodynamic chord
	mac = mean_aerodynamic_chord(wing)
	# Max Take Off Weight
	MTOW = 78100 # in Newton
	# Max Wing Loading
	max_wing_loading = MTOW/wing_area
	# Cruise Speed
	v_cruise = 0.77 * 316
	# Lift curve slope
	cl_alpha = 5.88
	# Max. Lift Coefficient
	cl_max = 2.08
	# Min. Lift Coefficient
	cl_min = -1.0 # -0.054735
	# % of the max wing loading
	factor 	= 0.6   		
	# Wing loading
	WbyS 	= factor * max_wing_loading
	# Altitude units ("ft" or "m")
	unit    = "ft" 
	# Altitude for analysis (not cruise here)
	alt 	= 35000
	g = 9.81 # m/s²
	
	function quadratic_roots(a, b, c)
    d = sqrt(b^2 - 4a*c)
    (-b .+ (d, -d)) ./ 2a
	end

	# Equivalent airspeed conversion
	EAS(V, ρ_alt, ρ_SL) = V * sqrt(ρ_alt / ρ_SL)
	# Atmospheric model from NASA
	function density(alt; units = "m")
		alt = ifelse(units == "m", 3.28084 * alt, alt)
		if alt <= 36152
			T = 59 - 3.56e-3 * alt
			p = 2116 * ((T + 459.7) / 518.6)^5.256
		elseif 36152 < alt < 82345
			T = -70
			p = 473.1 * exp(1.73 - 4.8e-5 * alt)
		else
			T = -205 + 0.00164 * alt
			p = 51.97 * ((T + 459.7) / 389.98)^(-11.388)
		end	
	ρ = p / (1718 * (T + 459.7)) * 515.379 # kg/m³
	end
	rho_SL = density(6000., units = "m") # kg/m³
	max_lift = 0.5 * rho_SL * v_cruise^2 * wing_area * cl_max 
	
	# Maximum Permissible Speeds
	# Density at altitude
	rho_alt = density(alt, units = unit)
	# Design Speed
	V_C = EAS(v_cruise, rho_alt, rho_SL)
	# Max Operating Speed
	V_MO = V_C * 1.06
	# Design Diving Speed
	V_D = V_MO * 1.07

	# Load Factor Limits
	# According to FAR-25 specifications for a jet aircraft
	n_max = 2.44  # Maximum load factor
	n_min = -1.0 # Minimum load factor
	ns = range(n_min - 0.4, n_max + 0.4, length = num)
	V_Ds = range(0, V_D, length = num) # knots

	# Maneuver Point
	load_factor(ρ, V, CL, WbyS) = ρ * V^2 * CL / 2WbyS
	speed_from_load(n, ρ, WbyS, CL) = sqrt(2WbyS / ρ * abs(n / CL))
	# Corner Speed Va
	V_A = speed_from_load(n_max, rho_SL, WbyS, cl_max) # m/s
		V_As  = range(0, V_A, length = num)
	n_pos = load_factor.(rho_SL, V_As, cl_max, WbyS)
	
	V_n_min = speed_from_load(n_min, rho_SL, WbyS, cl_min)
	V_n_max = speed_from_load(n_max, rho_SL, WbyS, cl_max)
	
	V_nminCs = range(V_n_min, V_C, length = num)
	V_nmaxDs = range(V_n_max, V_D, length = num)
	
	V_Amins = range(0, V_n_min, length = num)
	n_neg   = load_factor.(rho_SL, V_Amins, cl_min, WbyS)
	n_C_min = range(n_neg[end], n_min, length = num)
	
	V_CDs = range(V_C, V_D, length = num)
	n_CDs = range(n_min, 0., length = num)

	# Stall Speed
	n_stall = 1.0
	V_S = speed_from_load(n_stall, rho_SL, WbyS, cl_max)
	
	n_pos_stall = load_factor(rho_SL, V_S, cl_max, WbyS)
	n_neg_stall = load_factor(rho_SL, V_S, cl_min, WbyS)
	
	ns_stall = range(n_neg_stall, n_pos_stall, length = num)
	
	V_SAs = range(V_S, V_A, length = num)
	V_Snmins = range(V_S, V_n_min, length = num)
	
	ns_pos_stall = load_factor.(rho_SL, V_SAs, cl_max, WbyS)
	ns_neg_stall = load_factor.(rho_SL, V_Snmins, cl_min, WbyS)
	
	eps = 0.05 # Offset factor for annotations

	# Gust Speeds
	n_gust(Kg, CL_alpha, Ue, V, WbS; positive = true) = 1 + ifelse(positive, 1, -1) * Kg * CL_alpha * Ue * V / 2WbS
	mu(WbyS, ρ, c, CL_alpha, g) = 2WbyS / (ρ * c * CL_alpha * g)
	mu_case = mu(WbyS, rho_alt, mac, cl_alpha, g) # What is chord?
	Kg = (0.88 * mu_case) / (5.3 + mu_case)
	n_gust_VC = n_gust(Kg, cl_alpha, 50, 20, WbyS) 
	
	linear_interp(x, y0, y1, x0, x1) = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
	types = [:knots, :ftps, :mps]
	
	function gust_values(alt; units = "m")
		alt = ifelse(units == "m", 3.28084 * alt, alt)
		# Extrema 
		alt_min, alt_max = (20000, 50000)
		Ue_max, Ue_min   = (66, 50, 25) .* 0.3048, (38, 25, 12.5) .* 0.3048
		
		# Conditions
		if alt <= alt_min
			Ue_max
		elseif alt >= alt_max
			Ue_min
		else
			Ue_B, Ue_C, Ue_D = linear_interp.(alt, Ue_min, Ue_max, alt_max, alt_min)
		end
	end
	function gust_load_factor(Kg, CL_alpha, Ue, CL_max, rho, WbyS)
		# Define coefficients for quadratic equation
		a = rho * CL_max / 2WbyS
		b = - Kg * CL_alpha * Ue / (2WbyS)
		c = -1
			
		# Evaluate roots
		VB, VB_wrong = quadratic_roots(a, b, c)
	
		# Return correct root
		return VB
	end
	Ue_B, Ue_C, Ue_D = gust_values(alt; units = "ft") # m/s
	V_B = gust_load_factor(Kg, cl_alpha, Ue_B, cl_max, rho_SL, WbyS) # m/s

	V_Bs = range(0, V_B, length = num)
	ngust_B_pos = n_gust.(Kg, cl_alpha, Ue_B, V_Bs, WbyS)
	ngust_B_neg = n_gust.(Kg, cl_alpha, Ue_B, V_Bs, WbyS, positive = false)
	
	V_Cs = range(0, V_C, length = num)
	ngust_C_pos = n_gust.(Kg, cl_alpha, Ue_C, V_Cs, WbyS)
	ngust_C_neg = n_gust.(Kg, cl_alpha, Ue_C, V_Cs, WbyS, positive = false)
	
	ngust_D_pos = n_gust.(Kg, cl_alpha, Ue_D, V_Ds, WbyS)
	ngust_D_neg = n_gust.(Kg, cl_alpha, Ue_D, V_Ds, WbyS, positive = false)
	
	ngust_BCs_pos = range(ngust_B_pos[end], ngust_C_pos[end], length = num)
	ngust_BCs_neg = range(ngust_B_neg[end], ngust_C_neg[end], length = num)
	ngust_CDs_neg = range(ngust_C_neg[end], ngust_D_neg[end], length = num)
	ngust_CDs_pos = range(ngust_C_pos[end], ngust_D_pos[end], length = num)
	
	V_ABs = range(V_A, V_B, length = num)
	V_BCs = range(V_B, V_C, length = num)
	
	ngust_ABs_pos = load_factor.(rho_SL, V_ABs, cl_max, WbyS)
	ngust_VDs = range(ngust_D_pos[end], ngust_D_neg[end], length = num)
	
	# Intersection candidates
	V_nmin_gusts = range(V_Snmins[end], V_B, length = num)
	ngust_nminBs = range(n_neg[end], ngust_B_neg[end], length = num)
	ngust_BCs_neg = range(ngust_nminBs[end], ngust_C_neg[end], length = num)
	ngust_nminB2s = load_factor.(rho_SL, V_nmin_gusts, cl_min, WbyS)
	ngust_BCs_neg2 = range(ngust_nminB2s[end], ngust_C_neg[end], length = num)
	V_B_nmin = speed_from_load(n_min, rho_SL, WbyS, cl_min)
	n_neg_VB = load_factor(rho_SL, V_B, cl_min, WbyS)
end;

# ╔═╡ 9e0099cc-f587-4bc5-b6f5-f60912069739
begin	
	gust = plot(xlabel = "Equivalent Airspeed, V_EAS",
		 ylabel = "Load Factor, n",
		 legend = :bottomleft,
		 title = "V-n Diagram - Gust Envelope")
	
	# Vertical lines
	plot!(fill(V_C, num), ns,
		  line = :dash, color=:gray, linewidth = 1,
		  label = :none
		 ) # Design speed
	plot!(fill(V_MO, num), ns, 
		  line = :dash, color=:gray, linewidth = 1,
		  label = :none
		 ) # Maximum operating speed
	plot!(fill(V_D, num), ns, 
		  line = :dash, color=:gray, linewidth = 1,
		  label = :none
		 ) # Design driving speed
	plot!(fill(V_S, num), ns,
		  line = :dash, color=:gray, linewidth = 1,
		  label = :none
	     ) # Stall speed
	plot!(fill(V_A, num), ns,
			  line = :dash, color=:gray, linewidth = 1,
			  label = :none
		 ) # Corner speed
	plot!(fill(V_B, num), ns,
			  line = :dash, color=:gray, linewidth = 1,
			  label = :none
		 ) # Minimum gust speed
	
	# Horizontal lines
	plot!(V_Ds, fill(n_max, num), 
		  line = :dash, color=:gray, linewidth = 1,
		  label = :none
	     ) # Maximum load factor
	plot!(V_Ds, fill(n_min, num),
		  line = :dash, color=:gray, linewidth = 1,
		  label = :none
	     ) # Minimum load factor
	plot!(V_Ds, zeros(num),
		  line = :dash, color=:gray, linewidth = 1,
		  label = :none
	     ) # Zero load factor
	plot!(V_Ds, ones(num),
		  line = :dash, color=:gray, linewidth = 1,
		  label = :none
		 ) # Unity load factor
	
	# Maneuverable envelope
	plot!(V_SAs, ns_pos_stall,
		  label= :none, color = :blue
		 ) # Positive load factor
	plot!(V_Snmins, ns_neg_stall,
		  label= :none, color = :blue
		 ) # Negative load factor
	
	plot!(V_nminCs, n_C_min,
		  label = :none, color = :blue
		 ) # Minimum load factor
	plot!(V_nmaxDs, fill(n_max, num),
		  label = :none, color = :blue
	     ) # Maximum load factor
	
	plot!(fill(V_D, num), range(0, n_max, length = num),
		  label = :none, color = :blue
		 ) # Maximum operating pseed
	plot!(V_CDs, n_CDs,
		  label = :none, color = :blue
		 ) # Linearly increasing from VC to VD
	
	plot!(fill(V_S, num), ns_stall,
		  label = :none, color = :blue
		 ) # Stall speed
	
	# Gust envelope
	plot!(V_Bs, ngust_B_pos,
		  color = :red, line = :dot, linewidth = 1.5,
		  label = "n_gust, B") # Positive gust, B
	plot!(V_Bs, ngust_B_neg,
		  color = :red, line = :dot, linewidth = 1.5,
		  label = :none) # Negative gust, B 
	
	plot!(V_Cs, ngust_C_pos,
		  color = :green, line = :dot, linewidth = 1.5,
		  label = "n_gust, C") # Positive gust, C
	plot!(V_Cs, ngust_C_neg,
		  color = :green, line = :dot, linewidth = 1.5,
		  label = :none) # Negative gust, C
	
	plot!(V_Ds, ngust_D_pos,
		  color = :blue, line = :dot, linewidth = 1.5,
		  label = "n_gust, D") # Positive gust, D
	plot!(V_Ds, ngust_D_neg,
		  color = :blue, line = :dot, linewidth = 1.5,
		  label = :none) # Negative gust, D
	
	# Reds 
	plot!(V_BCs, ngust_BCs_pos,
		  color = :red,
		  label = :none
		 ) # Positive load factors from VB to VC

	plot!(V_CDs, ngust_CDs_pos,
		  color = :red,
		  label = :none
		 ) # Positive load factors from VC to VD
	plot!(V_CDs, ngust_CDs_neg,
		  color = :red,
		  label = :none
		 ) # Negative load factors from VC to VD

	plot!(fill(V_D, num), ngust_VDs,
		  color = :red,
		  label = :none
		 ) # Load factors on VD
	
	if V_B > V_A
		plot!(V_ABs, ngust_ABs_pos,
		  color = :red,
		  label = :none
		 ) # Positive load factors from VA to VB
	end
	
 	if V_B < V_B_nmin
		plot!(V_BCs, ngust_BCs_neg,
		  color = :red,
		  label = :none
		 ) # Negative load factors from VB to VC
	else 
		if n_neg_VB > ngust_B_neg[end] # If the minimum n is the VB gust load
			plot!(V_nmin_gusts, ngust_nminB2s,
			  color = :red,
			  label = :none
			 ) # Negative load factors from n_neg to VB
			plot!(V_BCs, ngust_BCs_neg2,
			  color = :red,
			  label = :none
			 ) # Negative load factors from VB to VC
		else # If the minimum n is the minimum negative load factor
			plot!(V_BCs, ngust_BCs_neg,
			  color = :red,
			  label = :none
			 ) # Negative load factors from VB to VC
		end
		
	end

	# Vertices
	vertices_2 = [ V_S n_stall 			"VS" ;
				   V_A n_max 			"VA" ;
				   V_B ngust_B_pos[end] "VB" ; 
		 		   V_C ngust_C_pos[end] "VC" ; 
				   V_D ngust_D_pos[end]	"VD" ]
	scatter!(vertices_2[:,1] , vertices_2[:,2], label = :none)
	annotate!(vertices_2[:,1] .* (1 - eps), vertices_2[:,2] .* (1 + eps), vertices_2[:,3])
		
	plot!()
end

# ╔═╡ a70b8c36-c715-47fc-bf5d-5b95e9fc2692
md"Simple Optimizer (Honestly Im not sure if it works)"

# ╔═╡ f2778e6d-1c56-449b-8236-5c8a464bfb01
begin
	model = Model(Ipopt.Optimizer)

	# Free variables
	@variable(model, V, start = 100.)   # Cruise speed
	@variable(model, W, start = 10000.) # Total weight
	@variable(model, A, start = 6.) 	# Wing aspect ratio
	@variable(model, S, start = 12.)    # Wing area
	@variable(model, Ww, start = 2000.) # Wing weight

	# Dependent variables
	@variable(model, Re, start = 1e6)   # Reynolds number
	@variable(model, Cf, start = 0.001) # Skin-friction coefficient
	@variable(model, CL, start = 1.0)   # Cruise lift coefficient
	@variable(model, CD, start = 0.01)  # Cruise drag coefficient

	# Constraints
	# Form factor
	k       = 1.2 
	# Oswald span efficiency factor
	e       = 0.96   
	# Fuselage drag area
	CDA₀    = 0.0306  
	# Wetted area ratio
	S_wet_S = 2.05    
	# Aircraft weight excluding wing (N)
	W0      = 4940   
	# Wing weight coefficient 1
	Ww_c1   = 8.71e-5  
	# Wing weight coefficient 2
	Ww_c2   = 45.42    
	# Ultimate load factor
	N       = 2.5      
	# Airfoil thickness-to-chord ratio
	τ       = 0.12    
	# Takeoff airspeed (m/s)
	v_to = 22 
	# Reynolds number
	@NLconstraint(model, Re == ρ * V * (S / A)^0.5 / μ)
	# Skin-friction coefficient
	@NLconstraint(model, Cf == (0.074 / Re^0.2))
	# Drag coefficient
	@NLconstraint(model, CD >= (CDA₀ / S + k * Cf * S_wet_S + CL^2 / (π * A * e)))
	# Wing weight
	@NLconstraint(model, Ww == Ww_c2 * S + Ww_c1 * N * A^1.5 * (W0 * W * S)^0.5 / τ)
	# Total weight
	@NLconstraint(model, W >= Ww + W0)
	# Weight-cruise lift
	@NLconstraint(model, W <= 0.5ρ * V^2 * CL * S)
	# Weight-takeoff lift
	@NLconstraint(model, W <= 0.5ρ * v_to^2 * cl_max * S)
	
	# Objective function (Minimization of drag)
	@NLobjective(model, Min, 1/2 * ρ * V^2 * S * CD)

	# Solve optimizer
	with_terminal() do
	optimize!(model)
	end
	
	# Optimum points
	res = Dict((name, value(var)) for (name, var) in model.obj_dict)
end;

# ╔═╡ cc27ce98-b848-4a72-b716-e65462bcf78d
md"Aerodynamic Optimization - Max L/D"

# ╔═╡ 97c50a2f-4b49-418a-a169-666cc3a840ba
begin
	function make_case(wing, htail, vtail, α)
	    # Meshing and assembly
	    wing_mesh = WingMesh(
	                    wing, 12, 6, 
	                    span_spacing = Cosine()
	                );
		htail_mesh = WingMesh(htail, [6], 4)
		vtail_mesh = WingMesh(vtail, [4], 3)
	
		aircraft = ComponentVector(wing = make_horseshoes(wing_mesh),
								htail = make_horseshoes(htail_mesh),
	                           vtail = make_horseshoes(vtail_mesh))
	
		# Freestream conditions
	    fs = Freestream(
	            alpha = α,  # Design variable: Angle of attack
	            beta  = 0.0,
	            omega = [0.,0.,0.]
	        )
	
	    # Reference values
	    refs = References(
		            speed     = airspeed,
		            density   = 1.225,
		            viscosity = 1.5e-5,
		            area      = projected_area(wing),
		            span      = span(wing),
		            chord     = mean_aerodynamic_chord(wing),
		            location  = mean_aerodynamic_center(wing)
		        )
		
	    # Solve system
	    system  = solve_case(aircraft, fs, refs);
	end
	function optimize_lift_to_drag(x, wing, htail, vtail)
		# Unpack arguments
	    α = x[1]
	
		# Solve inviscid case
		system = make_case(wing, htail, vtail, α)
	
		# Evaluate aerodynamic coefficients
	    CDi, CY, CL = farfield(system)
	
	    # Equivalent flat-plate skin-friction drag estimation
	    x_tr      = fill(0.7, length(spans(wing))) # Transition location ratios over sections
	    CDv_plate = profile_drag_coefficient(wing, x_tr_w, system.reference) + 			profile_drag_coefficient(htail, x_tr_h, system.reference) + profile_drag_coefficient(vtail, x_tr_h, system.reference)
	
		# Calculate total drag coefficient
	    CD = CDi + CDv_plate
	
		# Calculate lift-to-drag ratio
		LD = -CL / CD
		
	    return LD
	end
	optimize_lift_to_drag(x) = optimize_lift_to_drag(x, wing, htail, vtail)
	# Initial point
	x0 = [4.0]
	res_func = optimize( 
				optimize_lift_to_drag,  # Objective function
	            x0,                  	# Initial value
				Newton(),          	 	# Optimization algorithm
	            autodiff = :forward, 	# Automatic differentiation
	            Optim.Options( 			# Options
					  # extended_trace = true,
	 	              show_trace = true,
	                  store_trace = true
	            )
	        )
	# Minimum point
	α_opt = res_func.minimizer[1]
	# Equivalent flat-plate skin-friction drag estimation
	sys_o = make_case(wing, htail, vtail, α_opt)
	ffs_o = farfield(sys_o)
	cl_best = ffs_o.CL_ff
	cdi_o = ffs_o.CDi_ff
	CDv_plate = profile_drag_coefficient(wing, x_tr_w, sys_o.reference) + profile_drag_coefficient(htail, x_tr_h, sys_o.reference) + profile_drag_coefficient(vtail, x_tr_h, sys_o.reference)
end;

# ╔═╡ 3f73efc8-894d-43fb-a88f-aea367d74588
begin
	plot(CDs, CLs,
	     label  = "",
	     xlabel = L"C_{D}",
	     ylabel = L"C_L",
	     title  = "Drag Polar",
	     ls     = :solid)
	scatter!([CDv_plate + ffs_o.CDi_ff], [ffs_o.CL_ff], label = "Max (L/D)")
end

# ╔═╡ e14c6ab8-27bf-4a16-92e5-4643543ec5bb
begin
	ldmax = ffs_o.CL_ff/(CDv_plate + ffs_o.CDi_ff)
	md"Max L/D equals $ldmax"
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AeroMDAO = "ae1513eb-aba2-415b-9b88-80f66c7c7c76"
Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
AeroMDAO = "~0.3.9"
Ipopt = "~1.0.2"
JuMP = "~1.0.0"
LaTeXStrings = "~1.3.0"
Optim = "~1.6.2"
Plots = "~1.27.5"
PlutoUI = "~0.7.38"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.ASL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6252039f98492252f9e47c312c8ffda0e3b9e78d"
uuid = "ae81ac8f-d209-56e5-92de-9978fef736f9"
version = "0.1.3+0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.AeroMDAO]]
deps = ["ComponentArrays", "CoordinateTransformations", "DelimitedFiles", "DiffResults", "ForwardDiff", "Interpolations", "LabelledArrays", "LinearAlgebra", "PrettyTables", "Rotations", "Setfield", "SparseArrays", "SplitApplyCombine", "StaticArrays", "Statistics", "Test", "TimerOutputs"]
git-tree-sha1 = "9d5188a4ea0003b4813c623d9429c98774cc07ab"
uuid = "ae1513eb-aba2-415b-9b88-80f66c7c7c76"
version = "0.3.9"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "c933ce606f6535a7c7b98e1d86d5d1014f730596"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "5.0.7"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ComponentArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "Requires"]
git-tree-sha1 = "243d8b8afc829a6707bbb1cd00da868703c2ef42"
uuid = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
version = "0.11.15"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "681ea870b918e7cff7111da58791d7f718067a19"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Dictionaries]]
deps = ["Indexing", "Random"]
git-tree-sha1 = "0340cee29e3456a7de968736ceeb705d591875a2"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.20"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "dd933c4ef7b4c270aacd4eb88fa64c147492acf0"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.10.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "56956d1e4c1221000b7781104c58c34019792951"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "af237c08bda486b74318c8070adb96efa6952530"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.2"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "cd6efcf9dc746b06709df14e462f0a3fe0786b1e"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.2+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b15fc0a95c564ca2e0a7ae12c1f095ca848ceb31"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.5"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[deps.Ipopt]]
deps = ["Ipopt_jll", "MathOptInterface"]
git-tree-sha1 = "8b7b5fdbc71d8f88171865faa11d1c6669e96e32"
uuid = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
version = "1.0.2"

[[deps.Ipopt_jll]]
deps = ["ASL_jll", "Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "MUMPS_seq_jll", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "e3e202237d93f18856b6ff1016166b0f172a49a8"
uuid = "9cc047cb-c261-5740-88fc-0cf96f7bdcc7"
version = "300.1400.400+0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.JuMP]]
deps = ["Calculus", "DataStructures", "ForwardDiff", "LinearAlgebra", "MathOptInterface", "MutableArithmetics", "NaNMath", "OrderedCollections", "Printf", "SparseArrays", "SpecialFunctions"]
git-tree-sha1 = "936e7ebf6c84f0c0202b83bb22461f4ebc5c9969"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "1.0.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "fbd884a02f8bf98fd90c53c1c9d2b21f9f30f42a"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.8.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "6f14549f7760d84b2db7a9b10b88cd3cc3025730"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.14"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a970d55c2ad8084ca317a4658ba6ce99b7523571"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.12"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.METIS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1d31872bb9c5e7ec1f618e8c4a56c8b0d9bddc7e"
uuid = "d00139f3-1899-568f-a2f0-47f597d42d70"
version = "5.1.1+0"

[[deps.MUMPS_seq_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "METIS_jll", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "29de2841fa5aefe615dea179fcde48bb87b58f57"
uuid = "d7ed1dd3-d0ae-5e8e-bfb4-87a502085b8d"
version = "5.4.1+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "JSON", "LinearAlgebra", "MutableArithmetics", "OrderedCollections", "Printf", "SparseArrays", "Test", "Unicode"]
git-tree-sha1 = "779ad2ee78c4a24383887fdba177e9e5034ce207"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.1.2"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "ba8c0f8732a24facba709388c74ba99dcbfdda1e"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c6c2ed4b7acd2137b878eb96c68e63b76199d0f"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.17+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "bc0a748740e8bc5eeb9ea6031e6f050de1fc0ba2"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.6.2"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "621f4f3b4977325b9128d5fae7a8b4829a0c2222"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.4"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "88ee01b02fba3c771ac4dce0dfc4ecf0cb6fb772"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.27.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "670e559e5c8e191ded66fa9ea89c97f10376bb4c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.38"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.Quaternions]]
deps = ["DualNumbers", "LinearAlgebra", "Random"]
git-tree-sha1 = "b327e4db3f2202a4efafe7569fcbe409106a1f75"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.5.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "a167638e2cbd8ac41f9cd57282cab9b042fa26e6"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.SplitApplyCombine]]
deps = ["Dictionaries", "Indexing"]
git-tree-sha1 = "35efd62f6f8d9142052d9c7a84e35cd1f9d2db29"
uuid = "03a91e81-4c3e-53e1-a0a4-9c0c8f19dd66"
version = "1.2.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "87e9954dfa33fd145694e42337bdd3d5b07021a6"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.6.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "4f6ec5d99a28e1a749559ef7dd518663c5eca3d5"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "8d7530a38dbd2c397be7ddd01a424e4f411dcc41"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.2"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "d60b0c96a16aaa42138d5d38ad386df672cb8bd8"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.16"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─8a679de5-f08f-44c1-bd14-7bac0e58da90
# ╟─a550479a-df78-4b09-b213-d0c237d052b7
# ╟─dce5c906-6425-471a-9b89-114feaf02a38
# ╠═8b6d0ff5-c7c2-465f-8de3-779441fc13cc
# ╟─6e4db00a-e601-4519-aad2-0f456be3ba6a
# ╟─0fc8cd9d-d268-49a3-bef5-81bc205c9dac
# ╟─64b849b4-e761-4e50-95c3-996aea337e27
# ╟─06d14669-17dd-42e6-9287-02cc518a796c
# ╟─81bb09b6-ebbb-4e90-a119-a62ba10d7a68
# ╟─21394943-8cf5-4ee0-9027-7800aff65649
# ╟─3f73efc8-894d-43fb-a88f-aea367d74588
# ╟─e14c6ab8-27bf-4a16-92e5-4643543ec5bb
# ╟─30456a2d-904c-459a-a48b-92b5e2a5d030
# ╟─20b2fccf-6c9e-477a-b795-2dacd7d1d497
# ╟─9e0099cc-f587-4bc5-b6f5-f60912069739
# ╟─72fd45c4-ce91-4f49-9a4f-b6e695fc405b
# ╟─c6a480e4-a43d-4a1e-a0a3-01269a414272
# ╟─9b5900f3-de26-48e0-887b-e40d31c9ba77
# ╟─59da49fb-d572-4a29-94f7-1e9532cf2211
# ╟─2ba9a5d6-9b30-4b62-a8ad-e0271bf519bb
# ╟─3bdb49c0-de68-4b83-ba2f-02688c91a899
# ╟─176509dc-b4f0-11ec-08a6-317c4af234c4
# ╟─542efc07-818e-42ca-a667-3d13b7d40756
# ╟─2a1e745d-dce0-464a-ae23-cf94beb1767b
# ╟─18e0076a-c70e-452b-a1f6-d0dad0303a39
# ╟─65a3845f-c6c5-4d86-8be4-d83c34800d1f
# ╟─5644b1b4-12de-4628-8a07-a03a27244025
# ╟─136e6b80-0f1b-43d5-b35e-2a82b8777cf8
# ╟─24a8d3fc-61fe-42f8-8173-5e8557ae28ad
# ╟─fd92271f-cdfc-49c6-8487-ab7d57b04d28
# ╟─5ab3b285-47db-49a5-854f-c45d8023f8c4
# ╟─beaf63ea-08cd-4ad9-819e-2fb18de27bc0
# ╟─a2e12b5b-91d1-496b-a327-c3bbf2ab897f
# ╟─68ae084d-a483-491b-8f84-fa8eb445fa6a
# ╟─d40c579e-7bc2-4d3d-9473-80a966676d74
# ╟─15e7100e-8b71-40b8-9535-dd03f6c57ede
# ╟─a70b8c36-c715-47fc-bf5d-5b95e9fc2692
# ╟─f2778e6d-1c56-449b-8236-5c8a464bfb01
# ╟─cc27ce98-b848-4a72-b716-e65462bcf78d
# ╟─97c50a2f-4b49-418a-a169-666cc3a840ba
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
