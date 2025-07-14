using Roots, QuadGK

### Association model in the form of a function to get cluster concentration from total conentration. -->

# # Isodesmic
# function cn_from_ctot(
	# n::Int, # []
	# ctot::Float64; # [1/um^3]
	# K::Float64 = 5.0E-3) # [um^3]
	# #K::Float64 = 3.0/(1E-6*6.022E23*1E-15)) # [um^3]
	
	# c1(ctot) = ( 1 + 2*K*ctot - sqrt(1 + 4*K*ctot) ) / ( 2*K^2*ctot ) # [1/um^3]->[1/um^3]
	
	# cn = exp( (n-1)*log(K) + n*log(c1(ctot)) )
	
	# return cn::Float64 # [1/um^3]
	
# end

# No clustering
function cn_from_ctot(
	n::Int, # []
	ctot::Float64) # [1/um^3]
	
	if n == 1
		cn = ctot
	else
		cn = 0.0
	end
	
	return cn::Float64 # [1/um^3]
	
end

### <--

# The empirical relationship between surface pressure and adsorbed amount (gamma).
function surface_pressure(
	gamma::Float64; # [1/um^2]
	y0::Float64 = 72.0E-3, # [J/m^2]
	gamma0::Float64 = 24.0E3) # [1/um^2]
	
	p1 = 1.029
	p2 = -32.4
	p3 = 77.4
	p4 = 49.8
	
	f(gamma) = y0 - 1E-3*log(p1^(p2*gamma/gamma0 + p3) + p1^p4)/log(p1)
	
	pressure = max(0.0, f(gamma)) # [J/m^2]
	
	return pressure::Float64 # [J/m^2]
	
end

# The adsorption rate constant (dgamma/dt * (1/c)) to the interface given the already adsorbed amount (gamma).
function adsorption_rate_constant(
	n::Int, # []
	gamma::Float64; # [1/um^2]
	B::Float64 = 70.0E3, # [1/(s*um^2)]
	q::Float64 = 0.15, # []
	kT::Float64 = 4.11E-21, # [J]
	gamma0::Float64 = 24.0E3) # [1/um^2]
	
	function cavitation_free_energy(gamma::Float64) # [1/um^2]->[J]
		
		if gamma == 0
			G = 0.0 # [J]
		else
			G = (surface_pressure.(gamma; gamma0 = gamma0) ./ (1E12*gamma)) # [J]
		end
		
		return G::Float64 # [J]
		
	end
	
	function compression_free_energy(gamma::Float64) # [1/um^2]->[J]
		
		if gamma <= gamma0
			G = 0.0 # [J]
		else
			G = -quadgk(A -> surface_pressure(1/(1E12*A); gamma0 = gamma0), 1/(1E12*gamma0), 1/(1E12*gamma))[1] # [J]
		end
		
		return G::Float64 # [J]
		
	end
	
	rateconstant = n*B*exp(-q*n^(1/20)*(compression_free_energy(gamma) + cavitation_free_energy(gamma))/kT) # [um/s]
	
	return rateconstant::Float64 # [1/(s*um^2)]
	
end

# Numerically simulates the adsorption process (in a sphere) by combining the adsorption rate model and diffusion.
function simulate(
	cbulk::Float64, # [1/um^3]
	tmax::Float64, # [s]
	xmax::Float64; # [um] (drop radius)
	dt::Float64 = 0.01, # [s]
	tinterval::Float64 = 1.0, # [s]
	dx::Float64 = 5.0, # [um]
	y0::Float64 = 72.0E-3, # [J/m^2]
	diffcoeffmonomer::Float64 = 45.0, # [um^2/s]
	nmax::Int = 30) # [] (cluster size cutoff in simulation)
	
	# dc/dx
	function first_derivative_c_x(n)
		
		dcs_dx = Vector{Float64}(undef, length(cns[n, :]))
		
		dcs_dx[1] = (cns[n, 2] - ccenter[n]) / (2*dx)
		for i in eachindex(cns[n, :])[2:end-1]
			dcs_dx[i] = (cns[n, i+1] - cns[n, i-1]) / (2*dx)
		end
		dcs_dx[end] = (csurf[n] - cns[n, end-1]) / (2*dx)
		
		return dcs_dx::Vector{Float64}
		
	end
	
	# d^2c/dx^2
	function second_derivative_c_x(n)
		
		ddcs_dx2 = Vector{Float64}(undef, length(cns[n, :]))
		
		ddcs_dx2[1] = (cns[n, 2] - 2*cns[n, 1] + ccenter[n]) / dx^2
		for i in eachindex(cns[n, :])[2:end-1]
			ddcs_dx2[i] = (cns[n, i+1] - 2*cns[n, i] + cns[n, i-1]) / dx^2
		end
		ddcs_dx2[end] = (csurf[n] - 2*cns[n, end] + cns[n, end-1]) / dx^2
		
		return ddcs_dx2::Vector{Float64}
		
	end
	
	#diffcoeff(n) = diffcoeffmonomer*n^(-1/3) # [um^2/s] (diffusion coefficient for n-cluster)
	diffcoeff(n) = diffcoeffmonomer # [um^2/s] (diffusion coefficient for n-cluster)
	
	# Applies the physics of diffusion and adsorption to the system
	function physics!()
		
		for n in 1:nmax
			clusterflux = dt * diffcoeff(n) * (cns[n, end-1] - csurf[n]) / (2*dx) # [1/um^2]
			gamma += n * clusterflux # [1/um^2]
		end
		for n in 1:nmax
			# Sets the concentration next to the surface such that it simultaneously satisfies the diffusion flux and adsorption rate model
			csurf[n] = n * diffcoeff(n) * cns[n, end-1] / ( 2*dx * adsorption_rate_constant(n, gamma) + (n * diffcoeff(n)) ) # [1/um^3]
		end
		for n in 1:nmax
			nextccentern = dt * diffcoeff(n) * ( (cns[n, 1] - 2*ccenter[n] + cns[n, 1]) / dx^2 )
			cns[n, :] .+= dt .* diffcoeff(n) .* ( (2 ./ ((1:size(cns)[2]) .* dx)) .* first_derivative_c_x(n) .+ second_derivative_c_x(n) ) # [1/um^3]
			ccenter[n] = nextccentern # [1/um^3]
		end
		
		return nothing
		
	end
	
	ccenter::Vector{Float64} = fill(cbulk, nmax) # [1/um^3] (concentration in center of sphere for each n-cluster)
	csurf::Vector{Float64} = fill(0.0, nmax) # [1/um^3] (concentration next to the interface for each n-cluster)
	
	# Radial concentration profile for each n-cluster (rows)
	cns::Matrix{Float64} = zeros(nmax, round(Int, xmax/dx)-1)
	for n in 1:nmax
		cns[n, :] = fill(cn_from_ctot(n, cbulk), size(cns)[2]) # [1/um^3]
	end
	
	gamma::Float64 = 0.0 # [1/um^2]
	gammas::Vector{Tuple{Float64, Float64}} = []
	push!(gammas, (0.0, gamma)) # [s], [1/um^2]
	
	surfacetension::Float64 = y0 # [J/m^2]
	surfacetensions::Vector{Tuple{Float64, Float64}} = []
	push!(surfacetensions, (0.0, y0)) # [s], [J/m^2]
	
	for t in dt:dt:tmax
		
		physics!()
		
		if round(t/dt) % round(tinterval/dt) == 0
			push!(gammas, (t, gamma)) # [s], [1/um^2]
			push!(surfacetensions, (t, y0 - surface_pressure(gamma; y0 = y0))) # [s], [J/m^2]
		end
		
		### Total amount check (should stay constant) -->
		if round(t/dt) % round(1000/dt) == 0
			println(string(t)*" s of "*string(tmax)*" s")
			totalAmount = gamma*4*pi*xmax^2
			for n in 1:nmax
				totalAmount += csurf[n]*4*pi*xmax^2
				totalAmount += sum(dx*(n*cns[n, :] .* ( 4*pi*(dx*eachindex(cns[n, :])).^2 )))
			end
			println("Total protein amount: "*string(totalAmount))
		end
		### <--
		
	end
	
	return gammas::Vector{Tuple{Float64, Float64}}, surfacetensions::Vector{Tuple{Float64, Float64}}
	
end
