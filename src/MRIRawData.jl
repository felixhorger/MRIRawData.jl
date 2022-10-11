
module MRIRawData

	using PyCall

	function __init__()
		global siemens = pyimport("mapvbvd")
	end

	# Music of the future, praying that the ISMRMRD format is good
	#struct MRIRawDataSet{T} where T <: Number
	#	kspace::AbstractArray{T, 3} # Readout, channel, phase encoding
	#	shape::NTuple{4, Int64}
	#end

	function readtwix(path::AbstractString; quiet=true)
		!isfile(path) && error("File '$path' not found")
		twix = siemens.mapVBVD(path; quiet)
		return twix
	end
	function twix_kspace(twix; quiet::Bool=true, key::AbstractString="image")
		# k-space
		twix_obj = twix[key]
		kspace = twix_obj.unsorted(;quiet)
		return kspace
	end
	function twix_remove_oversampling(twix, b::Bool=true; key::AbstractString="image")
		twix[key].flagRemoveOS = b
		return
	end
	function twix_size(twix; key::AbstractString="image")
		twix_obj = twix[key]
		# Columns need extra care because of oversampling
		num_columns	= convert(Int, twix_obj["NCol"])
		if twix_obj.flagRemoveOS
			num_columns ÷= 2
		end
		# Read other dimensions
		#=
			Would be too easy if one could use

			num_lines		= convert(Int, twix_obj["NLin"])
			num_partitions	= convert(Int, twix_obj["NPar"])
			num_channels	= convert(Int, twix_obj["NCha"])

			but it happens that this doesn't match with the FOV that was selected, i.e.
			this gives you the sampled k-space size, not necessarily matching the reconstructed one
			I will refrain from using the language I want to use.
		=#
		bogus = twix["hdr"]["Meas"]
		num_lines		= convert(Int, bogus["lPhaseEncodingLines"])
		num_partitions	= convert(Int, bogus["lPartitions"])
		num_channels	= convert(Int, twix_obj["NCha"])
		return num_columns, num_lines, num_partitions, num_channels
	end
	function twix_sampling(twix; key::AbstractString="image")
		twix_obj = twix[key]
		return [CartesianIndex(Int.((l, p)) .+ 1) for (l, p) in zip(twix_obj.Lin, twix_obj.Par)]
	end
	function twix_index(twix, name::Symbol; key::AbstractString="image")
		twix_obj = twix[key]
		num_index = convert.(Int, twix_obj["N" * String(name)])
		indices = 1 .+ Int.(getproperty(twix_obj, name))
		return indices, num_index
	end

	for (param, twix_name, unit) in (("TE", "alTE", 0.001), ("TR", "alTR", 0.001), ("FA", "adFlipAngleDegree", π / 180.0))
		symb = Symbol(param)
		func_name = Symbol("twix_" * param)
		eval(Expr(:function,
			Expr(:call, func_name, :twix, :(i::Integer)),
			:(twix["hdr"]["MeasYaps"][$twix_name, string(i-1)] * $unit)
		))
		eval(Expr(:function, Expr(:call, func_name, :twix), quote
			bogus = twix["hdr"]["MeasYaps"]
			$symb = Vector{Float64}(undef, 0)
			i = 0
			while true
				key = ($twix_name, string(i))
				v = get(bogus, key, 0) * $unit
				v == 0 && break
				push!($symb, v)
				i += 1
			end
			return $symb
		end))
	end

	twix_afi_TR_ratio(twix) = twix["hdr"]["MeasYaps"][("sWipMemBlock", "alFree", "10")]

	
	"""
		Returns dc, β, fov, Δx
		fov[1] (readout) is set according to twix["image"].flagRemoveOS
		Δx is translation of the centre of volume.
		To get the origin, do origin = Δx .- R * (0.5 .* fov)
		where R is the rotation matrix of the volume obtained using MRICoordinates.ras (RAS coordinates)
		β is inplane rotation
	"""
	function twix_coordinates(twix; key::AbstractString="image")
		bogus = twix["hdr"]["MeasYaps"]
		keys = ("sSliceArray", "asSlice", "0") # Lord have mercy
		# Scaling
		fov = [ bogus[(keys..., dim)] for dim in ("dReadoutFOV", "dPhaseFOV", "dThickness") ]
		if !twix[key].flagRemoveOS
			fov[1] *= 2
		end
		# Rotation
		dc = Vector{Float64}(undef, 3); # Direction cosines
		orientations = ("dSag", "dCor", "dTra")
		for i = 1:3
			key = (keys..., "sNormal", orientations[i])
			dc[i] = get(bogus, key, 0.0)
		end
		# In-plane rotation
		β = let
			key = (keys..., "dInPlaneRot")
			get(bogus, key, 0.0)
		end
		# Translation
		Δx = Vector{Float64}(undef, 3)
		for i = 1:3
			key = (keys..., "sPosition", orientations[i])
			Δx[i] = get(bogus, key, 0.0)
		end
		return dc, β, fov, Δx
	end
end

