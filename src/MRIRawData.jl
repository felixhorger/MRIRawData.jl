
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
	function twix_kspace(twix; rm_oversampling=false, quiet=true)
		# k-space
		twix_image = twix["image"]
		twix_image.squeeze = true
		twix_image.flagRemoveOS = rm_oversampling
		kspace = twix_image.unsorted(;quiet)
		# Columns need extra care because of oversampling
		num_columns	= convert(Int, twix_image["NCol"])
		if rm_oversampling
			num_columns ÷= 2
		end
		# Read other dimensions
		num_lines		= convert.(Int, twix_image["NLin"])
		num_partitions	= convert.(Int, twix_image["NPar"])
		num_channels	= convert.(Int, twix_image["NCha"])
		num_repetitions	= convert.(Int, twix_image["NRep"])
		# Sampling indices
		sampling = [CartesianIndex(Int.((l+1, p+1))) for (l, p) in zip(twix_image.Lin, twix_image.Par)]
		return kspace, sampling, num_columns, num_lines, num_partitions, num_channels, num_repetitions
	end
	
	"""
		Returns dc, β, fov, Δx
		Δx is translation of the centre of volume.
		To get the origin, do origin = Δx .- R * (0.5 .* fov)
		where R is the rotation matrix of the volume obtained using MRICoordinates.ras (RAS coordinates)
	"""
	function twix_coordinates(twix)
		twix_image = twix["image"]
		bogus = twix["hdr"]["MeasYaps"]
		keys = ("sSliceArray", "asSlice", "0") # Lord have mercy
		# Scaling
		fov = [ bogus[(keys..., dim)] for dim in ("dReadoutFOV", "dPhaseFOV", "dThickness") ]
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

