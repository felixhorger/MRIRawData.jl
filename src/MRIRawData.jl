
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
		twix_image = twix["image"]
		twix_image.squeeze = true
		twix_image.flagRemoveOS = true
		kspace = twix_image.unsorted(;quiet)
		num_columns		= convert(Int, twix_image["NCol"]) ÷ 2 # F me
		num_lines		= convert.(Int, twix_image["NLin"])
		num_partitions	= convert.(Int, twix_image["NPar"])
		num_channels	= convert.(Int, twix_image["NCha"])
		sampling = [CartesianIndex(Int.((l+1, p+1))) for (l, p) in zip(twix_image.Lin, twix_image.Par)]
		return twix, kspace, sampling, num_columns, num_lines, num_partitions, num_channels
	end
	
	"""
		Returns dc, β, δx, Δx
		Δx is translation of the centre of volume.
		To get the origin, do origin = Δx .- R * (0.5 * δx .* [num_columns, num_lines, num_partitions])
		where R is the rotation matrix of the volume obtained using MRICoordinates.ras (RAS coordinates)
	"""
	function twix_coordinates(twix)
		twix_image = twix["image"]
		bogus = twix["hdr"]["MeasYaps"]
		keys = ("sSliceArray", "asSlice", "0") # Lord have mercy
		# Scaling
		δx = [
			(bogus[(keys..., "dReadoutFOV")] * 2) / convert(Int, twix_image["NCol"]),
			bogus[(keys..., "dPhaseFOV")]	/ convert(Int, twix_image["NLin"]),
			bogus[(keys..., "dThickness")]	/ convert(Int, twix_image["NPar"])
		]
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
		return dc, β, δx, Δx
	end

end

