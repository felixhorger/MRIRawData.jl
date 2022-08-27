
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

	function twixread(path::AbstractString)
		!isfile(path) && error("File '$path' not found")
		local twix, twix_image, kspace
		let
			the_stdout = stdout
			the_stderr = stderr
			redirect_stdout(open("/dev/null", "w")) # Deedleduuuu Waaa
			redirect_stderr(open("/dev/null", "w")) # Print the progress to stderr, really? Amazing
			twix = siemens.mapVBVD(path)
			twix_image = twix["image"]
			twix_image.squeeze = true
			twix_image.flagRemoveOS = true
			kspace = twix_image.unsorted()
			redirect_stdout(the_stdout)
			redirect_stderr(the_stderr)
		end
		num_columns = convert(Int, twix_image["NCol"]) ÷ 2 # F me
		other_dims = convert.(
			Int,
			twix_image[key] for key in ("NLin", "NPar", "NCha")
		)
		sampling = [CartesianIndex(Int.((l, p))) for (l, p) in zip(twix_image.Lin, twix_image.Par)]
		return twix, kspace, sampling, num_columns, other_dims...
	end
	
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
		φ = atan(dc[2], dc[1])
		θ = atan(sqrt(dc[1]^2 + dc[2]^2), dc[3])
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
		return φ, θ, β, δx, Δx
	end

end

