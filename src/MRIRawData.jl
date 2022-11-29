
module MRIRawData

	import PyCall
	import Base: size, show
	import MRICoordinates

	function __init__()
		global siemens = PyCall.pyimport("mapvbvd")
	end

	"""
		Wrapper for Siemens raw data.
		This name is more informative than "twix".
		Some of the internal code still uses the name twix to separate it from other formats.
	"""
	struct SiemensRawData
		data::Dict{String, Any}
	end

	show(io::IO, _::SiemensRawData) = print("Siemens MRI raw data")
	show(io::IO, ::MIME"text/plain", raw::SiemensRawData) = show(io, raw)
	load_siemens(path::AbstractString; quiet=true) = SiemensRawData(siemens.mapVBVD(path; quiet))

	function get_kspace(raw::SiemensRawData; quiet::Bool=true, key::String="image")
		# k-space
		twix_obj = raw.data[key]
		kspace = twix_obj.unsorted(;quiet)
		return kspace
	end
	function remove_oversampling!(raw::SiemensRawData, b::Bool=true; key::String="image")
		raw.data[key].flagRemoveOS = b
		return
	end
	function size(raw::SiemensRawData; key::String="image")
		twix = raw.data
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
			this gives you the sampled k-space size, not necessarily matching the reconstructed one.
		=#
		bogus = twix["hdr"]["Meas"]
		num_lines		= convert(Int, bogus["lPhaseEncodingLines"])
		num_partitions	= convert(Int, bogus["lPartitions"])
		num_channels	= convert(Int, twix_obj["NCha"])
		return num_columns, num_lines, num_partitions, num_channels
	end
	function get_sampling(raw::SiemensRawData; key::String="image")
		twix_obj = raw.data[key]
		return [CartesianIndex(Int.((l, p)) .+ 1) for (l, p) in zip(twix_obj.Lin, twix_obj.Par)]
	end
	function get_sampling_index(raw::SiemensRawData, name::Symbol; key::String="image")
		twix_obj = raw.data[key]
		num_index = convert.(Int, twix_obj["N" * String(name)])
		indices = 1 .+ Int.(getproperty(twix_obj, name))
		return indices, num_index
	end

	for (param, twix_name, unit) in (
		("TE", "alTE", 0.001),
		("TR", "alTR", 0.001),
		("FA", "adFlipAngleDegree", π / 180.0)
	)
		symb = Symbol(param)
		func_name = Symbol("get_" * param)
		eval(Expr(:function,
			Expr(:call, func_name, :(raw::SiemensRawData), :(i::Integer)),
			:(raw.data["hdr"]["MeasYaps"][$twix_name, string(i-1)] * $unit)
		))
		eval(Expr(:function, Expr(:call, func_name, :(raw::SiemensRawData)), quote
			bogus = raw.data["hdr"]["MeasYaps"]
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

	"""
		Dwell time in μs
	"""
	function get_dwelltime(raw::SiemensRawData)
		# Get array of dwell times (why array?)
		bogus = raw.data["hdr"]["Meas"]["alDwellTime"] # I don't want to know why this is a string
		# Find the first blank and get first "word"
		bogus = bogus[1:findfirst(" ", bogus).start-1]
		# Convert to number and choose proper unit
		return parse(Float64, bogus) / 1000.0
	end

	
	"""
		Returns normal, β, fov, Δx
		normal ≡ partition direction in patient coordinate system (LPS coordinates)
		fov[1] (readout) is set according to twix["image"].flagRemoveOS
		Δx is translation of the centre of volume in patient coordinates.
		β is inplane rotation, clockwise
	"""
	function get_coordinates(raw::SiemensRawData; key::String="image")
		twix = raw.data
		bogus = twix["hdr"]["MeasYaps"]
		keys = ("sSliceArray", "asSlice", "0") # Lord have mercy
		# Scaling
		fov = [ bogus[(keys..., dim)] for dim in ("dReadoutFOV", "dPhaseFOV", "dThickness") ]
		if !twix[key].flagRemoveOS
			fov[1] *= 2
		end
		# Rotation
		normal = Vector{Float64}(undef, 3); # Direction cosines
		orientations = ("dSag", "dCor", "dTra")
		for i = 1:3
			key = (keys..., "sNormal", orientations[i])
			normal[i] = get(bogus, key, 0.0)
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
		return normal, β, fov, Δx
	end

	"""
		Patient position(ing)
	"""
	function get_patient_position(raw::SiemensRawData)
		pos_raw = raw.data["hdr"]["Config"]["PatientPosition"]
		local pos::Orientation
		if pos_raw == "HFS"
			pos = MRICoordinates.HeadFirstSupine
		elseif pos_raw == "HFP"
			pos = MRICoordinates.HeadFirstProne
		elseif pos_raw == "HFLR"
			pos = MRICoordinates.HeadFirstLateralRight
		elseif pos_raw == "HFLL"
			pos = MRICoordinates.HeadFirstLateralLeft
		elseif pos_raw == "FFS"
			pos = MRICoordinates.FeetFirstSupine
		elseif pos_raw == "FFP"
			pos = MRICoordinates.FeetFirstProne
		elseif pos_raw == "FFLR"
			pos = MRICoordinates.FeetFirstLateralRight
		elseif pos_raw == "FFLL"
			pos = MRICoordinates.FeetFirstLateralLeft
		else
			error("Could not determine patient position from \"$raw\"")
		end
		return pos
	end

	"""
		types: 0 for integer, 1 for double
		The so called "WIP" parameters
	"""
	function get_special_parameters(raw::SiemensRawData)
		num_WIPparameters = 14
		bogus = raw.data["hdr"]["MeasYaps"]
		params = Vector{Union{Int64, Float64}}(undef, num_WIPparameters)
		for i = 1:num_WIPparameters
			j = i - 1
			key = ("sWipMemBlock", "alFree", "$j")
			if haskey(bogus, key)
				params[i] = bogus[key]
			else
				params[i] = get(bogus, ("sWipMemBlock", "adFree", "$j"), 0)
			end
		end
		return params
	end
	# Special case, convenient
	get_afi_TR_ratio(raw::SiemensRawData) = raw.data["hdr"]["MeasYaps"][("sWipMemBlock", "alFree", "10")]

end

