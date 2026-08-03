const ControlType = Symbol
const CONTROL_NORMAL = :normal
const CONTROL_NO_MOLECULE = :no_molecule
const CONTROL_IDENTICAL_MOLDS = :identical_molds
const CONTROL_SWAPPED_TYPES = :swapped_types
const CONTROL_MISSING_BWD = :missing_bwd
const CONTROL_CORRUPTED_VIEW = :corrupted_view
const ALL_CONTROLS = (
    CONTROL_NORMAL, CONTROL_NO_MOLECULE, CONTROL_IDENTICAL_MOLDS,
    CONTROL_SWAPPED_TYPES, CONTROL_MISSING_BWD, CONTROL_CORRUPTED_VIEW,
)

struct ProxySite
    t::Float64
    u::Float64
    weight::Float64
    sigma_t::Float64
    sigma_u::Float64
end

Base.:(==)(a::ProxySite, b::ProxySite) =
    a.t == b.t && a.u == b.u && a.weight == b.weight &&
    a.sigma_t == b.sigma_t && a.sigma_u == b.sigma_u

struct ProxyEnsemble
    sites0::Vector{ProxySite}
    sites1::Vector{ProxySite}
end

Base.:(==)(a::ProxyEnsemble, b::ProxyEnsemble) = a.sites0 == b.sites0 && a.sites1 == b.sites1

Base.@kwdef struct SimulatorConfig
    n_min::Int = 2
    n_max::Int = 12
    width::Int = 128
    height::Int = 64
    range_nm::Tuple{Float64,Float64} = (8.0, 4.0)
    spacing_nm::Float64 = 0.52
    sigma_par_nm::Float64 = 0.20
    sigma_perp_nm::Float64 = 0.14
    amplitude_nm::Float64 = 0.35
    contrast_strength::Float64 = 0.08
    noise_sigma::Float64 = 0.03
    correlated_noise_frac::Float64 = 0.30
    blur_sigma_px::Float64 = 0.8
    affine_scale_jitter::Float64 = 0.03
    drift_strength::Float64 = 0.02
    row_offset_sigma::Float64 = 0.01
    seed::Int = 20260710
end

Base.:(==)(a::SimulatorConfig, b::SimulatorConfig) =
    a.n_min == b.n_min && a.n_max == b.n_max && a.width == b.width &&
    a.height == b.height && a.range_nm == b.range_nm &&
    a.spacing_nm == b.spacing_nm && a.sigma_par_nm == b.sigma_par_nm &&
    a.sigma_perp_nm == b.sigma_perp_nm && a.amplitude_nm == b.amplitude_nm &&
    a.contrast_strength == b.contrast_strength && a.noise_sigma == b.noise_sigma &&
    a.correlated_noise_frac == b.correlated_noise_frac &&
    a.blur_sigma_px == b.blur_sigma_px && a.affine_scale_jitter == b.affine_scale_jitter &&
    a.drift_strength == b.drift_strength && a.row_offset_sigma == b.row_offset_sigma &&
    a.seed == b.seed

struct ViewNuisanceTruth
    scale_x::Float64
    scale_y::Float64
    shear::Float64
    blur_sigma_px::Float64
end

Base.:(==)(a::ViewNuisanceTruth, b::ViewNuisanceTruth) =
    a.scale_x == b.scale_x && a.scale_y == b.scale_y &&
    a.shear == b.shear && a.blur_sigma_px == b.blur_sigma_px

struct SyntheticTruth
    N::Int
    sequence::Vector{Int}
    lobe_x_nm::Vector{Float64}
    lobe_y_nm::Vector{Float64}
    amplitudes::Vector{Float64}
    sigma_par::Vector{Float64}
    sigma_perp::Vector{Float64}
    orientation_deg::Float64
    direction::Int
    phase::Int
    mirror::Int
    fwd_nuisance::ViewNuisanceTruth
    bwd_nuisance::ViewNuisanceTruth
    control::ControlType
end

Base.:(==)(a::SyntheticTruth, b::SyntheticTruth) =
    a.N == b.N && a.sequence == b.sequence && a.lobe_x_nm == b.lobe_x_nm &&
    a.lobe_y_nm == b.lobe_y_nm && a.amplitudes == b.amplitudes &&
    a.sigma_par == b.sigma_par && a.sigma_perp == b.sigma_perp &&
    a.orientation_deg == b.orientation_deg && a.direction == b.direction &&
    a.phase == b.phase && a.mirror == b.mirror &&
    a.fwd_nuisance == b.fwd_nuisance && a.bwd_nuisance == b.bwd_nuisance &&
    a.control == b.control

struct SyntheticCase
    case_id::String
    img::SXMImage
    truth::SyntheticTruth
    control::ControlType
    case_seed::UInt64
end

Base.:(==)(a::SyntheticCase, b::SyntheticCase) =
    a.case_id == b.case_id && a.control == b.control && a.case_seed == b.case_seed &&
    image_data_bytes(a) == image_data_bytes(b) && a.truth == b.truth

function image_data_bytes(case::SyntheticCase)
    io = IOBuffer()
    write(io, codeunits(case.case_id))
    write(io, codeunits(String(case.control)))
    write(io, case.case_seed)
    write(io, codeunits(case.img.filepath))
    write(io, case.img.width)
    write(io, case.img.height)
    write(io, case.img.range_nm[1])
    write(io, case.img.range_nm[2])
    write(io, case.img.offset_nm[1])
    write(io, case.img.offset_nm[2])
    header_keys = sort!(collect(keys(case.img.header)))
    write(io, length(header_keys))
    for key in header_keys
        write(io, codeunits(key))
        write(io, codeunits(case.img.header[key]))
    end
    write(io, length(case.img.channels))
    for ch in case.img.channels
        write(io, codeunits(ch.name))
        write(io, codeunits(ch.unit))
        write(io, codeunits(ch.direction))
        write(io, size(ch.data, 1))
        write(io, size(ch.data, 2))
        for v in ch.data
            write(io, v)
        end
    end
    return take!(io)
end

case_checksum(case::SyntheticCase) = bytes2hex(sha256(image_data_bytes(case)))
