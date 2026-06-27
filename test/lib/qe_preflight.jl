module QEPreflight

export PreflightRow, check_dir!

using Printf

const PreflightRow = NamedTuple

function _read(path::String)
    isfile(path) || error("Missing required file: $path")
    return read(path, String)
end

function _match1(re, text::String, label::String)
    m = match(re, text)
    m === nothing && error("Could not parse $label")
    return strip(m.captures[1])
end

function _parse_atomic_positions(text::String)
    lines = split(text, '\n')
    idx = findfirst(line -> startswith(strip(line), "ATOMIC_POSITIONS"), lines)
    idx === nothing && error("Missing ATOMIC_POSITIONS")
    atoms = String[]
    flags = 0
    for line in lines[(idx+1):end]
        s = strip(line)
        isempty(s) && break
        startswith(s, "K_POINTS") && break
        startswith(s, "CELL_PARAMETERS") && break
        startswith(s, "ATOMIC_SPECIES") && break
        parts = split(s)
        length(parts) >= 4 || break
        push!(atoms, parts[1])
        length(parts) >= 7 && parts[end-2:end] == ["0", "0", "0"] && (flags += 1)
    end
    isempty(atoms) && error("No atom rows after ATOMIC_POSITIONS")
    return atoms, flags
end

function _parse_pseudo_dir(text::String)
    value = _match1(r"(?m)^\s*pseudo_dir\s*=\s*'([^']+)'", text, "pseudo_dir")
    return strip(value)
end

function _parse_atomic_species(text::String)
    lines = split(text, '\n')
    idx = findfirst(line -> startswith(strip(line), "ATOMIC_SPECIES"), lines)
    idx === nothing && error("Missing ATOMIC_SPECIES")
    pseudos = Dict{String,String}()
    for line in lines[(idx+1):end]
        s = strip(line)
        isempty(s) && break
        startswith(s, "CELL_PARAMETERS") && break
        startswith(s, "ATOMIC_POSITIONS") && break
        startswith(s, "K_POINTS") && break
        parts = split(s)
        length(parts) >= 3 || break
        pseudos[parts[1]] = parts[3]
    end
    isempty(pseudos) && error("No species rows after ATOMIC_SPECIES")
    return pseudos
end

function _parse_sbatch_tasks(text::String)
    m = match(r"(?m)^#SBATCH\s+--ntasks-per-node=(\d+)", text)
    m === nothing && error("Missing #SBATCH --ntasks-per-node")
    return parse(Int, m.captures[1])
end

function _parse_sbatch_cpus_per_task(text::String)
    m = match(r"(?m)^#SBATCH\s+--cpus-per-task=(\d+)", text)
    m === nothing && return 1
    return parse(Int, m.captures[1])
end

function _parse_sbatch_mem_mb(text::String)
    m = match(r"(?m)^#SBATCH\s+--mem(?:=|\s+)([0-9.]+)\s*([KMGTP]?B?)?", text)
    m === nothing && error("Missing #SBATCH --mem")
    value = parse(Float64, m.captures[1])
    unit = uppercase(something(m.captures[2], "M"))
    unit = isempty(unit) ? "M" : unit
    factor = if unit in ("K", "KB")
        1 / 1024
    elseif unit in ("M", "MB")
        1
    elseif unit in ("G", "GB")
        1024
    elseif unit in ("T", "TB")
        1024^2
    else
        error("Unsupported #SBATCH --mem unit: $unit")
    end
    return ceil(Int, value * factor)
end

function _report!(rows, run::String, key::String, value)
    push!(rows, (run=run, key=key, value=string(value)))
end

function _run_script(dir::String)
    full = joinpath(dir, "run_qe_mold.sbatch")
    prelim = joinpath(dir, "run_scf_pp.sbatch")
    if isfile(full)
        return full, "relax_scf_pp"
    elseif isfile(prelim)
        return prelim, "scf_pp_only"
    end
    error("Missing run script: expected $full or $prelim")
end

function check_dir!(rows, dir::String; min_mem_mb::Int=0)
    run = basename(normpath(dir))
    script_path, mode = _run_script(dir)
    required = ["pw_scf.in", "pp_ldos.in", basename(script_path)]
    mode == "relax_scf_pp" && push!(required, "pw_relax.in")
    for file in required
        isfile(joinpath(dir, file)) || error("Missing $(joinpath(dir, file))")
    end
    scf = _read(joinpath(dir, "pw_scf.in"))
    pp = _read(joinpath(dir, "pp_ldos.in"))
    sbatch = _read(script_path)

    source = mode == "relax_scf_pp" ? _read(joinpath(dir, "pw_relax.in")) : scf
    nat = parse(Int, _match1(r"(?m)^\s*nat\s*=\s*(\d+)", source, "nat"))
    ntyp = parse(Int, _match1(r"(?m)^\s*ntyp\s*=\s*(\d+)", source, "ntyp"))
    prefix = _match1(r"(?m)^\s*prefix\s*=\s*'([^']+)'", source, "prefix")
    pseudo_dir = _parse_pseudo_dir(source)
    pseudo_root = isabspath(pseudo_dir) ? pseudo_dir : normpath(joinpath(dir, pseudo_dir))
    pseudo_files = _parse_atomic_species(source)
    atoms, frozen = _parse_atomic_positions(source)
    scf_atoms, scf_frozen = _parse_atomic_positions(scf)
    species = unique(atoms)
    emin = parse(Float64, _match1(r"(?m)^\s*emin\s*=\s*([-+0-9.eEdD]+)", pp, "emin"))
    emax = parse(Float64, _match1(r"(?m)^\s*emax\s*=\s*([-+0-9.eEdD]+)", pp, "emax"))
    tasks = _parse_sbatch_tasks(sbatch)
    cpus_per_task = _parse_sbatch_cpus_per_task(sbatch)
    mem_mb = _parse_sbatch_mem_mb(sbatch)

    length(atoms) == nat || error("$dir: nat=$nat but ATOMIC_POSITIONS has $(length(atoms)) rows")
    length(scf_atoms) == nat || error("$dir: SCF ATOMIC_POSITIONS has $(length(scf_atoms)) rows, expected $nat")
    length(species) == ntyp || error("$dir: ntyp=$ntyp but positions contain $(length(species)) species")
    for elem in species
        haskey(pseudo_files, elem) || error("$dir: ATOMIC_SPECIES missing pseudo for $elem")
        path = joinpath(pseudo_root, pseudo_files[elem])
        isfile(path) || error("$dir: missing pseudopotential for $elem: $path")
    end
    scf_frozen == 0 || error("$dir: pw_scf.in should not contain relaxation flags")
    emin < emax || error("$dir: expected emin < emax, got $emin >= $emax")
    mem_mb >= min_mem_mb || error("$dir: #SBATCH --mem=$(mem_mb)MB is below required minimum $(min_mem_mb)MB")
    cpus_per_task == 1 || error("$dir: expected #SBATCH --cpus-per-task=1 for MPI-only QE, got $cpus_per_task")
    occursin("QE_NTASKS=", sbatch) || error("$dir: sbatch missing explicit QE_NTASKS")
    occursin("srun -n \"\$QE_NTASKS\"", sbatch) || error("$dir: sbatch missing explicit srun -n QE_NTASKS")
    if mode == "relax_scf_pp"
        occursin("extract_qe_relaxed_xyz.jl", sbatch) || error("$dir: sbatch missing extract_qe_relaxed_xyz.jl handoff")
        occursin("update_qe_positions_from_xyz.jl", sbatch) || error("$dir: sbatch missing update_qe_positions_from_xyz.jl handoff")
        occursin("pw.x -in pw_relax.in", sbatch) || error("$dir: sbatch missing relax pw.x command")
    else
        occursin("pw.x -in pw_scf.in", sbatch) || error("$dir: preliminary sbatch missing SCF pw.x command")
        !occursin("pw.x -in pw_relax.in", sbatch) || error("$dir: preliminary sbatch must not run pw_relax.in")
        !occursin("extract_qe_relaxed_xyz.jl", sbatch) || error("$dir: preliminary sbatch must not extract relaxed geometry")
        !occursin("update_qe_positions_from_xyz.jl", sbatch) || error("$dir: preliminary sbatch must not update from relaxed geometry")
    end
    occursin("pp.x -in pp_ldos.in", sbatch) || error("$dir: sbatch missing pp.x command")

    _report!(rows, run, "dir", dir)
    _report!(rows, run, "mode", mode)
    _report!(rows, run, "sbatch", basename(script_path))
    _report!(rows, run, "prefix", prefix)
    _report!(rows, run, "nat", nat)
    _report!(rows, run, "ntyp", ntyp)
    _report!(rows, run, "species", join(sort(species), ","))
    _report!(rows, run, "pseudo_dir", pseudo_dir)
    _report!(rows, run, "pseudo_files", join(["$elem=$(pseudo_files[elem])" for elem in sort(species)], ","))
    _report!(rows, run, "frozen_relax_atoms", frozen)
    _report!(rows, run, "emin_ev", @sprintf("%.6g", emin))
    _report!(rows, run, "emax_ev", @sprintf("%.6g", emax))
    _report!(rows, run, "ntasks_per_node", tasks)
    _report!(rows, run, "cpus_per_task", cpus_per_task)
    _report!(rows, run, "mem_mb", mem_mb)
    _report!(rows, run, "mem_per_task_mb", @sprintf("%.6g", mem_mb / max(tasks, 1)))
    _report!(rows, run, "status", "OK")
    return tasks
end

end
