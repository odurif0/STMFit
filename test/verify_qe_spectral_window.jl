#!/usr/bin/env julia
# Read-only saved-output verification, deliberately independent of the cube,
# frame and spectral helpers used by qe_spectral_window.jl. No fit or labels.
module VerifyQESpectralWindow
using LinearAlgebra, Statistics, SHA, TOML, Test, Printf

digest(file) = open(io -> bytes2hex(sha256(io)), file)
value(s) = parse(Float64, replace(s, 'D'=>'E', 'd'=>'e'))
function tsv(file)
    rows = readlines(file)
    keys = split(first(rows), '\t')
    [Dict(zip(keys, split(row, '\t'))) for row in rows[2:end]]
end
function cube(file)
    open(file) do io
        readline(io); readline(io)
        header = split(readline(io)); atoms = parse(Int, header[1])
        @assert atoms >= 0 && length(header) == 4
        origin = value.(header[2:4])
        dims = Int[]; axes = zeros(3,3)
        for k in 1:3
            row = split(readline(io)); push!(dims, parse(Int,row[1]))
            axes[:,k] = value.(row[2:4])
        end
        @assert all(>(0), dims)
        for _ in 1:atoms; @assert length(split(readline(io)))==5; end
        raw = Vector{Float64}(undef,prod(dims)); i = 0
        for line in eachline(io), token in split(line)
            i += 1; @assert i <= length(raw)
            raw[i] = value(token)
        end
        @assert i==length(raw) && all(isfinite,raw)
        # QE advances the third coordinate first. Keep this native order.
        (;origin,axes,dims,raw,grid=reshape(raw,dims[3],dims[2],dims[1]))
    end
end
function sample(c, xyz_bohr)
    q = c.axes \ (xyz_bohr-c.origin)
    lo = floor.(Int,q); fractions=q-lo
    @assert all(0 .<= lo) && all(lo .+ 1 .< c.dims)
    result = 0.0
    for dx in 0:1, dy in 0:1, dz in 0:1
        ds = [dx,dy,dz]
        w = prod(ds[k]==0 ? 1-fractions[k] : fractions[k] for k in 1:3)
        result += w*c.grid[lo[3]+dz+1,lo[2]+dy+1,lo[1]+dx+1]
    end
    result
end
function element(text, tag)
    hits = collect(eachmatch(Regex("<"*tag*">(.*?)</"*tag*">","s"),text))
    @assert length(hits)==1
    only(hits).captures[1]
end

function verify(run)
    @assert VERSION.major==1 && VERSION.minor==13
    config=TOML.parsefile(joinpath(run,"settings.toml"))
    @assert config["model"]["spectral_window"]=="sharp_zero_temperature"
    @assert config["preprocessing"]["cube_order"]=="qe_last_axis_fast"
    @assert !config["preprocessing"]["clip_negative_values"]
    # Independent conversion from the SI constants in QE 7.4.1.
    hartree_ev = 4.3597447222071e-18 / 1.602176634e-19
    @testset "Independent saved spectral outputs" begin
        for molecule in ("glcn","glcnac")
            dir=joinpath(run,molecule); meta=TOML.parsefile(joinpath(dir,"metadata.toml"))
            xml=read(joinpath(dir,"data-file-schema.xml"),String)
            convergence=element(element(xml,"output"),"convergence_info")
            scf=element(convergence,"scf_conv")
            @test strip(element(convergence,"wf_collected"))=="true"
            @test strip(element(scf,"convergence_achieved"))=="true"
            @test 0<=2value(element(scf,"scf_error"))<=config["model"]["scf_acceptance_ry"]
            bands=element(element(xml,"output"),"band_structure")
            eigen=value.(split(only(collect(eachmatch(r"<eigenvalues[^>]*>(.*?)</eigenvalues>"s,bands))).captures[1]))
            fermi=value(element(bands,"fermi_energy"))
            low,high=minmax(fermi,fermi+config["model"]["sample_bias_ry"]/2)
            selected=findall(e->low<=e<=high,eigen)
            @test length(selected)==meta["ildos_bands"]
            @test 2length(selected)==meta["ildos_integral_expected"]
            @test meta["emin_ev"]≈low*hartree_ev rtol=1e-14
            @test meta["emax_ev"]≈high*hartree_ev rtol=1e-14
            @test meta["xml_sha256"]==digest(joinpath(dir,"data-file-schema.xml"))
            @test meta["config_sha256"]==digest(joinpath(run,"settings.toml"))
            weights=tsv(joinpath(dir,"bands.tsv"))
            @test length(weights)==length(eigen)
            @test all(parse(Int,r["band"])==i && value(r["ildos_weight"])==(i in selected ? 2 : 0)
                for (i,r) in enumerate(weights))
            input=read(joinpath(dir,"ildos.in"),String)
            for (key,expected) in (("emin",low*hartree_ev),("emax",high*hartree_ev))
                actual=value(match(Regex(key*"\\s*=\\s*([^\\s]+)"),input).captures[1])
                @test actual≈expected rtol=1e-14
            end
            @test digest(joinpath(dir,"control.cube"))==digest(joinpath(run,molecule*"_reference.cube"))
            fields=Dict(r["key"]=>value.(split(r["value"],',')) for r in tsv(joinpath(run,molecule*"_frame.tsv")))
            normal=cross(fields["t_axis"],fields["u_axis"]); normal/=norm(normal)
            grid=collect(-8:8)*config["preprocessing"]["step_nm"]
            @test first(grid)==-config["preprocessing"]["half_nm"]
            stats=tsv(joinpath(dir,"comparison","cube_summary.tsv"))
            planes=tsv(joinpath(dir,"comparison","planes.tsv"))
            @test length(planes)==3*17*17
            for (k,observable) in enumerate(("control","ildos"))
                c=cube(joinpath(dir,observable*".cube")); row=stats[k]
                voxel=abs(det(c.axes)); integral=sum(c.raw)*voxel
                negative=count(<(0),c.raw)
                @test row["observable"]==observable
                @test parse(Int,row["values"])==length(c.raw)
                @test parse(Int,row["negative"])==negative
                @test value(row["minimum"])==minimum(c.raw)
                @test value(row["maximum"])==maximum(c.raw)
                @test value(row["integral"])≈integral rtol=1e-12
                @test value(row["negative_abs_integral"])≈-sum(min(0,v) for v in c.raw)*voxel rtol=1e-12
                observable=="ildos" && @test abs(integral-2length(selected))<=config["selection"]["ildos_integral_rtol"]*2length(selected)
                j=0
                for h in config["preprocessing"]["diagnostic_heights_nm"], u in grid, t in grid
                    j+=1; row=planes[j]
                    point=(fields["origin_nm"]+t*fields["t_axis"]+u*fields["u_axis"]+h*normal)/0.05291772109
                    @test value(row["height_nm"])==h
                    @test parse(Int,row["pixel"])==1+mod(j-1,289)
                    @test value(row[observable])≈sample(c,point) rtol=1e-10 atol=1e-15
                end
                @printf("%s %s: integral %.12g; negatives %d/%d; min %.9g; max %.9g\n",
                    molecule,observable,integral,negative,length(c.raw),minimum(c.raw),maximum(c.raw))
            end
        end
    end
end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    if ARGS==["--help"]
        println("verify_qe_spectral_window.jl RUN_DIR  # read-only, Julia 1.13")
    else
        length(ARGS)==1 || error("Expected one saved run directory")
        VerifyQESpectralWindow.verify(only(ARGS))
    end
end
