using Test, TOML, LinearAlgebra, Statistics
include(joinpath(@__DIR__, "qe_cutoff_wavefunctions.jl"))
const X=QECutoffWavefunctions
const D=X.D; const S=X.S; const W=X.W
const ROOT=normpath(joinpath(@__DIR__,".."))
BLAS.set_num_threads(1)
num(r,k)=parse(Float64,r[k])

function verify_saved(run)
    s=X.settings(joinpath(run,"settings.toml")); heights=s["preprocessing"]["diagnostic_heights_nm"]
    per=(round(Int,2s["preprocessing"]["half_nm"]/s["preprocessing"]["step_nm"])+1)^2
    data=Dict()
    for m in X.MOLECULES, cutoff in (360,720)
        dir=joinpath(run,m,string(cutoff)); out=joinpath(dir,"analysis")
        report=TOML.parsefile(joinpath(out,"summary.toml")); meta=TOML.parsefile(joinpath(dir,"metadata.toml"))
        @test all(values(report["checks"]))
        @test !report["mold_adopted"] && !report["new_cube_validation"]
        @test report["config_sha256"]==S.sha(joinpath(run,"settings.toml"))
        @test report["xml_sha256"]==S.sha(joinpath(dir,"data-file-schema.xml"))
        @test report["wfc_sha256"]==meta["wfc_sha256"]
        @test report["point_clearance_nm"]>0 && report["vertex_clearance_nm"]>0
        @test report["parseval_smooth_norm"]>0 && isfinite(report["parseval_smooth_norm"])
        sp=S.spectrum(joinpath(dir,"data-file-schema.xml"),s)
        @test report["selected_bands"]==[r.band for r in sp.rows if r.ildos_weight>0]
        @test report["selected_weights"]==[r.ildos_weight for r in sp.rows if r.ildos_weight>0]
        bands=D.V.tsv(joinpath(out,"bands.tsv"))
        @test length(bands)==length(sp.rows)
        for (r,b) in zip(bands,sp.rows)
            @test parse(Int,r["band"])==b.band
            @test num(r,"energy_ev")==b.energy_ev && num(r,"ildos_weight")==b.ildos_weight
        end
        planes=D.V.tsv(joinpath(out,"planes.tsv")); queries=D.V.tsv(joinpath(out,"queries.tsv"))
        nv=report["native_vertices"]; data[(m,cutoff)]=planes
        @test length(planes)==per*length(heights)==report["plane_points"]
        @test length(queries)==nv+length(planes)
        @test [parse(Int,r["query"]) for r in queries]==collect(1:length(queries))
        for (i,r) in enumerate(queries)
            @test num(r,"density")==num(r,"repeat")
            @test isfinite(num(r,"density")) && num(r,"density")>=0
            @test r["kind"]==(i<=nv ? "native_vertex" : "plane")
        end
        frame=D.G.C.read_frame(joinpath(dir,"frame.tsv"))
        pts=reduce(vcat,[S.plane_points(frame,h,s) for h in heights])
        for (i,r) in enumerate(planes)
            @test (num(r,"height_nm"),parse(Int,r["pixel"]))==(heights[1+div(i-1,per)],1+mod(i-1,per))
            @test [num(r,k) for k in ("x_nm","y_nm","z_nm")]==collect(pts[i,:])
            @test num(r,"density")==num(queries[nv+i],"density")
        end
        stats=D.V.tsv(joinpath(out,"plane_summary.tsv")); tails=D.V.tsv(joinpath(out,"decay.tsv"))
        @test length(stats)==length(heights) && length(tails)==per*(length(heights)-1)
        for (j,h) in enumerate(heights)
            vals=[num(r,"density") for r in planes[(j-1)*per+1:j*per]]; r=stats[j]
            @test num(r,"height_nm")==h && parse(Int,r["points"])==per
            @test parse(Int,r["zeros"])==count(iszero,vals)
            @test num(r,"minimum")==minimum(vals) && num(r,"maximum")==maximum(vals)
            @test num(r,"median")==median(vals)
            @test isapprox(num(r,"density_sum"),sum(vals);rtol=16eps(),atol=0)
        end
        for (i,r) in enumerate(tails)
            j=1+div(i-1,per); pixel=1+mod(i-1,per)
            a=num(planes[(j-1)*per+pixel],"density"); b=num(planes[j*per+pixel],"density")
            @test num(r,"height_from_nm")==heights[j] && num(r,"height_to_nm")==heights[j+1]
            @test parse(Int,r["pixel"])==pixel && num(r,"lower_density")==a && num(r,"upper_density")==b
            @test parse(Bool,r["rising"])==(b>a)
            @test isequal(num(r,"ratio"),a>0 ? b/a : NaN)
            expected=a>0 && b>0 ? log(a/b)/(heights[j+1]-heights[j]) : NaN
            @test isequal(num(r,"log_decay_per_nm"),expected) || isapprox(num(r,"log_decay_per_nm"),expected;rtol=16eps(),atol=0)
        end
        if cutoff==360
            old=D.V.tsv(joinpath(dir,"reference_planes.tsv")); nodes=D.V.tsv(joinpath(dir,"reference_vertices.tsv"))
            @test length(old)==length(planes) && length(nodes)==nv
            @test all(num(a,"direct")==num(b,"density") for (a,b) in zip(old,planes))
            @test all(num(a,"direct")==num(b,"density") for (a,b) in zip(nodes,queries[1:nv]))
        end
    end
    pairs=D.V.tsv(joinpath(run,"comparison/paired_planes.tsv"))
    stats=D.V.tsv(joinpath(run,"comparison/summary.tsv"))
    @test length(pairs)==2per*length(heights) && length(stats)==2length(heights)
    for (midx,m) in enumerate(X.MOLECULES)
        a=data[(m,360)]; b=data[(m,720)]
        for i in eachindex(a)
            r=pairs[(midx-1)*length(a)+i]; va=num(a[i],"density"); vb=num(b[i],"density")
            @test r["molecule"]==m && r["height_nm"]==a[i]["height_nm"] && r["pixel"]==a[i]["pixel"]
            @test num(r,"density_360")==va && num(r,"density_720")==vb && num(r,"difference")==vb-va
            @test isequal(num(r,"ratio"),va>0 ? vb/va : NaN)
        end
        for (j,h) in enumerate(heights)
            ix=(j-1)*per+1:j*per; va=num.(a[ix],"density"); vb=num.(b[ix],"density")
            r=stats[(midx-1)*length(heights)+j]; ratios=[x>0 ? y/x : NaN for (x,y) in zip(va,vb)]
            @test r["molecule"]==m && num(r,"height_nm")==h && parse(Int,r["points"])==per
            @test isequal(num(r,"minimum_ratio"),minimum(ratios)) && isequal(num(r,"maximum_ratio"),maximum(ratios))
            @test isequal(num(r,"median_ratio"),median(ratios))
            expected_l2=sqrt(sum(abs2,vb-va)/sum(abs2,va))
            expected_sum=sum(va)>0 ? sum(vb)/sum(va) : NaN
            @test isequal(num(r,"relative_l2"),expected_l2) || isapprox(num(r,"relative_l2"),expected_l2;rtol=16eps(),atol=0)
            @test isequal(num(r,"sum_ratio"),expected_sum) || isapprox(num(r,"sum_ratio"),expected_sum;rtol=16eps(),atol=0)
        end
    end
end

function record(io,items...)
    b=IOBuffer(); foreach(x->write(b,x),items); bytes=take!(b); n=Int32(length(bytes))
    write(io,n); write(io,bytes); write(io,n)
end
function toy_wfc(file,coeff)
    miller=Int32[0 1 0; 0 0 0; 0 0 1]; cell=20Matrix{Float64}(I,3,3)
    open(file,"w") do io
        record(io,Int32(1),zeros(3),Int32(1),Int32(-1),1.)
        record(io,Int32.([3,3,1,3])); record(io,2pi*inv(transpose(cell))); record(io,miller)
        for col in eachcol(coeff); record(io,collect(col)); end
    end
    open(W.header,file),miller
end

function toy_xml(cutoff,dims)
    """
    <qes Units="Hartree atomic units"><creator VERSION="7.4.1"/>
    <input><atomic_species><species name="C"><mass>12.011</mass><pseudo_file>C.UPF</pseudo_file></species></atomic_species>
    <dft><functional>PBE</functional><vdw_corr>grimme-d3</vdw_corr><dftd3_version>3</dftd3_version><dftd3_threebody>true</dftd3_threebody></dft>
    <bands><tot_charge>0</tot_charge><occupations>smearing</occupations></bands>
    <spin><lsda>false</lsda><noncolin>false</noncolin><spinorbit>false</spinorbit></spin></input>
    <output><convergence_info><wf_collected>true</wf_collected><scf_conv><convergence_achieved>true</convergence_achieved><scf_error>1e-10</scf_error></scf_conv></convergence_info>
    <atomic_structure nat="1" alat="20"><atomic_positions><atom name="C" index="1">0.01 0.01 0.01</atom></atomic_positions>
    <cell><a1>20 0 0</a1><a2>0 20 0</a2><a3>0 0 20</a3></cell></atomic_structure>
    <basis_set><gamma_only>true</gamma_only><ecutwfc>25</ecutwfc><ecutrho>$(cutoff/2)</ecutrho><fft_grid nr1="$dims" nr2="$dims" nr3="$dims"></fft_grid></basis_set>
    <band_structure><nks>1</nks><nbnd>3</nbnd><nelec>2</nelec><lsda>false</lsda><noncolin>false</noncolin><spinorbit>false</spinorbit>
    <ks_energies><k_point weight="2">0 0 0</k_point><eigenvalues size="3">-0.1 0 0.1</eigenvalues></ks_energies>
    <smearing degauss="0.01">mv</smearing><fermi_energy>0.005</fermi_energy></band_structure>
    <total_energy><etot>-1</etot></total_energy></output></qes>
    """
end

@testset "Frozen scope and nonselected zero/rising samples" begin
    config=joinpath(ROOT,"config/qe_cutoff_wavefunctions.toml"); s=X.settings(config)
    @test s["model"]["ecutrho_ry"]==[360.,720.] && !s["selection"]["adopt_reference"]
    @test X.main(["--help"])===nothing
    mktempdir() do dir
        changed=deepcopy(s); changed["preprocessing"]["diagnostic_heights_nm"]=[.3,.4,.5]
        bad=joinpath(dir,"bad.toml"); open(io->TOML.print(io,changed),bad,"w")
        @test_throws ErrorException X.settings(bad)
    end
    product=X.plane_products([1.,0.,2.,0.,4.,0.],zeros(6,3),[.4,.5,.6])
    @test length(product.rows)==6 && length(product.tails)==4
    @test count(r->r.rising,product.tails)==2
    @test count(r->isnan(r.ratio),product.tails)==2
    @test all(r->r.zeros==1,product.stats)
    @test_throws ErrorException X.plane_products([1.],zeros(1,3),[.4,.5])
    job=read(joinpath(ROOT,"hpc/qe_cutoff_wavefunctions.sbatch"),String)
    @test occursin("#SBATCH --cpus-per-task=4",job) && occursin("#SBATCH --mem=32000MB",job)
    @test occursin("#SBATCH --time=01:00:00",job) && occursin("#SBATCH --no-requeue",job)
    @test !occursin("pw.x",job) && !occursin("pp.x",job) && !occursin("sbatch ",job)
end

@testset "End-to-end synthetic four-state Fourier comparison" begin
    mktempdir() do run
        cp(joinpath(ROOT,"config/qe_cutoff_wavefunctions.toml"),joinpath(run,"settings.toml"))
        s=X.settings(joinpath(run,"settings.toml")); mkpath(joinpath(run,"pseudo"))
        write(joinpath(run,"pseudo/C.UPF"),"""
        <PP_HEADER element="C" is_paw="true" mesh_size="3"/><PP_R>0.01 0.02 0.1</PP_R>
        <PP_AUGMENTATION cutoff_r_index="3"></PP_AUGMENTATION><PP_BETA.1 cutoff_radius_index="3"></PP_BETA.1>
        """)
        for m in X.MOLECULES, cutoff in (360,720)
            dir=joinpath(run,m,string(cutoff)); mkpath(dir)
            write(joinpath(dir,"data-file-schema.xml"),toy_xml(cutoff,cutoff==360 ? 32 : 40))
            write(joinpath(dir,"frame.tsv"),"key\tvalue\norigin_nm\t0.5,0.5,0.2\nt_axis\t1,0,0\nu_axis\t0,1,0\n")
            coeff=zeros(ComplexF64,3,3); coeff[1,2]=.25
            cutoff==720 && (coeff[2,2]=.02+.01im; coeff[3,2]=.03)
            wfc=joinpath(dir,"wfc.dat"); header,miller=toy_wfc(wfc,coeff)
            geo=D.xml_geometry(joinpath(dir,"data-file-schema.xml"))
            frame=D.G.C.read_frame(joinpath(dir,"frame.tsv")); heights=s["preprocessing"]["diagnostic_heights_nm"]
            pts=reduce(vcat,[S.plane_points(frame,h,s) for h in heights]); ns=D.native_vertices(geo.cell,geo.dims,pts)
            baseline=2*.25^2/header.volume
            S.table(joinpath(dir,"reference_planes.tsv"),[(;height_nm=h,pixel,direct=baseline) for h in heights for pixel in 1:289])
            S.table(joinpath(dir,"reference_vertices.tsv"),[(;ix=x,iy=y,iz=z,direct=baseline) for (x,y,z) in ns.indices])
            _,rs=D.radii(run,dir)
            meta=Dict("molecule"=>m,"ecutrho_ry"=>cutoff,"wfc_sha256"=>S.sha(wfc),"config_sha256"=>S.sha(joinpath(run,"settings.toml")),
                "pseudo_sha256"=>Dict(r.pseudo=>r.sha256 for r in rs),"input_sha256"=>Dict(f=>S.sha(joinpath(dir,f)) for f in
                    ("data-file-schema.xml","frame.tsv","reference_planes.tsv","reference_vertices.tsv")))
            open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
            # Exercise the hash guard before any analysis directory exists.
            original_frame=read(joinpath(dir,"frame.tsv"),String)
            write(joinpath(dir,"frame.tsv"),"changed")
            @test_throws ErrorException X.analyze(run,m,cutoff,wfc)
            @test !ispath(joinpath(dir,"analysis"))
            write(joinpath(dir,"frame.tsv"),original_frame)
            @test X.analyze(run,m,cutoff,wfc)===nothing
            @test_throws ErrorException X.analyze(run,m,cutoff,wfc)
            rows=D.V.tsv(joinpath(dir,"analysis/planes.tsv"))
            for (i,r) in enumerate(rows)
                f=geo.cell\(pts[i,:]/D.G.C.BOHR_NM)
                truth=2abs2(coeff[1,2]+sum(coeff[g,2]*cis(2pi*dot(miller[:,g],f))+
                    conj(coeff[g,2])*cis(-2pi*dot(miller[:,g],f)) for g in 2:3))/header.volume
                @test isapprox(num(r,"density"),truth;rtol=1e-13,atol=1e-18)
            end
        end
        X.summarize(run); verify_saved(run)
        @test_throws ErrorException X.summarize(run)
    end
end

if !isempty(ARGS)
    length(ARGS)==1 || error("Expected one saved run")
    @testset "Independent saved real-orbital tables" begin
        verify_saved(only(ARGS))
    end
end
