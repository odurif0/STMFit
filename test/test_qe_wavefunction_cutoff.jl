using Test, TOML, LinearAlgebra, Statistics
include(joinpath(@__DIR__,"qe_wavefunction_cutoff.jl"))
const B=QEWavefunctionCutoff
const D=B.D; const S=B.S; const W=B.W
const ROOT=normpath(joinpath(@__DIR__,".."))
const CONFIG=joinpath(ROOT,"config/qe_wavefunction_cutoff.toml")
num(r,k)=parse(Float64,r[k])
BLAS.set_num_threads(1)

function verify_saved(run,m)
    s,meta=B.verify_inputs(run,m); p=s["preprocessing"]; hs=p["diagnostic_heights_nm"]
    per=(round(Int,2p["half_nm"]/p["step_nm"])+1)^2
    planes=Dict(); reports=Dict()
    for cutoff in (50,60)
        dir=joinpath(run,m,string(cutoff)); out=joinpath(dir,"analysis")
        report=TOML.parsefile(joinpath(out,"summary.toml")); reports[cutoff]=report
        @test all(values(report["checks"])) && !report["mold_adopted"] && !report["new_cube_validation"]
        @test report["ecutwfc_ry"]==cutoff && report["ecutrho_ry"]==720
        @test report["config_sha256"]==S.sha(joinpath(run,"settings.toml"))
        @test report["xml_sha256"]==S.sha(joinpath(dir,"data-file-schema.xml"))
        expected=cutoff==50 ? meta["baseline_wfc_sha256"] : B.X.manifest_hash(
            joinpath(run,m,"converged_before_analysis.sha256"),"work/$(m)_central.save/wfc1.dat")
        @test report["wfc_sha256"]==expected
        @test 0<=report["scf_error_ry"]<=s["model"]["scf_acceptance_ry"]
        @test min(report["point_clearance_nm"],report["vertex_clearance_nm"])>0
        @test report["parseval_smooth_norm"]>0 && isfinite(report["parseval_smooth_norm"])
        a=B.state(joinpath(dir,"data-file-schema.xml"),s,cutoff)
        @test report["native_grid"]==a.geo.dims
        @test report["selected_bands"]==[r.band for r in a.sp.rows if r.ildos_weight>0]
        @test report["selected_weights"]==[r.ildos_weight for r in a.sp.rows if r.ildos_weight>0]
        bands=D.V.tsv(joinpath(out,"bands.tsv"))
        @test length(bands)==length(a.sp.rows)
        @test all(parse(Int,r["band"])==v.band && num(r,"energy_ev")==v.energy_ev &&
            num(r,"ildos_weight")==v.ildos_weight for (r,v) in zip(bands,a.sp.rows))
        frame=D.G.C.read_frame(joinpath(dir,"frame.tsv"))
        points=reduce(vcat,[S.plane_points(frame,h,s) for h in hs]); ns=D.native_vertices(a.geo.cell,a.geo.dims,points)
        fractions=hcat(ns.vertices,ns.fractions); nv=length(ns.indices)
        q=D.V.tsv(joinpath(out,"queries.tsv")); rows=D.V.tsv(joinpath(out,"planes.tsv")); planes[cutoff]=rows
        @test length(rows)==report["plane_points"]==length(hs)*per
        @test nv==report["native_vertices"] && length(q)==nv+length(rows)
        for (i,r) in enumerate(q)
            @test parse(Int,r["query"])==i && r["kind"]==(i<=nv ? "native_vertex" : "plane")
            @test isfinite(num(r,"density")) && num(r,"density")>=0 && num(r,"density")==num(r,"repeat")
            @test [num(r,k) for k in ("fractional_x","fractional_y","fractional_z")]==collect(fractions[:,i])
        end
        for (i,r) in enumerate(rows)
            @test num(r,"height_nm")==hs[1+div(i-1,per)] && parse(Int,r["pixel"])==1+mod(i-1,per)
            @test [num(r,k) for k in ("x_nm","y_nm","z_nm")]==collect(points[i,:])
            @test num(r,"density")==num(q[nv+i],"density")
        end
        if cutoff==50
            ref=D.V.tsv(joinpath(run,m,"reference_queries.tsv"))
            old=D.V.tsv(joinpath(run,m,"reference_planes.tsv"))
            @test length(ref)==length(q) && length(old)==length(rows)
            @test all(num(u,"density")==num(v,"density") for (u,v) in zip(ref,q))
            @test all(num(u,"density")==num(v,"density") for (u,v) in zip(old,rows))
        end
        stats=D.V.tsv(joinpath(out,"plane_summary.tsv")); tails=D.V.tsv(joinpath(out,"decay.tsv"))
        @test length(stats)==length(hs) && length(tails)==per*(length(hs)-1)
        for (j,h) in enumerate(hs)
            vs=num.(rows[(j-1)*per+1:j*per],"density"); r=stats[j]
            @test num(r,"height_nm")==h && parse(Int,r["points"])==per && parse(Int,r["zeros"])==count(iszero,vs)
            @test num(r,"minimum")==minimum(vs) && num(r,"median")==median(vs) && num(r,"maximum")==maximum(vs)
            @test isapprox(num(r,"density_sum"),sum(vs);rtol=16eps(),atol=0)
        end
        for (i,r) in enumerate(tails)
            j=1+div(i-1,per); pixel=1+mod(i-1,per)
            x=num(rows[(j-1)*per+pixel],"density"); y=num(rows[j*per+pixel],"density")
            @test num(r,"height_from_nm")==hs[j] && num(r,"height_to_nm")==hs[j+1] && parse(Int,r["pixel"])==pixel
            @test num(r,"lower_density")==x && num(r,"upper_density")==y && parse(Bool,r["rising"])==(y>x)
            @test isequal(num(r,"ratio"),x>0 ? y/x : NaN)
            expected=x>0 && y>0 ? log(x/y)/(hs[j+1]-hs[j]) : NaN
            @test isequal(num(r,"log_decay_per_nm"),expected) || isapprox(num(r,"log_decay_per_nm"),expected;rtol=16eps(),atol=0)
        end
    end
    @test reports[60]["plane_waves"]>reports[50]["plane_waves"]
    pairs=D.V.tsv(joinpath(run,m,"comparison/paired_planes.tsv"))
    stats=D.V.tsv(joinpath(run,m,"comparison/summary.tsv"))
    @test length(pairs)==per*length(hs) && length(stats)==length(hs)
    for (i,r) in enumerate(pairs)
        x=num(planes[50][i],"density"); y=num(planes[60][i],"density")
        @test r["molecule"]==m && r["height_nm"]==planes[50][i]["height_nm"] && r["pixel"]==planes[50][i]["pixel"]
        @test num(r,"density_50")==x && num(r,"density_60")==y && num(r,"difference")==y-x
        @test isequal(num(r,"ratio"),x>0 ? y/x : NaN)
    end
    for (j,h) in enumerate(hs)
        ix=(j-1)*per+1:j*per; x=num.(planes[50][ix],"density"); y=num.(planes[60][ix],"density")
        ratio=[a>0 ? b/a : NaN for (a,b) in zip(x,y)]; r=stats[j]
        @test r["molecule"]==m && num(r,"height_nm")==h && parse(Int,r["points"])==per
        @test isequal(num(r,"minimum_ratio"),minimum(ratio)) && isequal(num(r,"median_ratio"),median(ratio)) && isequal(num(r,"maximum_ratio"),maximum(ratio))
        expected=sum(abs2,x)>0 ? sqrt(sum(abs2,y-x)/sum(abs2,x)) : NaN
        @test isequal(num(r,"relative_l2"),expected) || isapprox(num(r,"relative_l2"),expected;rtol=16eps(),atol=0)
        expected=sum(x)>0 ? sum(y)/sum(x) : NaN
        @test isequal(num(r,"sum_ratio"),expected) || isapprox(num(r,"sum_ratio"),expected;rtol=16eps(),atol=0)
    end
end

function record(io,items...)
    b=IOBuffer(); foreach(x->write(b,x),items); bytes=take!(b); n=Int32(length(bytes))
    write(io,n); write(io,bytes); write(io,n)
end
function toy_wfc(file,cutoff)
    miller=cutoff==50 ? Int32[0 1 0;0 0 0;0 0 1] : Int32[0 1 0 2;0 0 0 0;0 0 1 0]
    coeff=zeros(ComplexF64,size(miller,2),3); coeff[1,2]=.25
    cutoff==60 && (coeff[2,2]=.02+.01im; coeff[3,2]=.03; coeff[4,2]=-.015im)
    open(file,"w") do io
        record(io,Int32(1),zeros(3),Int32(1),Int32(-1),1.)
        record(io,Int32.([size(miller,2),size(miller,2),1,3])); record(io,2pi*inv(20Matrix{Float64}(I,3,3))); record(io,miller)
        for col in eachcol(coeff); record(io,collect(col)); end
    end
    open(W.header,file),miller,coeff
end
function toy_xml(cutoff)
    """
    <qes Units="Hartree atomic units"><creator VERSION="7.4.1"/>
    <input><atomic_species><species name="C"><mass>12.011</mass><pseudo_file>C.UPF</pseudo_file></species></atomic_species>
    <basis><gamma_only>true</gamma_only><ecutwfc>$(cutoff/2)</ecutwfc><ecutrho>360</ecutrho></basis>
    <electron_control><conv_thr>5e-8</conv_thr><max_nstep>100</max_nstep><mixing_mode>plain</mixing_mode><mixing_beta>0.3</mixing_beta><mixing_ndim>8</mixing_ndim><diagonalization>davidson</diagonalization></electron_control>
    <dft><functional>PBE</functional><vdw_corr>grimme-d3</vdw_corr><dftd3_version>3</dftd3_version><dftd3_threebody>true</dftd3_threebody></dft>
    <bands><tot_charge>0</tot_charge><occupations>smearing</occupations></bands>
    <spin><lsda>false</lsda><noncolin>false</noncolin><spinorbit>false</spinorbit></spin></input>
    <output><convergence_info><wf_collected>true</wf_collected><scf_conv><convergence_achieved>true</convergence_achieved><scf_error>1e-10</scf_error></scf_conv></convergence_info>
    <atomic_structure nat="1" alat="20"><atomic_positions><atom name="C" index="1">0.01 0.01 0.01</atom></atomic_positions>
    <cell><a1>20 0 0</a1><a2>0 20 0</a2><a3>0 0 20</a3></cell></atomic_structure>
    <basis_set><gamma_only>true</gamma_only><ecutwfc>$(cutoff/2)</ecutwfc><ecutrho>360</ecutrho><fft_grid nr1="40" nr2="40" nr3="40"></fft_grid></basis_set>
    <band_structure><nks>1</nks><nbnd>3</nbnd><nelec>2</nelec><lsda>false</lsda><noncolin>false</noncolin><spinorbit>false</spinorbit>
    <ks_energies><k_point weight="2">0 0 0</k_point><eigenvalues size="3">-0.1 0 0.1</eigenvalues></ks_energies>
    <smearing degauss="0.01">mv</smearing><fermi_energy>0.005</fermi_energy></band_structure><total_energy><etot>-1</etot></total_energy></output></qes>
    """
end

@testset "Bounded scope, input generation and job limits" begin
    s=B.settings(CONFIG)
    @test B.main(["--help"])===nothing
    mktempdir() do dir
        for (section,key,value) in (("model","ecutwfc_ry",[50.,70.]),("model","ecutrho_ry",800.),
                ("model","scf_acceptance_ry",1e-4),("model","startingwfc","file"),("model","scf_max_seconds",14000),
                ("selection","adopt_reference",true),("preprocessing","normalize_density",true),("preprocessing","half_nm",.64))
            bad=deepcopy(s); bad[section][key]=value
            file=joinpath(dir,"bad.toml"); open(io->TOML.print(io,bad),file,"w")
            @test_throws ErrorException B.settings(file)
        end
        for m in B.MOLECULES
            xml=joinpath(ROOT,s["model"]["baseline_run"],m,"720/data-file-schema.xml")
            a=B.state(xml,s,50); input=B.input_text(xml,m,s)
            @test occursin("ecutwfc=60.0",input) && occursin("ecutrho=720.0",input)
            @test occursin("startingwfc='atomic+random'",input) && occursin("startingpot='atomic'",input)
            @test occursin("restart_mode='from_scratch'",input) && occursin("max_seconds=12600",input)
            @test occursin("conv_thr="*(m=="glcn" ? "1.0e-7" : "5.0e-5"),input)
            lines=split(input,'\n'); ia=findfirst(==("ATOMIC_POSITIONS bohr"),lines); ic=findfirst(==("CELL_PARAMETERS bohr"),lines)
            @test findfirst(==("K_POINTS gamma"),lines)-ia-1==length(a.geo.geometry.atoms)
            for k in 1:3; @test parse.(Float64,split(lines[ic+k]))==a.geo.cell[:,k]; end
            for (k,atom) in enumerate(a.geo.geometry.atoms)
                cols=split(lines[ia+k]); @test D.P.Z[cols[1]]==atom.z
                @test parse.(Float64,cols[2:4])*D.G.C.BOHR_NM==atom.position_nm
            end
        end
    end
    job=read(joinpath(ROOT,"hpc/qe_wavefunction_cutoff.sbatch"),String)
    @test occursin("#SBATCH --ntasks-per-node=8",job) && occursin("#SBATCH --mem=96000MB",job)
    @test occursin("#SBATCH --time=04:00:00",job) && occursin("#SBATCH --no-requeue",job)
    @test count("srun -n 8",job)==1 && !occursin("pp.x",job) && !occursin("sbatch ",job)
    run(`bash -n $(joinpath(ROOT,"hpc/qe_wavefunction_cutoff.sbatch"))`)
end

@testset "Synthetic expanded basis, SCF rejection and complete paired outputs" begin
    mktempdir() do run
        cp(CONFIG,joinpath(run,"settings.toml")); s=B.settings(CONFIG); mkpath(joinpath(run,"pseudo"))
        write(joinpath(run,"pseudo/C.UPF"),"<PP_HEADER element=\"C\" is_paw=\"true\" mesh_size=\"3\"/><PP_R>0.01 0.02 0.1</PP_R><PP_AUGMENTATION cutoff_r_index=\"3\"></PP_AUGMENTATION><PP_BETA.1 cutoff_radius_index=\"3\"></PP_BETA.1>")
        for m in B.MOLECULES
            dir=joinpath(run,m); mkpath(joinpath(dir,"50")); mkpath(joinpath(dir,"60"))
            work=joinpath(dir,"work",m*"_central.save"); mkpath(work)
            write(joinpath(dir,"50/data-file-schema.xml"),toy_xml(50))
            candidate=joinpath(work,"data-file-schema.xml"); write(candidate,toy_xml(60))
            frame="key\tvalue\norigin_nm\t0.5,0.5,0.2\nt_axis\t1,0,0\nu_axis\t0,1,0\n"
            for c in (50,60); write(joinpath(dir,string(c),"frame.tsv"),frame); end
            h,_,_=toy_wfc(joinpath(dir,"baseline.dat"),50)
            toy_wfc(joinpath(work,"wfc1.dat"),60)
            geo=D.xml_geometry(joinpath(dir,"50/data-file-schema.xml")); f=D.G.C.read_frame(joinpath(dir,"50/frame.tsv"))
            hs=s["preprocessing"]["diagnostic_heights_nm"]; pts=reduce(vcat,[S.plane_points(f,z,s) for z in hs]); ns=D.native_vertices(geo.cell,geo.dims,pts)
            fractions=hcat(ns.vertices,ns.fractions); rho=2*.25^2/h.volume
            S.table(joinpath(dir,"reference_queries.tsv"),[(;density=rho,fractional_x=fractions[1,i],fractional_y=fractions[2,i],fractional_z=fractions[3,i]) for i in axes(fractions,2)])
            S.table(joinpath(dir,"reference_planes.tsv"),[(;height_nm=z,pixel,density=rho) for z in hs for pixel in 1:289])
            files=("50/data-file-schema.xml","50/frame.tsv","60/frame.tsv","reference_queries.tsv","reference_planes.tsv")
            meta=Dict("molecule"=>m,"config_sha256"=>S.sha(CONFIG),"baseline_wfc_sha256"=>S.sha(joinpath(dir,"baseline.dat")),
                "pseudo_sha256"=>Dict("C.UPF"=>S.sha(joinpath(run,"pseudo/C.UPF"))),"input_sha256"=>Dict(x=>S.sha(joinpath(dir,x)) for x in files))
            open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
            write(joinpath(dir,"pw_scf.out"),"Initial potential from superposition of free atoms\nStarting wfcs are 3 randomized atomic wfcs\nconvergence has been achieved\nJOB DONE.\n")
            for (before,after) in (("<convergence_achieved>true","<convergence_achieved>false"),
                    ("<scf_error>1e-10","<scf_error>1e-3"),("<functional>PBE<","<functional>PBE0<"),
                    ("<ecutwfc>30.0<","<ecutwfc>25.0<"),("<mixing_ndim>8<","<mixing_ndim>20<"),
                    ("0.01 0.01 0.01","0.01 0.01 0.02"))
                @test occursin(before,toy_xml(60))
                write(candidate,replace(toy_xml(60),before=>after))
                @test_throws ErrorException B.check_scf(run,m)
                @test !ispath(joinpath(dir,"60/data-file-schema.xml"))
            end
            write(candidate,toy_xml(60)); @test B.check_scf(run,m)===nothing
            write(joinpath(dir,"converged_before_analysis.sha256"),S.sha(joinpath(work,"wfc1.dat"))*"  work/$(m)_central.save/wfc1.dat\n")
            write(joinpath(dir,"60/frame.tsv"),"tampered")
            @test_throws ErrorException B.analyze(run,m,60,joinpath(work,"wfc1.dat"))
            @test !ispath(joinpath(dir,"60/analysis"))
            write(joinpath(dir,"60/frame.tsv"),frame)
            for c in (50,60)
                file=c==50 ? joinpath(dir,"baseline.dat") : joinpath(work,"wfc1.dat")
                @test B.analyze(run,m,c,file)===nothing
                @test_throws ErrorException B.analyze(run,m,c,file)
                w=W.read_wfc(file,[2]); rows=D.V.tsv(joinpath(dir,string(c),"analysis/planes.tsv"))
                for (i,r) in enumerate(rows)
                    v=geo.cell\(pts[i,:]/D.G.C.BOHR_NM)
                    truth=2abs2(w.coeff[1,1]+sum(w.coeff[g,1]*cis(2pi*dot(w.miller[:,g],v))+
                        conj(w.coeff[g,1])*cis(-2pi*dot(w.miller[:,g],v)) for g in 2:w.ng))/w.volume
                    @test isapprox(num(r,"density"),truth;rtol=1e-13,atol=1e-18)
                end
            end
            @test B.summarize(run,m)===nothing
            @test_throws ErrorException B.summarize(run,m)
            verify_saved(run,m)
        end
    end
end

if !isempty(ARGS)
    length(ARGS)==2 || error("Expected SAVED_RUN MOLECULE")
    @testset "Independent saved wavefunction-cutoff tables" begin
        verify_saved(ARGS...)
    end
end
