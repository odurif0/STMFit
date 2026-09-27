using Test, TOML, LinearAlgebra, Statistics
include(joinpath(@__DIR__,"qe_diagonal_precision.jl"))
const P=QEDiagonalPrecision
const D=P.D; const S=P.S; const W=P.W
const ROOT=normpath(joinpath(@__DIR__,".."))
const CONFIG=joinpath(ROOT,"config/qe_diagonal_precision.toml")
const CG_CONFIG=joinpath(ROOT,"config/qe_cg_precision.toml")
num(r,k)=parse(Float64,r[k])
BLAS.set_num_threads(1)

function verify_saved(run,m; stages=P.STAGES, comparison=true)
    @test stages in (P.STAGES,("reference",),("reference","lo"))
    @test !comparison || stages==P.STAGES
    s,meta=P.verify_inputs(run,m); hs=s["preprocessing"]["diagnostic_heights_nm"]
    per=(round(Int,2s["preprocessing"]["half_nm"]/s["preprocessing"]["step_nm"])+1)^2
    old=P.reference_state(joinpath(run,m,"reference/data-file-schema.xml"),s)
    planes=Dict(); reports=Dict()
    for stage in stages
        dir=joinpath(run,m,stage); out=joinpath(dir,"analysis")
        r=TOML.parsefile(joinpath(out,"summary.toml")); reports[stage]=r
        @test all(values(r["checks"])) && !r["new_scf"] && !r["mold_adopted"] && !r["new_cube_validation"]
        @test r["stage"]==stage && r["molecule"]==m
        @test r["source_scf_error_ry"]==old.sp.scf_error_ry<=5e-5
        @test r["config_sha256"]==S.sha(joinpath(run,"settings.toml"))
        @test r["xml_sha256"]==S.sha(joinpath(dir,"data-file-schema.xml"))
        @test r["source_density_sha256"]==meta["source_sha256"]["charge-density.dat"]
        @test r["source_paw_sha256"]==meta["source_sha256"]["paw.txt"]
        expected=stage=="reference" ? meta["source_sha256"]["wfc1.dat"] : P.X.manifest_hash(
            joinpath(dir,"converged_before_analysis.sha256"),"work/$(m)_central.save/wfc1.dat")
        @test r["wfc_sha256"]==expected
        @test r["diago_thr_init_ry"]==(stage=="reference" ? 0. : P.threshold(s,stage))
        a=stage=="reference" ? old : P.nscf_state(joinpath(dir,"data-file-schema.xml"),s,old,stage)
        if P.iscg(s) || haskey(r,"diagonalization")
            @test r["diagonalization"]==S.tag(S.tag(a.inp,"electron_control"),"diagonalization")
        end
        @test r["state_kind"]==(stage=="reference" ? "scf_reference" : "nscf")
        @test a.sp.accepted==(stage=="reference")
        @test r["native_grid"]==a.geo.dims && r["native_grid"]==old.geo.dims
        @test r["selected_bands"]==[b.band for b in a.sp.rows if b.ildos_weight>0]
        @test r["selected_weights"]==[b.ildos_weight for b in a.sp.rows if b.ildos_weight>0]
        @test r["fermi_ev"]==a.sp.fermi_ev && r["emin_ev"]==a.sp.emin_ev && r["emax_ev"]==a.sp.emax_ev
        @test r["parseval_smooth_norm"]>0 && isfinite(r["parseval_smooth_norm"])
        @test min(r["point_clearance_nm"],r["vertex_clearance_nm"])>0
        bands=D.V.tsv(joinpath(out,"bands.tsv"))
        @test length(bands)==length(a.sp.rows)
        @test all(parse(Int,b["band"])==v.band && num(b,"energy_ev")==v.energy_ev &&
            num(b,"ildos_weight")==v.ildos_weight for (b,v) in zip(bands,a.sp.rows))
        frame=D.G.C.read_frame(joinpath(dir,"frame.tsv"))
        points=reduce(vcat,[S.plane_points(frame,h,s) for h in hs]); ns=D.native_vertices(a.geo.cell,a.geo.dims,points)
        fractions=hcat(ns.vertices,ns.fractions); nv=length(ns.indices)
        qs=D.V.tsv(joinpath(out,"queries.tsv")); rows=D.V.tsv(joinpath(out,"planes.tsv")); planes[stage]=rows
        @test length(rows)==r["plane_points"]==per*length(hs)
        @test nv==r["native_vertices"] && length(qs)==nv+length(rows)
        for (i,q) in enumerate(qs)
            @test parse(Int,q["query"])==i && q["kind"]==(i<=nv ? "native_vertex" : "plane")
            @test isfinite(num(q,"density")) && num(q,"density")>=0 && num(q,"density")==num(q,"repeat")
            @test [num(q,k) for k in ("fractional_x","fractional_y","fractional_z")]==collect(fractions[:,i])
        end
        for (i,q) in enumerate(rows)
            @test num(q,"height_nm")==hs[1+div(i-1,per)] && parse(Int,q["pixel"])==1+mod(i-1,per)
            @test [num(q,k) for k in ("x_nm","y_nm","z_nm")]==collect(points[i,:])
            @test num(q,"density")==num(qs[nv+i],"density")
        end
        if stage=="reference"
            ref=D.V.tsv(joinpath(run,m,"reference_queries.tsv")); refp=D.V.tsv(joinpath(run,m,"reference_planes.tsv"))
            @test length(ref)==length(qs) && length(refp)==length(rows)
            @test all(num(u,"density")==num(v,"density") for (u,v) in zip(ref,qs))
            @test all(num(u,"density")==num(v,"density") for (u,v) in zip(refp,rows))
        end
        stats=D.V.tsv(joinpath(out,"plane_summary.tsv")); tails=D.V.tsv(joinpath(out,"decay.tsv"))
        @test length(stats)==length(hs) && length(tails)==per*(length(hs)-1)
        for (j,h) in enumerate(hs)
            vs=num.(rows[(j-1)*per+1:j*per],"density"); q=stats[j]
            @test num(q,"height_nm")==h && parse(Int,q["points"])==per && parse(Int,q["zeros"])==count(iszero,vs)
            @test num(q,"minimum")==minimum(vs) && num(q,"median")==median(vs) && num(q,"maximum")==maximum(vs)
            @test isapprox(num(q,"density_sum"),sum(vs);rtol=16eps(),atol=0)
        end
        for (i,q) in enumerate(tails)
            j=1+div(i-1,per); pixel=1+mod(i-1,per)
            x=num(rows[(j-1)*per+pixel],"density"); y=num(rows[j*per+pixel],"density")
            @test num(q,"height_from_nm")==hs[j] && num(q,"height_to_nm")==hs[j+1] && parse(Int,q["pixel"])==pixel
            @test num(q,"lower_density")==x && num(q,"upper_density")==y && parse(Bool,q["rising"])==(y>x)
            @test isequal(num(q,"ratio"),x>0 ? y/x : NaN)
            expected=x>0 && y>0 ? log(x/y)/(hs[j+1]-hs[j]) : NaN
            @test isequal(num(q,"log_decay_per_nm"),expected) || isapprox(num(q,"log_decay_per_nm"),expected;rtol=16eps(),atol=0)
        end
    end
    for stage in intersect(P.TRIALS,stages), key in ("plane_waves","miller_sha256","source_density_sha256","source_paw_sha256")
        @test reports[stage][key]==reports["reference"][key]
    end
    comparison || return nothing
    pairs=D.V.tsv(joinpath(run,m,"comparison/paired_planes.tsv")); stats=D.V.tsv(joinpath(run,m,"comparison/summary.tsv"))
    @test length(pairs)==3per*length(hs) && length(stats)==3length(hs)
    for (from,to) in (("reference","lo"),("reference","hi"),("lo","hi"))
        rows=filter(r->r["from"]==from && r["to"]==to,pairs)
        @test length(rows)==per*length(hs)
        for (i,q) in enumerate(rows)
            x=num(planes[from][i],"density"); y=num(planes[to][i],"density")
            @test q["molecule"]==m && q["height_nm"]==planes[from][i]["height_nm"] && q["pixel"]==planes[from][i]["pixel"]
            @test all(planes[from][i][k]==planes[to][i][k] for k in ("x_nm","y_nm","z_nm"))
            @test num(q,"density_from")==x && num(q,"density_to")==y && num(q,"difference")==y-x
            @test isequal(num(q,"ratio"),x>0 ? y/x : NaN)
        end
        for (j,h) in enumerate(hs)
            q=only(filter(r->r["from"]==from && r["to"]==to && num(r,"height_nm")==h,stats))
            ix=(j-1)*per+1:j*per; x=num.(planes[from][ix],"density"); y=num.(planes[to][ix],"density")
            ratios=[a>0 ? b/a : NaN for (a,b) in zip(x,y)]
            @test parse(Int,q["points"])==per && q["molecule"]==m
            @test isequal(num(q,"minimum_ratio"),minimum(ratios)) && isequal(num(q,"median_ratio"),median(ratios)) && isequal(num(q,"maximum_ratio"),maximum(ratios))
            expected=sum(abs2,x)>0 ? sqrt(sum(abs2,y-x)/sum(abs2,x)) : NaN
            @test isequal(num(q,"relative_l2"),expected) || isapprox(num(q,"relative_l2"),expected;rtol=16eps(),atol=0)
            expected=sum(x)>0 ? sum(y)/sum(x) : NaN
            @test isequal(num(q,"sum_ratio"),expected) || isapprox(num(q,"sum_ratio"),expected;rtol=16eps(),atol=0)
        end
    end
end

function verify_rejected(run,m)
    s,meta=P.verify_inputs(run,m)
    stages=Tuple(stage for stage in P.STAGES if isdir(joinpath(run,m,stage,"analysis")))
    @test stages in (("reference",),("reference","lo"))
    @test !ispath(joinpath(run,m,"comparison"))
    verify_saved(run,m;stages,comparison=false)
    trial=length(stages)==1 ? "lo" : "hi"
    dir=joinpath(run,m,trial); work=P.workdir(run,m,trial)
    @test !ispath(joinpath(dir,"analysis")) && !ispath(joinpath(dir,"data-file-schema.xml"))
    @test !ispath(joinpath(dir,"converged_before_analysis.sha256"))
    old=P.reference_state(joinpath(run,m,"reference/data-file-schema.xml"),s)
    a=P.nscf_state(joinpath(work,"data-file-schema.xml"),s,old,trial)
    @test !a.sp.accepted
    for f in ("charge-density.dat","paw.txt",keys(meta["pseudo_sha256"])...)
        @test S.sha(joinpath(work,f))==meta["source_sha256"][f]
    end
    log=read(joinpath(dir,"pw_nscf.out"),String)
    @test count("JOB DONE.",log)==1
    @test occursin(r"c_bands:\s+[1-9][0-9]* eigenvalues not converged",log)
    # Reproduce the actual rejection, not an unrelated missing-file exception.
    err=try
        P.check_nscf(run,m,trial)
        nothing
    catch e
        e
    end
    @test err isa ErrorException && err.msg=="Unfinished or unexpected NSCF"
    @test !ispath(joinpath(dir,"data-file-schema.xml"))
    @test !ispath(joinpath(run,m,"comparison"))
    println(m,": verified retained stages ",stages," and rejected ",trial,"; no paired result")
end

function record(io,items...)
    b=IOBuffer(); foreach(x->write(b,x),items); bytes=take!(b); n=Int32(length(bytes))
    write(io,n); write(io,bytes); write(io,n)
end
function toy_wfc(file,stage)
    miller=Int32[0 1 0;0 0 0;0 0 1]; coeff=zeros(ComplexF64,3,3); coeff[1,2]=.25
    stage=="lo" && (coeff[2,2]=.02+.01im; coeff[3,2]=.03)
    stage=="hi" && (coeff[2,2]=.021+.009im; coeff[3,2]=.029)
    open(file,"w") do io
        record(io,Int32(1),zeros(3),Int32(1),Int32(-1),1.)
        record(io,Int32.([3,3,1,3])); record(io,2pi*inv(20Matrix{Float64}(I,3,3))); record(io,miller)
        for col in eachcol(coeff); record(io,collect(col)); end
    end
    open(W.header,file)
end
function toy_xml(stage; cg=false)
    ref=stage=="reference"; thr=ref ? 0. : stage=="lo" ? 1e-10 : 1e-12
    """
    <qes Units="Hartree atomic units"><creator VERSION="7.4.1"/>
    <input><control_variables><calculation>$(ref ? "scf" : "nscf")</calculation><restart_mode>from_scratch</restart_mode><max_seconds>3000</max_seconds><forces>false</forces><stress>false</stress></control_variables>
    <atomic_species><species name="C"><mass>12.011</mass><pseudo_file>C.UPF</pseudo_file></species></atomic_species>
    <basis><gamma_only>true</gamma_only><ecutwfc>30</ecutwfc><ecutrho>360</ecutrho></basis>
    <electron_control><conv_thr>5e-8</conv_thr><max_nstep>100</max_nstep><mixing_mode>plain</mixing_mode><mixing_beta>0.3</mixing_beta><mixing_ndim>8</mixing_ndim><diagonalization>$(cg && !ref ? "cg" : "davidson")</diagonalization><diago_cg_maxiter>20</diago_cg_maxiter><diago_thr_init>$thr</diago_thr_init><diago_full_acc>$(ref ? "false" : "true")</diago_full_acc></electron_control>
    <dft><functional>PBE</functional><vdw_corr>grimme-d3</vdw_corr><dftd3_version>3</dftd3_version><dftd3_threebody>true</dftd3_threebody></dft>
    <bands><tot_charge>0</tot_charge><occupations>smearing</occupations></bands>
    <spin><lsda>false</lsda><noncolin>false</noncolin><spinorbit>false</spinorbit></spin>
    <k_points_IBZ><nk>1</nk><k_point weight="1">0 0 0</k_point></k_points_IBZ></input>
    <output><convergence_info><wf_collected>true</wf_collected><scf_conv><convergence_achieved>$(ref ? "true" : "false")</convergence_achieved><n_scf_steps>$(ref ? 39 : 1)</n_scf_steps><scf_error>$(ref ? "1e-10" : "0")</scf_error></scf_conv></convergence_info>
    <atomic_structure nat="1" alat="20"><atomic_positions><atom name="C" index="1">0.01 0.01 0.01</atom></atomic_positions>
    <cell><a1>20 0 0</a1><a2>0 20 0</a2><a3>0 0 20</a3></cell></atomic_structure>
    <basis_set><gamma_only>true</gamma_only><ecutwfc>30</ecutwfc><ecutrho>360</ecutrho><fft_grid nr1="40" nr2="40" nr3="40"></fft_grid></basis_set>
    <band_structure><nks>1</nks><nbnd>3</nbnd><nelec>2</nelec><lsda>false</lsda><noncolin>false</noncolin><spinorbit>false</spinorbit>
    <ks_energies><k_point weight="2">0 0 0</k_point><eigenvalues size="3">-0.1 0 0.1</eigenvalues></ks_energies>
    <smearing degauss="0.01">mv</smearing><fermi_energy>0.005</fermi_energy></band_structure><total_energy><etot>-1</etot></total_energy></output></qes>
    """
end

@testset "Fixed density scope and inherited real inputs" begin
    s=P.settings(CONFIG); @test P.main(["--help"])===nothing
    cg=P.settings(CG_CONFIG)
    @test !P.iscg(s) && P.iscg(cg)
    inherited=deepcopy(cg)
    delete!(inherited["model"],"diagonalization"); delete!(inherited["model"],"diago_cg_maxiter")
    @test inherited==s
    mktempdir() do dir
        for (section,key,value) in (("model","ecutwfc_ry",50.),("model","ecutrho_ry",360.),
                ("model","diago_thresholds_ry",[1e-9,1e-12]),("model","diago_full_acc",false),
                ("model","startingpot","atomic"),("model","startingwfc","atomic+random"),
                ("model","scf_acceptance_ry",1e-4),("model","nscf_max_seconds",7200),
                ("selection","adopt_reference",true),("preprocessing","normalize_density",true))
            bad=deepcopy(s); bad[section][key]=value
            path=joinpath(dir,"bad.toml"); open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException P.settings(path)
        end
        for (key,value) in (("diagonalization","david"),("diagonalization","ppcg"),("diago_cg_maxiter",100))
            bad=deepcopy(cg); bad["model"][key]=value
            path=joinpath(dir,"bad.toml"); open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException P.settings(path)
        end
        for key in ("diagonalization","diago_cg_maxiter")
            bad=deepcopy(cg); delete!(bad["model"],key)
            path=joinpath(dir,"bad.toml"); open(io->TOML.print(io,bad),path,"w")
            @test_throws ErrorException P.settings(path)
        end
        source=joinpath(dir,"source.xml")
        write(source,replace(toy_xml("reference"),"<diago_cg_maxiter>20<"=>"<diago_cg_maxiter>100<"))
        @test_throws ErrorException P.reference_state(source,cg)
        for m in P.MOLECULES
            xml=joinpath(ROOT,s["model"]["baseline_run"],m,"60/data-file-schema.xml")
            a=P.reference_state(xml,s)
            for trial in P.TRIALS
                input=P.input_text(xml,m,s,trial)
                @test input==read(joinpath(ROOT,"qe/diagonal_precision_20260927",m,trial,"pw_nscf.in"),String)
                @test P.input_text(xml,m,cg,trial)==replace(input," diagonalization='david'"=>" diagonalization='cg'\n diago_cg_maxiter=20")
                @test occursin("calculation='nscf'",input) && !occursin("calculation='scf'",input)
                @test occursin("ecutwfc=60.0",input) && occursin("ecutrho=720.0",input)
                @test occursin("startingpot='file'",input) && occursin("startingwfc='file'",input)
                @test occursin("restart_mode='from_scratch'",input) && occursin("max_seconds=3000",input)
                @test occursin("diago_thr_init=$(P.threshold(s,trial))",input) && occursin("diago_full_acc=.true.",input)
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
    end
    job=read(joinpath(ROOT,"hpc/qe_diagonal_precision.sbatch"),String)
    @test occursin("#SBATCH --ntasks-per-node=8",job) && occursin("#SBATCH --mem=96000MB",job)
    @test occursin("#SBATCH --time=02:00:00",job) && occursin("#SBATCH --no-requeue",job)
    @test occursin("for STMFIT_DIAG_TRIAL in lo hi",job) && count("srun -n 8",job)==1
    @test !occursin("pp.x",job) && !occursin("sbatch ",job) && !occursin("pw_scf.in",job)
    run(`bash -n $(joinpath(ROOT,"hpc/qe_diagonal_precision.sbatch"))`)
end

@testset "Synthetic fixed density: $(basename(config))" for config in (CONFIG,CG_CONFIG)
    mktempdir() do run
        cp(config,joinpath(run,"settings.toml")); s=P.settings(config); mkpath(joinpath(run,"pseudo"))
        write(joinpath(run,"pseudo/C.UPF"),"<PP_HEADER element=\"C\" is_paw=\"true\" mesh_size=\"3\"/><PP_R>0.01 0.02 0.1</PP_R><PP_AUGMENTATION cutoff_r_index=\"3\"></PP_AUGMENTATION><PP_BETA.1 cutoff_radius_index=\"3\"></PP_BETA.1>")
        for m in P.MOLECULES
            dir=joinpath(run,m); source=joinpath(dir,"source"); mkpath(source)
            write(joinpath(source,"data-file-schema.xml"),toy_xml("reference"))
            h=toy_wfc(joinpath(source,"wfc1.dat"),"reference")
            write(joinpath(source,"charge-density.dat"),"frozen density"); write(joinpath(source,"paw.txt"),"frozen PAW")
            cp(joinpath(run,"pseudo/C.UPF"),joinpath(source,"C.UPF"))
            frame="key\tvalue\norigin_nm\t0.5,0.5,0.2\nt_axis\t1,0,0\nu_axis\t0,1,0\n"
            for stage in P.STAGES
                mkpath(joinpath(dir,stage)); write(joinpath(dir,stage,"frame.tsv"),frame)
            end
            cp(joinpath(source,"data-file-schema.xml"),joinpath(dir,"reference/data-file-schema.xml"))
            for trial in P.TRIALS
                write(joinpath(dir,trial,"pw_nscf.in"),P.input_text(joinpath(source,"data-file-schema.xml"),m,s,trial))
            end
            geo=D.xml_geometry(joinpath(source,"data-file-schema.xml")); f=D.G.C.read_frame(joinpath(dir,"reference/frame.tsv"))
            hs=s["preprocessing"]["diagnostic_heights_nm"]; pts=reduce(vcat,[S.plane_points(f,z,s) for z in hs]); ns=D.native_vertices(geo.cell,geo.dims,pts)
            fractions=hcat(ns.vertices,ns.fractions); rho=2*.25^2/h.volume
            S.table(joinpath(dir,"reference_queries.tsv"),[(;density=rho,fractional_x=fractions[1,i],fractional_y=fractions[2,i],fractional_z=fractions[3,i]) for i in axes(fractions,2)])
            S.table(joinpath(dir,"reference_planes.tsv"),[(;height_nm=z,pixel,density=rho) for z in hs for pixel in 1:289])
            files=("reference/data-file-schema.xml","reference/frame.tsv","lo/frame.tsv","hi/frame.tsv","lo/pw_nscf.in","hi/pw_nscf.in","reference_queries.tsv","reference_planes.tsv")
            meta=Dict("molecule"=>m,"config_sha256"=>S.sha(config),"source_sha256"=>Dict(f=>S.sha(joinpath(source,f)) for f in readdir(source)),
                "pseudo_sha256"=>Dict("C.UPF"=>S.sha(joinpath(run,"pseudo/C.UPF"))),"input_sha256"=>Dict(f=>S.sha(joinpath(dir,f)) for f in files))
            open(io->TOML.print(io,meta),joinpath(dir,"metadata.toml"),"w")
            @test P.analyze(run,m,"reference",joinpath(source,"wfc1.dat"))===nothing
            verify_saved(run,m;stages=("reference",),comparison=false)
            for trial in P.TRIALS
                work=P.workdir(run,m,trial); mkpath(dirname(work)); cp(source,work)
                @test P.check_start(run,m,trial)===nothing
                write(joinpath(work,"charge-density.dat"),"changed")
                @test_throws ErrorException P.check_start(run,m,trial)
                write(joinpath(work,"charge-density.dat"),"frozen density")
                candidate=joinpath(work,"data-file-schema.xml"); xml=toy_xml(trial;cg=P.iscg(s))
                marker=P.iscg(s) ? "CG style diagonalization" : "Davidson diagonalization with overlap"
                log="The potential is recalculated from file\nStarting wfcs from file\nBand Structure Calculation\n$marker\nethr = $(P.threshold(s,trial))\nEnd of band structure calculation\nJOB DONE.\n"
                logfile=joinpath(dir,trial,"pw_nscf.out"); write(logfile,log)
                for (before,after) in (("<calculation>nscf","<calculation>scf"),("<convergence_achieved>false","<convergence_achieved>true"),
                        ("<wf_collected>true","<wf_collected>false"),("<functional>PBE<","<functional>PBE0<"),
                        ("<ecutwfc>30<","<ecutwfc>31<"),("<diago_full_acc>true","<diago_full_acc>false"),
                        ("<diago_thr_init>$(P.threshold(s,trial))<","<diago_thr_init>$(P.threshold(s,trial)/2)<"),
                        ("<conv_thr>5e-8","<conv_thr>1e-3"),("0.01 0.01 0.01","0.01 0.01 0.02"))
                    @test occursin(before,xml); write(candidate,replace(xml,before=>after))
                    @test_throws ErrorException P.check_nscf(run,m,trial)
                    @test !ispath(joinpath(dir,trial,"data-file-schema.xml"))
                end
                for solver in (P.iscg(s) ? ("davidson","ppcg") : ("cg","ppcg"))
                    tag="<diagonalization>$(P.iscg(s) ? "cg" : "davidson")<"
                    @test occursin(tag,xml)
                    write(candidate,replace(xml,tag=>"<diagonalization>$solver<"))
                    @test_throws ErrorException P.check_nscf(run,m,trial)
                end
                if P.iscg(s)
                    for replacement in ("<diago_cg_maxiter>100</diago_cg_maxiter>","")
                        write(candidate,replace(xml,"<diago_cg_maxiter>20</diago_cg_maxiter>"=>replacement))
                        @test_throws ErrorException P.check_nscf(run,m,trial)
                    end
                    write(candidate,xml)
                    for replacement in ("","PPCG style diagonalization","Davidson diagonalization with overlap",marker*"\n"*marker)
                        write(logfile,replace(log,marker=>replacement))
                        @test_throws ErrorException P.check_nscf(run,m,trial)
                        @test !ispath(joinpath(dir,trial,"data-file-schema.xml"))
                    end
                end
                write(candidate,xml)
                @test_throws ErrorException P.reference_state(candidate,s)
                for bad in (replace(log,"JOB DONE."=>""),log*"c_bands: 1 eigenvalues not converged\n",replace(log,"Starting wfcs from file"=>"Starting wfcs are random"))
                    write(logfile,bad); @test_throws ErrorException P.check_nscf(run,m,trial)
                    occursin("eigenvalues not converged",bad) && verify_rejected(run,m)
                end
                write(logfile,log); write(joinpath(work,"paw.txt"),"changed")
                @test_throws ErrorException P.check_nscf(run,m,trial)
                write(joinpath(work,"paw.txt"),"frozen PAW")
                toy_wfc(joinpath(work,"wfc1.dat"),trial)
                @test P.check_nscf(run,m,trial)===nothing
                write(joinpath(dir,trial,"converged_before_analysis.sha256"),S.sha(joinpath(work,"wfc1.dat"))*"  work/$(m)_central.save/wfc1.dat\n")
                @test P.analyze(run,m,trial,joinpath(work,"wfc1.dat"))===nothing
                @test_throws ErrorException P.analyze(run,m,trial,joinpath(work,"wfc1.dat"))
                trial=="lo" && verify_saved(run,m;stages=("reference","lo"),comparison=false)
                w=W.read_wfc(joinpath(work,"wfc1.dat"),[2]); rows=D.V.tsv(joinpath(dir,trial,"analysis/planes.tsv"))
                for (i,r) in enumerate(rows)
                    v=geo.cell\(pts[i,:]/D.G.C.BOHR_NM)
                    truth=2abs2(w.coeff[1,1]+sum(w.coeff[g,1]*cis(2pi*dot(w.miller[:,g],v))+
                        conj(w.coeff[g,1])*cis(-2pi*dot(w.miller[:,g],v)) for g in 2:w.ng))/w.volume
                    @test isapprox(num(r,"density"),truth;rtol=1e-13,atol=1e-18)
                end
            end
            @test P.summarize(run,m)===nothing
            @test_throws ErrorException P.summarize(run,m)
            verify_saved(run,m)
        end
    end
end

if !isempty(ARGS)
    if length(ARGS)==3 && ARGS[1]=="--incomplete"
        @testset "Incomplete real result: retained tables and strict rejection" begin
            verify_rejected(ARGS[2:end]...)
        end
    else
        length(ARGS)==2 || error("Expected [--incomplete] SAVED_RUN MOLECULE")
        @testset "Independent saved fixed-density precision tables" begin
            verify_saved(ARGS...)
        end
    end
end
