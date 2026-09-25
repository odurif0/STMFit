#!/usr/bin/env julia
# Geometry-only preparation, followed by conditional analysis of saved roots.
module DiagnoseMolecularExterior
using TOML, Statistics
include(joinpath(@__DIR__,"verify_qe_substrate_reference.jl"))
include(joinpath(@__DIR__,"lib/molecular_exterior.jl"))
const V=VerifyQESubstrateReference
const R=V.R
const E=MolecularExterior
const n=V.n
const i=V.i
const b=V.b
const tsv=R.V.tsv

function settings(file)
    s=TOML.parsefile(file)
    Set(keys(s))==Set(["model","selection","preprocessing"]) || error("Unexpected sections")
    s["selection"]==Dict("policy"=>"report_all_original_domains_isovalues_and_complement",
        "footprint"=>"convex_hull_molecular_centers_plus_max_active_molecular_paw_radius",
        "support"=>"strict_common_exterior_both_periodic_footprints",
        "groups"=>["all_prior","common_exterior","removed_by_envelope"],
        "domains"=>["full_cu_gap","cu_half_gap","saved_molecular_half"],"adopt_reference"=>false) || error("Changed selection")
    s["preprocessing"]==Dict("coordinates"=>"unchanged_accepted_nm_float64",
        "predicates"=>"exact_rational_of_float64","require_offset_bbox_inside_primary_cell"=>true,
        "periodic_images"=>[-1,0,1],"additional_clearance_nm"=>0.,"read_density"=>false) || error("Changed geometry policy")
    s
end

function sources(root,s)
    m=s["model"]; run=joinpath(root,m["run"]); report=joinpath(root,m["report"])
    files=Dict("settings_sha256"=>joinpath(run,"settings.toml"),"grid_sha256"=>joinpath(run,"grid.tsv"),
        "summary_sha256"=>joinpath(report,"summary.tsv"),"paired_sha256"=>joinpath(report,"paired.tsv"))
    for mol in R.MOLECULES; files[mol*"_roots_sha256"]=joinpath(run,mol,"diagnostic/roots.tsv"); end
    Set(keys(m))==union(Set(keys(files)),Set(["run","report"])) || error("Unexpected model fields")
    for (key,path) in files; R.S.sha(path)==m[key] || error("Changed source: $path"); end
    (;run,report)
end

function geometry(run)
    fs=Dict{String,Any}()
    for m in R.MOLECULES
        dir=joinpath(run,m); meta=TOML.parsefile(joinpath(dir,"metadata.toml"))
        R.S.sha(joinpath(dir,"data-file-schema.xml"))==meta["xml_sha256"] || error("Changed XML")
        radii,rs=R.D.radii(run,dir)
        Dict(r.pseudo=>r.sha256 for r in rs)==meta["pseudo_sha256"] || error("Changed active PAW data")
        fs[m]=E.footprint(R.D.xml_geometry(joinpath(dir,"data-file-schema.xml")).geometry,radii)
    end
    fs["glcn"].cell==fs["glcnac"].cell || error("Different lateral cells")
    grid=tsv(joinpath(run,"grid.tsv")); rows=NamedTuple[]
    i.(grid,"pixel")==collect(1:length(grid)) || error("Unordered grid")
    for r in grid
        xy=(n(r,"x_nm"),n(r,"y_nm")); a=E.clearance(xy,fs["glcn"]); c=E.clearance(xy,fs["glcnac"])
        prior=b(r,"eligible"); exterior=a.exterior && c.exterior
        exterior && !prior && error("Exterior is not a subset of old PAW support")
        push!(rows,(;pixel=i(r,"pixel"),x_nm=xy[1],y_nm=xy[2],prior_eligible=prior,
            glcn_exterior=a.exterior,glcnac_exterior=c.exterior,
            glcn_boundary=a.boundary,glcnac_boundary=c.boundary,
            glcn_margin_nm=a.margin_nm,glcnac_margin_nm=c.margin_nm,
            common_exterior=exterior,removed_by_envelope=prior && !exterior))
    end
    (;fs,rows)
end

function prepare(root,config,out)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    s=settings(config); src=sources(root,s); g=geometry(src.run)
    R.S.newdir(out); cp(config,joinpath(out,"settings.toml"))
    R.table(joinpath(out,"membership.tsv"),g.rows)
    vertices=[(;molecule=m,vertex=k,x_nm=Float64(p[1]),y_nm=Float64(p[2]),
        radius_nm=Float64(g.fs[m].radius),cell_x_nm=Float64(g.fs[m].cell[1]),cell_y_nm=Float64(g.fs[m].cell[2]))
        for m in R.MOLECULES for (k,p) in enumerate(g.fs[m].poly)]
    R.table(joinpath(out,"hulls.tsv"),vertices)
    meta=Dict("membership_sha256"=>R.S.sha(joinpath(out,"membership.tsv")),
        "hulls_sha256"=>R.S.sha(joinpath(out,"hulls.tsv")),"config_sha256"=>R.S.sha(config),
        "sites"=>length(g.rows),"prior"=>count(p->p.prior_eligible,g.rows),
        "exterior"=>count(p->p.common_exterior,g.rows),"removed"=>count(p->p.removed_by_envelope,g.rows),
        "density_read"=>false,"reference_adopted"=>false)
    open(io->TOML.print(io,meta),joinpath(out,"geometry.toml"),"w")
    println(meta)
end

function saved_geometry(root,prepared)
    s=settings(joinpath(prepared,"settings.toml")); src=sources(root,s)
    meta=TOML.parsefile(joinpath(prepared,"geometry.toml"))
    for (file,key) in (("settings.toml","config_sha256"),("membership.tsv","membership_sha256"),("hulls.tsv","hulls_sha256"))
        R.S.sha(joinpath(prepared,file))==meta[key] || error("Changed prepared $file")
    end
    rows=tsv(joinpath(prepared,"membership.tsv"))
    !meta["density_read"] && !meta["reference_adopted"] || error("Not geometry-only preparation")
    (;s,src,meta,rows)
end

function analyze(root,prepared,out)
    VERSION.major==1 && VERSION.minor==13 || error("Julia 1.13 required")
    g=saved_geometry(root,prepared); ids=Dict(
        group=>Set(i(r,"pixel") for r in g.rows if b(r,group=="all_prior" ? "prior_eligible" : group))
        for group in g.s["selection"]["groups"])
    isempty(intersect(ids["common_exterior"],ids["removed_by_envelope"])) || error("Overlapping groups")
    union(ids["common_exterior"],ids["removed_by_envelope"])==ids["all_prior"] || error("Missing complement")
    roots=Dict(m=>tsv(joinpath(g.src.run,m,"diagnostic/roots.tsv")) for m in R.MOLECULES)
    old=tsv(joinpath(g.src.report,"summary.tsv")); oldpaired=tsv(joinpath(g.src.report,"paired.tsv"))
    domains=g.s["selection"]["domains"]; isos=sort!(unique(i.(old,"interval")))
    for m in R.MOLECULES
        Set((i(r,"pixel"),r["domain"],i(r,"interval")) for r in roots[m])==
            Set((p,d,j) for p in ids["all_prior"] for d in domains for j in isos) || error("Missing root outcomes")
        length(roots[m])==length(ids["all_prior"])*length(domains)*length(isos) || error("Duplicate roots")
    end
    state_rows=NamedTuple[]; pair_rows=NamedTuple[]; pair_summaries=NamedTuple[]
    for group in g.s["selection"]["groups"], domain in domains, interval in isos
        selected=Dict(m=>Dict(i(r,"pixel")=>r for r in roots[m] if
            r["domain"]==domain && i(r,"interval")==interval && i(r,"pixel") in ids[group]) for m in R.MOLECULES)
        molecular=Dict(m=>n(only(filter(r->r["molecule"]==m && r["domain"]==domain && i(r,"interval")==interval,old)),
            "molecular_mean_nm") for m in R.MOLECULES)
        for m in R.MOLECULES
            rr=[selected[m][p] for p in sort!(collect(keys(selected[m])))]; valid=filter(r->b(r,"valid"),rr)
            push!(state_rows,merge((;group,molecule=m,domain,interval,eligible=length(rr),
                absent=count(r->i(r,"root_count")==0 && i(r,"unresolved")==0,rr),
                multiple=count(r->i(r,"root_count")>1,rr),unresolved=count(r->i(r,"unresolved")>0,rr),
                one_descending=count(r->i(r,"descending")==1,rr)),V.stats(n.(valid,"z_nm"))))
        end
        pairs=NamedTuple[]
        for pixel in sort!(collect(ids[group]))
            a=selected["glcn"][pixel]; c=selected["glcnac"][pixel]
            b(a,"valid") && b(c,"valid") || continue
            background=n(c,"z_nm")-n(a,"z_nm"); original=molecular["glcnac"]-molecular["glcn"]
            push!(pairs,(;group,domain,interval,pixel,background_contrast_nm=background,
                original_molecular_contrast_nm=original,referenced_contrast_nm=original-background))
        end
        append!(pair_rows,pairs)
        push!(pair_summaries,merge((;group,domain,interval,eligible=length(ids[group]),
            original_molecular_contrast_nm=molecular["glcnac"]-molecular["glcn"],
            negative=count(r->r.referenced_contrast_nm<0,pairs),positive=count(r->r.referenced_contrast_nm>0,pairs),
            zero=count(r->r.referenced_contrast_nm==0,pairs)),V.stats([r.referenced_contrast_nm for r in pairs])))
    end
    # The unfiltered arm must exactly replay every previously verified paired value.
    prior=filter(r->r.group=="all_prior",pair_rows)
    length(prior)==length(oldpaired) || error("Changed old paired support")
    original_by_key=Dict((r["domain"],i(r,"interval"),i(r,"pixel"))=>r for r in oldpaired)
    for r in prior
        oldrow=original_by_key[(r.domain,r.interval,r.pixel)]
        all(getproperty(r,Symbol(k))==n(oldrow,k) for k in
            ("background_contrast_nm","original_molecular_contrast_nm","referenced_contrast_nm")) || error("Changed paired arithmetic")
    end
    R.S.newdir(out); R.table(joinpath(out,"states.tsv"),state_rows)
    R.table(joinpath(out,"paired.tsv"),pair_rows); R.table(joinpath(out,"paired_summary.tsv"),pair_summaries)
    open(io->TOML.print(io,Dict("membership_sha256"=>g.meta["membership_sha256"],
        "config_sha256"=>g.meta["config_sha256"],"old_paired_replay_exact"=>true,
        "new_density_calculation"=>false,"reference_adopted"=>false)),joinpath(out,"summary.toml"),"w")
    foreach(println,filter(r->r.group=="common_exterior",pair_summaries))
end

function main(args)
    if args==["--help"]
        println("diagnose_molecular_exterior.jl prepare ROOT CONFIG NEW_GEOMETRY\n",
            "diagnose_molecular_exterior.jl analyze ROOT SAVED_GEOMETRY NEW_REPORT")
    elseif length(args)==4 && args[1]=="prepare"; prepare(args[2:4]...)
    elseif length(args)==4 && args[1]=="analyze"; analyze(args[2:4]...)
    else; error("Use --help"); end
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && DiagnoseMolecularExterior.main(ARGS)
