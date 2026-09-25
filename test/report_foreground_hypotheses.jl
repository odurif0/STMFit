#!/usr/bin/env julia
# Descriptive old hypotheses AFTER foreground construction, never composition fitting.
module ReportForegroundHypotheses
include(joinpath(@__DIR__,"verify_image_foreground.jl"))
const V=VerifyImageForeground
const t=V.t
const i=V.i
const b=V.b

function main(args=ARGS)
    length(args)==4 || error("report_foreground_hypotheses.jl REPO SAVED_RUN VERIFIED_REPORT NEW_TSV")
    root,run,report,out=args; ispath(out) && error("Output exists")
    c=V.D.F.settings(joinpath(run,"settings.toml")); hashes=t(joinpath(report,"prediction_input_hashes.tsv"))
    hashby=Dict(r["input"]=>r["sha256"] for r in hashes)
    counts=Dict{Tuple{String,String,Int,String,String,Int},Int}()
    for file in getindex.(t(joinpath(run,"files.tsv")),"file"),view in ("fwd","bwd")
        stem=splitext(file)[1]; old=joinpath(root,c["model"]["transfer_run"],stem,view)
        isfile(joinpath(old,"patches.tsv")) || continue
        identified=V.TOML.parsefile(joinpath(old,"calibration.toml"))["identified"]
        anchors=t(joinpath(run,stem,view,"anchor_membership.tsv")); byanchor=Dict(i(r,"anchor")=>r for r in anchors)
        for (arm,input,key) in (("original_chemical",old,"chemical"),
            ("plane_chemical",joinpath(root,c["model"]["controls_run"],stem,view),"plane_chemical"))
            path=joinpath(input,"patches.tsv")
            V.R.sha(path)==hashby[path] || error("Prediction file changed after verification")
            patches=filter(r->r["arm"]==key,t(path))
            for patch in patches
                a=byanchor[i(patch,"anchor")]; type=i(patch,"type"); type in (0,1) || error("Invalid physical hypothesis")
                for (mask,prefix) in (("initial","initial"),("joint_background","joint"))
                    group=!b(a,prefix*"_available") ? "unavailable" : b(a,prefix*"_core") ? "retained" : "rejected"
                    for subset in (identified ? ("all_scored","identified") : ("all_scored",)),g in ("all",group)
                        k=(subset,mask,i(patch,"interval"),g,arm,type); counts[k]=get(counts,k,0)+1
                    end
                end
            end
        end
    end
    rows=NamedTuple[]
    for subset in ("all_scored","identified"),mask in ("initial","joint_background"),interval in 1:3,
        group in ("all","retained","rejected","unavailable"),arm in ("original_chemical","plane_chemical")
        zero=get(counts,(subset,mask,interval,group,arm,0),0); one=get(counts,(subset,mask,interval,group,arm,1),0)
        push!(rows,(;subset,mask,interval,group,arm,glcn_hypotheses=zero,glcnac_hypotheses=one,total=zero+one,
            glcn_fraction=zero+one>0 ? zero/(zero+one) : NaN))
    end
    # Counts must partition precisely the independently verified prediction support.
    reference=t(joinpath(report,"prediction_summary.tsv"))
    for row in rows
        r=only(filter(r->r["subset"]==row.subset && r["mask"]==row.mask && i(r,"interval")==row.interval &&
            r["group"]==row.group && r["arm"]==row.arm,reference))
        row.total==i(r,"patches") || error("Changed patch population")
    end
    V.R.table(out,rows)
    println("Saved old physical hypothesis counts; not units, composition estimates or label truth.")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && ReportForegroundHypotheses.main()
