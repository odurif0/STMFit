#!/usr/bin/env julia
# Aggregate verified saved tables; no mask, fit, parameter or label input.
module SummarizeImageForeground
using Statistics
include(joinpath(@__DIR__,"verify_image_foreground.jl"))
const V=VerifyImageForeground
const t=V.t
const n=V.n
const i=V.i
const b=V.b

function main(args=ARGS)
    length(args)==2 || error("summarize_image_foreground.jl VERIFIED_REPORT NEW_SUMMARY_DIR")
    input,out=args; V.R.S.newdir(out)
    support=t(joinpath(input,"support.tsv")); pairs=t(joinpath(input,"pair_overlap.tsv"))
    cases=t(joinpath(input,"prediction_cases.tsv")); statuses=t(joinpath(input,"statuses.tsv"))
    support_rows=NamedTuple[]; overlap_rows=NamedTuple[]; change_rows=NamedTuple[]; comparisons=NamedTuple[]
    for arm in ("initial","joint_background")
        rr=filter(r->r["arm"]==arm,support); available=sum(i.(rr,"available_pixels")); foreground=sum(i.(rr,"foreground_pixels"))
        push!(support_rows,(;arm,views=length(rr),files=length(unique(getindex.(rr,"file"))),raw_pixels=sum(i.(rr,"raw_pixels")),
            available_pixels=available,foreground_pixels=foreground,foreground_fraction=available>0 ? foreground/available : NaN,
            unavailable_views=count(r->i(r,"available_pixels")==0,rr),empty_views=count(r->i(r,"foreground_pixels")==0,rr),
            components=sum(i.(rr,"components")),edge_components=sum(i.(rr,"edge_components")),missing_components=sum(i.(rr,"missing_components"))))
    end
    for subset in ("all_views","identified"),kind in ("own","common_methods")
        selected=filter(r->r["support"]==kind && (subset=="all_views" || b(r,"identified")),pairs)
        for arm in ("initial","joint_background")
            rr=filter(r->r["arm"]==arm,selected); intersection=sum(i.(rr,"intersection")); union_pixels=sum(i.(rr,"union_pixels"))
            valid=filter(r->isfinite(n(r,"iou")),rr)
            push!(overlap_rows,(;subset,support=kind,arm,views=length(rr),files=length(unique(getindex.(rr,"file"))),
                undefined_iou_views=length(rr)-length(valid),comparison_pixels=sum(i.(rr,"pixels")),intersection,union_pixels,
                pooled_iou=union_pixels>0 ? intersection/union_pixels : NaN,median_iou=isempty(valid) ? NaN : median(n.(valid,"iou"))))
        end
        initial=Dict((r["file"],r["view"])=>r for r in selected if r["arm"]=="initial")
        joint=Dict((r["file"],r["view"])=>r for r in selected if r["arm"]=="joint_background")
        Set(keys(initial))==Set(keys(joint)) || error("Different pair populations")
        differences=[n(joint[k],"iou")-n(initial[k],"iou") for k in sort(collect(keys(initial))) if isfinite(n(joint[k],"iou")) && isfinite(n(initial[k],"iou"))]
        kind=="common_methods" && !all(i(initial[k],"pixels")==i(joint[k],"pixels") for k in keys(initial)) && error("Different common coverage")
        push!(change_rows,(;subset,support=kind,paired_views=length(differences),undefined_pairs=length(initial)-length(differences),
            improved=count(>(0),differences),worsened=count(<(0),differences),unchanged=count(==(0),differences),
            median_iou_change=isempty(differences) ? NaN : median(differences)))
    end
    for subset in ("all_scored","identified"),mask in ("initial","joint_background"),interval in 1:3,group in ("all","retained","rejected","unavailable")
        selected=filter(r->r["mask"]==mask && i(r,"interval")==interval && r["group"]==group &&
            (subset=="all_scored" || b(r,"identified")),cases)
        base=Dict((r["file"],r["view"])=>r for r in selected if r["arm"]=="local_plane")
        chemical=Dict((r["file"],r["view"])=>r for r in selected if r["arm"]=="plane_chemical")
        Set(keys(base))==Set(keys(chemical)) || error("Different prediction populations")
        all(i(base[k],"test_pixels")==i(chemical[k],"test_pixels") for k in keys(base)) || error("Different prediction coverage")
        keys_used=sort([k for k in keys(base) if i(base[k],"test_pixels")>0])
        byfile=Dict{String,Float64}()
        for k in keys_used; byfile[k[1]]=get(byfile,k[1],0.)+n(chemical[k],"test_sse_nm2")-n(base[k],"test_sse_nm2"); end
        push!(comparisons,(;subset,mask,interval,group,scored_views=length(keys_used),scored_files=length(byfile),
            chemical_better_views=count(k->n(chemical[k],"test_sse_nm2")<n(base[k],"test_sse_nm2"),keys_used),
            chemical_better_files=count(<(0),values(byfile))))
    end
    V.R.table(joinpath(out,"support_summary.tsv"),support_rows)
    V.R.table(joinpath(out,"overlap_summary.tsv"),overlap_rows)
    V.R.table(joinpath(out,"overlap_changes.tsv"),change_rows)
    V.R.table(joinpath(out,"chemical_plane_comparisons.tsv"),comparisons)
    V.R.table(joinpath(out,"status_counts.tsv"),[(;status,views=count(r->r["status"]==status,statuses)) for status in sort(unique(getindex.(statuses,"status")))])
    println("Saved descriptive summaries; no choice of mask, isovalue, count or chemical class.")
end
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && SummarizeImageForeground.main()
