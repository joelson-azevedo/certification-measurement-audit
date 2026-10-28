"""Correspondence, content and cross-anchor diagnostics without drift p-values."""
from __future__ import annotations

import json
import numpy as np
import pandas as pd

from .linking_stocking_lord import fit_link
from .linking_transformations import fixed_link_covariance,information,probabilities,stable_difficulty,transform_parameters


def anchor_sets(metadata,items_a,items_b):
    common=[item for item in items_a if item in set(items_b)]
    if len(common)!=len(set(common)):
        raise ValueError("Duplicate common-item correspondence")
    meta=metadata.set_index("item_id")
    designated=meta.anchor_item.astype(str).str.lower().eq("true")
    if set(common)!=set(meta.loc[designated].index):raise ValueError("Common IDs and designated anchors differ")
    if not {"ANC_14","ANC_18"}<=set(common):raise ValueError("Prespecified exclusion items missing")
    return {"S20":common,"S19":[i for i in common if i!="ANC_14"],
            "S18":[i for i in common if i not in ["ANC_14","ANC_18"]]}


def align_anchors(model_a,model_b,members):
    if len(set(members))!=len(members) or not members:raise ValueError("Unique nonempty anchors required")
    try:
        ia=[model_a["items"].index(item) for item in members]
        ib=[model_b["items"].index(item) for item in members]
    except ValueError as exc:raise ValueError("Anchor missing from one form") from exc
    return model_a["d"][ia],model_a["a"][ia],model_b["d"][ib],model_b["a"][ib]


def content_coverage(metadata,models,sets):
    meta=metadata.set_index("item_id");domains=sorted(set(meta.domain));rows=[]
    for specification,members in sets.items():
        anchor_counts=meta.loc[members,"domain"].value_counts()
        for form in ["A","B"]:
            full_counts=meta.loc[models[form]["items"],"domain"].value_counts()
            for domain in domains:
                n=int(anchor_counts.get(domain,0));full=int(full_counts.get(domain,0))
                rows.append({"anchor_specification":specification,"form":form,"domain":domain,"anchor_count":n,
                             "anchor_share":n/len(members),"full_form_count":full,"full_form_share":full/len(models[form]["items"]),
                             "share_difference":n/len(members)-full/len(models[form]["items"]),
                             "number_of_anchors":len(members),"fraction_original_pool_retained":len(members)/len(sets["S20"]),
                             "scope":"supplied broad-domain labels only; detailed-content representativeness not established"})
    return pd.DataFrame(rows)


def supported_range(model_a,model_b,scale,shift,quantiles=(.01,.99),grid_limits=(-3,3)):
    a=np.quantile(model_a["eap"],quantiles)
    b=scale*np.quantile(model_b["eap"],quantiles)+shift
    lower=max(a[0],b[0],grid_limits[0]);upper=min(a[1],b[1],grid_limits[1])
    if upper<=lower:raise ValueError("No observed central proficiency overlap after preliminary linking")
    return float(lower),float(upper)


def curve_diagnostics(model_a,model_b,sets,link,theta,weights,support):
    members=sets["S20"];da,aa,db,ab=align_anchors(model_a,model_b,members)
    dt,at=transform_parameters(db,ab,link["scale_A"],link["shift_B"])
    pa,pb=probabilities(theta,da,aa),probabilities(theta,dt,at)
    delta=pb-pa;supported=(theta>=support[0])&(theta<=support[1])
    sw=weights*supported;sw/=sw.sum()
    rows=[];curves=[]
    covariance=fixed_link_covariance(model_b["covariance"],link["scale_A"],link["shift_B"])
    ka,kb=len(model_a["items"]),len(model_b["items"])
    for j,item in enumerate(members):
        ia=model_a["items"].index(item);ib=model_b["items"].index(item)
        row={"item_id":item,"a_A":aa[j],"d_A":da[j],"a_B_original":ab[j],"d_B_original":db[j],
             "a_B_transformed":at[j],"d_B_transformed":dt[j],"slope_difference_BminusA":at[j]-aa[j],
             "intercept_difference_BminusA":dt[j]-da[j],
             "se_a_A":np.sqrt(model_a["covariance"][ka+ia,ka+ia]),"se_d_A":np.sqrt(model_a["covariance"][ia,ia]),
             "se_a_B_original":np.sqrt(model_b["covariance"][kb+ib,kb+ib]),"se_d_B_original":np.sqrt(model_b["covariance"][ib,ib]),
             "se_a_B_transformed_fixed_link":np.sqrt(covariance[kb+ib,kb+ib]),"se_d_B_transformed_fixed_link":np.sqrt(covariance[ib,ib]),
             "integrated_absolute_probability_difference":float(weights@abs(delta[:,j])),
             "integrated_absolute_probability_difference_supported":float(sw@abs(delta[:,j])),
             "maximum_absolute_probability_difference":float(abs(delta[:,j]).max()),
             "maximum_absolute_probability_difference_supported":float(abs(delta[supported,j]).max()),
             "included_in_linking_objective":item in sets[link["anchor_specification"]],
             "uncertainty_scope":"source SEs and transformed SEs conditional on fixed link; estimated-link dependence addressed by uncertainty replications, not naive difference tests"}
        row.update({"member_"+name:item in anchors for name,anchors in sets.items()});rows.append(row)
        curves.append(pd.DataFrame({"item_id":item,"theta_A":theta,"P_A":pa[:,j],"P_B_linked":pb[:,j],
                                    "probability_difference_BminusA":delta[:,j],"within_observed_common_support":supported}))
    included=[members.index(item) for item in sets[link["anchor_specification"]]]
    tcc_difference=delta[:,included].sum(axis=1)
    summary={"anchor_tcc_rmse":float(np.sqrt(weights@(tcc_difference*tcc_difference))),
             "anchor_tcc_maximum_absolute_difference":float(abs(tcc_difference).max()),
             "anchor_tcc_rmse_supported":float(np.sqrt(sw@(tcc_difference*tcc_difference))),
             "anchor_tcc_maximum_absolute_difference_supported":float(abs(tcc_difference[supported]).max())}
    return pd.DataFrame(rows),pd.concat(curves,ignore_index=True),summary


def leave_one_out(model_a,model_b,members,metadata,point_links,theta,weights,support,config,safeguards):
    rows=[];attempts=[];meta=metadata.set_index("item_id")
    evaluate=(theta>=support[0])&(theta<=support[1]);w=weights[evaluate];w/=w.sum()
    for heldout in members:
        retained=[item for item in members if item!=heldout]
        counts=meta.loc[retained,"domain"].value_counts()
        adequate=len(retained)>=safeguards["minimum_anchors"] and all(counts.get(domain,0)>=safeguards["minimum_anchors_per_domain"] for domain in set(meta.loc[model_a["items"],"domain"]))
        da,aa,db,ab=align_anchors(model_a,model_b,retained)
        d0,a0,d1,a1=align_anchors(model_a,model_b,[heldout])
        for method,point in point_links.items():
            link=fit_link(da,aa,db,ab,theta,weights,config,method)
            dt,at=transform_parameters(d1,a1,link["scale_A"],link["shift_B"])
            difference=(probabilities(theta[evaluate],dt,at)-probabilities(theta[evaluate],d0,a0))[:,0]
            rows.append({"held_out_item":heldout,"method":method,"n_retained_anchors":len(retained),
                         "retained_anchor_ids":json.dumps(retained),"content_supported":adequate,
                         "interpretation":"cross-anchor consistency; same calibration samples, no independent validation or drift significance test",
                         **{key:value for key,value in link.items() if key!="attempts"},
                         "scale_change_from_primary":link["scale_A"]-point["scale_A"],
                         "shift_change_from_primary":link["shift_B"]-point["shift_B"],
                         "held_out_integrated_absolute_probability_difference":float(w@abs(difference)),
                         "held_out_maximum_absolute_probability_difference":float(abs(difference).max())})
            log=link["attempts"].copy();log.insert(0,"held_out_item",heldout);log.insert(0,"method",method);attempts.append(log)
    return pd.DataFrame(rows),pd.concat(attempts,ignore_index=True)


def anchor_information(model_a,model_b,members,link,theta,cut_theta):
    da,aa,db,ab=align_anchors(model_a,model_b,members)
    dt,at=transform_parameters(db,ab,link["scale_A"],link["shift_B"])
    full_d,full_a=transform_parameters(model_b["d"],model_b["a"],link["scale_A"],link["shift_B"])
    ia,ib=information(theta,da,aa),information(theta,dt,at)
    return pd.DataFrame({"theta_A":theta,"anchor_information_A":ia,"anchor_information_B_linked":ib,
                         "anchor_expected_A":probabilities(theta,da,aa).sum(axis=1),
                         "anchor_expected_B_linked":probabilities(theta,dt,at).sum(axis=1)}),{
        "anchor_information_A_at_A_expected42":float(information(cut_theta,da,aa)),
        "anchor_information_B_at_A_expected42":float(information(cut_theta,dt,at)),
        "anchor_to_full_information_ratio_A_at_A_expected42":float(information(cut_theta,da,aa)/information(cut_theta,model_a["d"],model_a["a"])),
        "anchor_to_full_information_ratio_B_at_A_expected42":float(information(cut_theta,dt,at)/information(cut_theta,full_d,full_a))}
