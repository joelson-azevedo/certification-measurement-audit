"""Directional A-grid Haebara item-characteristic-curve linking."""
from .linking_stocking_lord import fit_link,objective_gradient


def haebara(d_a,a_a,d_b,a_b,theta,weights,config=None,starts=None):
    """Minimize weighted SUM of item squared discrepancies (directional A metric).

    No reversed-direction term or implicit symmetric objective is included.
    Normalization by anchor count changes numerical units, not its minimizer.
    """
    return fit_link(d_a,a_a,d_b,a_b,theta,weights,config,"Haebara-directional-A",starts)


def haebara_objective(parameters,d_a,a_a,d_b,a_b,theta,weights):
    return objective_gradient(parameters,d_a,a_a,d_b,a_b,theta,weights,"Haebara-directional-A")
