The bucket percolation scheme
-------------------------------

Overview
========

The bucket scheme (``bucket()`` in ``melt.py``) routes surface meltwater and
rainfall through a firn column for a single model time step. It is one of
several percolation schemes available in the CFM for representing how liquid
water moves downward through firn, refreezes, is retained against gravity
drainage, or runs off, selected with the **liquid** json key (see
:doc:`../running/json`).

Unlike more physically detailed approaches (the Darcy-flow scheme, see
:doc:`darcy_scheme`), the bucket scheme treats each firn layer ("node") as a
discrete reservoir with a finite holding capacity. Water in excess of that
capacity is passed downward to the next node, refrozen if cold content
allows, or blocked/ponded above impermeable ice lenses. This makes the
scheme computationally efficient, at the cost of not resolving sub-timestep
flow dynamics.

Originally coded by Vincent Verjans, with later edits by Max S.

Physical motivation
====================

Meltwater percolating through firn interacts with three limiting factors:

1. **Available pore space** -- firn has a finite air-filled volume that
   liquid water can occupy without changing the layer's thickness. Because
   liquid water (1000 kg/m\ :sup:`3`) is denser than ice (917 kg/m\
   :sup:`3`), the amount of water a layer can hold without exceeding
   solid-ice density is *less* than the raw pore volume would suggest.
2. **Retention against gravity** -- a portion of that pore space holds
   water irreducibly (via capillary/grain-scale forces) even when no
   temperature gradient or barrier is present.
3. **Refreezing** -- cold firn has a finite energy budget (cold content)
   available to freeze incoming water, releasing latent heat and warming
   the layer toward the melting point.

Water that cannot be accommodated by any of the above, and that encounters
an impermeable ice lens, either ponds above the barrier (if enabled) or runs
off.

Algorithm outline
==================

For each model time step, ``bucket()`` performs the following steps in
order:

1. **Surface melting and regridding** -- mass lost to melt is removed from
   the top of the column. Nodes are consumed from the surface downward
   until the cumulative melted mass is satisfied; the last node consumed is
   only partially melted and becomes the new surface node. The grid is then
   "refilled" at the bottom with copies of the deepest existing node,
   keeping the total node count constant (a placeholder, corrected later by
   accumulation/densification modules elsewhere in the CFM).

2. **Storage capacity calculation** -- for each node, compute:

   - *Refreezing capacity*: how much additional liquid water could be
     frozen, limited by both cold content (energy) and available pore
     space (volume).
   - *Retention capacity*: how much additional liquid water could be held
     irreducibly, using either the Coléou & Lesaffre (1998)
     density-based formulation or a fixed user-set fraction of pore space.
   - Total storage capacity = refreezing capacity + retention capacity.

3. **Impermeable node detection** -- ice lenses are identified as
   impermeable barriers to percolation, based on a density threshold
   (**RhoImp**) and, optionally, a minimum lens thickness (**ThickImp**) or
   a "bypass shallow lenses" mode (**DownToIce**). The domain's bottom node
   is always treated as impermeable.

4. **Distribution of surface liquid input** -- melt and rain input at the
   surface is distributed downward into available storage capacity, node by
   node, until either the input is fully accommodated or an impermeable
   barrier is reached (in which case the remainder is recorded as blocked
   LWC, ``LWCblocked``).

5. **Redistribution of excess LWC** -- liquid water already present in a
   node beyond its irreducible retention value ("excess" LWC) is passed
   downward into any remaining storage capacity, or blocked above the
   nearest impermeable barrier if none is found.

6. **Refreezing (pass 1)** -- excess/redistributed LWC refreezes according
   to each node's refreezing capacity; mass, density, LWC, and temperature
   are updated accordingly.

7. **Runoff and ponding** -- a user-set fraction of blocked LWC can run off
   directly; the remainder either ponds in available pore space above
   impermeable barriers (forming a perched saturated zone) or runs off
   entirely, depending on the **Ponding** setting.

8. **Zuo & Oerlemans (1996) lateral runoff** *(optional)* -- an additional
   time-dependent lateral drainage term is applied to any LWC exceeding
   irreducible retention, using a slope-dependent characteristic drainage
   timescale.

9. **Refreezing (pass 2)** -- a second refreeze pass catches any additional
   freezing made possible by the redistribution/ponding steps, which may
   have moved LWC into nodes with leftover cold content not accounted for
   in pass 1.

10. **Mass conservation checks** -- liquid water mass is checked for
    conservation (liquid remaining + total refrozen + total runoff, versus
    the initial liquid input) twice: once after pass 1 (loose tolerance),
    and once at the very end (tight tolerance).

Configuration options
======================

These are the bucket-scheme-relevant json keys; see :doc:`../running/json`
for types and defaults.

.. list-table::
   :header-rows: 1

   * - Key
     - Description
   * - **ColeouLesaffre**
     - If true, irreducible water content is computed via the Coléou &
       Lesaffre (1998) density-based formulation (Langen et al. 2017, Eqs.
       3-4). If false, a fixed fraction (**IrrVal**) of available pore space
       is used instead.
   * - **IrrVal**
     - Fixed irreducible water fraction of pore space (used only if
       **ColeouLesaffre** is false).
   * - **RhoImp**
     - Density threshold above which a node is considered part of an
       impermeable ice lens.
   * - **DownToIce**
     - If true, water bypasses all shallower ice lenses and is only blocked
       once it reaches the last node before "true" ice-sheet density is
       reached. If false, impermeability is instead governed by
       **ThickImp**.
   * - **ThickImp**
     - Minimum ice-lens thickness required for a lens to be impermeable
       (used only if **DownToIce** is false). If 0, every node exceeding
       **RhoImp** is impermeable regardless of lens thickness. The
       bottom-most lens in the domain is always impermeable regardless of
       this setting.
   * - **Ponding**
     - If true, LWC blocked by an impermeable barrier can pond in available
       pore space above the barrier. If false, all blocked LWC runs off.
   * - **DirectRunoff**
     - Fraction of blocked LWC that runs off immediately, without
       attempting storage or ponding.
   * - **RunoffZuoOerlemans**
     - If true, an additional lateral runoff term (Zuo & Oerlemans, 1996)
       is applied to LWC exceeding irreducible retention, after refreezing.
   * - **Slope**
     - Surface slope used in the Zuo & Oerlemans (1996) drainage-timescale
       parameterization (Eq. 22). Best-supported interpretation is
       rise/run; the source paper does not give an explicit defining
       formula for this quantity, so this is an inference from context, not
       a directly quoted definition -- see inline code comments for
       details.
   * - **keep_firnthickness**
     - If true, new bottom nodes added during regridding are sized to
       approximately conserve total column thickness (rather than
       shrinking as surface nodes are melted away).
   * - **debug_bucket**
     - If true, enables verbose diagnostic printouts of intermediate
       LWC/runoff bookkeeping values throughout the function (see
       ``_debug_report`` below).

Helper functions
=================

Several calculations that recur multiple times within ``bucket()`` are
factored into standalone helper functions, defined above ``bucket()`` in
``melt.py``. These are shared with :doc:`darcy_scheme`.

``_available_pore_space(rho, dz, rho_i=RHO_I, cap_density=RHO_I)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Computes the pore space available for liquid water storage per node,
following the "saturated water content" formulation of Wever et al. (2014),
Eq. 9 (see also discussion in Yamaguchi et al., 2010). The result is capped
so that filling the available space with water cannot cause a node's
density to exceed ``cap_density`` -- this enforces the physical constraint
that, since liquid water is denser than ice, a fixed-volume layer cannot
hold more mass (firn + water, assuming any of that water later refreezes)
than if it were pure solid ice.

``_irreducible_lwc(rho, phivol_av, RhoImp, rho_i=RHO_I, coleou_lesaffre=True, irr_val=0.02)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Computes the irreducible liquid water content each node can hold against
gravity drainage, using either the Coléou & Lesaffre (1998)
density-dependent formulation (via Langen et al., 2017, Eqs. 3-4), or a
fixed fraction of available pore space. Nodes at or above the
impermeability density threshold (``RhoImp``) are assigned zero irreducible
LWC, since they are treated as solid ice incapable of holding retained
water.

``_find_impermeable_layers(rho, dz, RhoImp, ThickImp=0., DownToIce=False, force_bottom_impermeable=False)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Identifies which nodes in the column are impermeable to percolation, and
the index of the top node of each impermeable lens/block (``imptop``).
Three modes are supported (see the Configuration Options table above for
**DownToIce**/**ThickImp** behavior). If no layer qualifies as impermeable
under any mode, the last layer in the domain is forced impermeable so the
percolation schemes still have a terminating barrier.

Shared between the bucket and Darcy schemes, but with one deliberate
difference, controlled by ``force_bottom_impermeable``: ``bucket()`` passes
``True``, which additionally forces the bottom node impermeable even when
other layers already qualify -- this reproduces the scheme's historical
behavior and is relied on by the unguarded storage-capacity lookups further
down in ``bucket()``. ``darcyscheme()`` leaves it ``False``. See
:doc:`darcy_scheme` for the Darcy side of this.

``_debug_report(debug, label, value)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A lightweight print wrapper used to optionally report intermediate
diagnostic values (e.g. blocked LWC at various stages of the redistribution
logic) without cluttering the main routine with conditional print
statements. Controlled by the **debug_bucket** config flag. Note: this is
distinct from the mass-conservation and dry-cold-firn warning prints
elsewhere in ``bucket()``, which always print regardless of this flag,
since they indicate potential bugs rather than routine diagnostics.

Key variables reference
=========================

A selection of variable names used throughout ``bucket()``, for readers
cross-referencing the code:

.. list-table::
   :header-rows: 1

   * - Variable
     - Meaning
   * - ``mass_sum``
     - Cumulative mass from the surface downward [kg]; used to locate the
       new surface node after melting. Local to this function (not stored
       on ``self``).
   * - ``i_surf_new`` *(formerly ``ind1``)*
     - Index of the partially-melted node that becomes the new surface
       after melting.
   * - ``dh_melt``
     - Total surface height lost to melting this time step [m] (negative
       by convention).
   * - ``phivol_av``
     - Available pore space per node [m] (see ``_available_pore_space``).
   * - ``LWCirr``
     - Irreducible LWC per node [m] (see ``_irreducible_lwc``).
   * - ``LWC_excess``
     - LWC beyond the irreducible retention value in each node [m] --
       water that is mobile and must be redistributed, refrozen, or
       blocked.
   * - ``refr_cap_from_cc``
     - Refreezing capacity from cold content alone, before intersecting
       with the pore-space limit [m w.e.].
   * - ``refr_cap``
     - Total refreezing capacity per node, after intersecting
       ``refr_cap_from_cc`` with available pore space [m w.e.].
   * - ``stcap``
     - Total storage capacity per node for *additional* incoming LWC:
       refreezing capacity + retention capacity [m].
   * - ``imp``
     - Indices of impermeable nodes.
   * - ``imptop``
     - Index of the top node of each impermeable lens/block (returned by
       ``_find_impermeable_layers``, largely unused in ``bucket()``).
   * - ``i_input_bottom`` *(formerly ``ii0``)*
     - Bottom-most node index that receives distributed surface liquid
       input.
   * - ``LWCblocked``
     - LWC blocked above an impermeable barrier, pending direct runoff,
       ponding, or full runoff.
   * - ``LWC_redistributed`` *(formerly ``LWC1``)*
     - Working copy of LWC as it is modified by excess-LWC
       redistribution/blocking.
   * - ``storage_excess`` *(formerly ``storage1``)*
     - LWC stored in each node specifically via excess-LWC redistribution
       (as opposed to ``storageinp``, surface-input storage).
   * - ``stcap_remaining`` *(formerly ``stcap1``)*
     - Storage capacity remaining after surface melt input has already
       consumed some.
   * - ``i_scan`` *(formerly ``jj0``)*
     - Moving cursor index in the excess-LWC redistribution loop; advances
       downward through the column.
   * - ``i_target`` *(formerly ``jj1``)*
     - Node that the current redistribution segment's excess LWC is being
       routed to.
   * - ``refrozentot``
     - Running total of refrozen water across both refreeze passes [m
       w.e.].
   * - ``runofftot``
     - Running total of runoff (direct + blocked + lateral) [m w.e.].
   * - ``liqmcinit`` / ``liqmcfinal``
     - Total liquid water mass conservation bookkeeping, initial and final
       [m].

Known limitations & caveats
=============================

- **Regridding placeholder nodes**: new nodes appended at the bottom of the
  column during regridding are copies of the deepest existing node's
  properties, not physically derived fresh accumulation. This is expected
  to be corrected by other CFM modules (density/accumulation) in
  subsequent steps, not by ``bucket()`` itself.

- **Fixed-volume assumption**: ``dz`` is never modified by refreezing
  within ``bucket()`` -- refreezing only increases density (``mass/dz``) in
  place. This is why pore-space calculations are capped at ice density (see
  ``_available_pore_space``): a fixed-volume layer cannot physically hold
  more mass than solid ice would represent at that volume.

- **Two refreeze passes**: refreezing is applied twice per time step
  (steps 6 and 9 in the Algorithm Outline above) because redistribution,
  ponding, and lateral runoff can move LWC into nodes with unused cold
  content left over from the first pass. This is a deliberate design
  choice, not a redundant duplication.

- **Zuo & Oerlemans (1996) slope parameter**: the **Slope** config option
  is best understood as a dimensionless rise/run value, but this is
  inferred from context (default values, dimensional consistency with
  Eq. 22's constants) rather than an explicit definition quoted from the
  source paper. See inline code comments near its definition for details.

- **Debug output**: routine diagnostic printouts are gated behind the
  **debug_bucket** config flag (see ``_debug_report``), but
  mass-conservation and dry-cold-firn warnings always print regardless of
  this flag, since they indicate potential bugs rather than expected
  behavior.

- **Forcing the bottom node impermeable**: ``bucket()`` always treats the
  domain's bottom node as impermeable (via
  ``force_bottom_impermeable=True``, see Helper Functions above), even if
  its density is well below **RhoImp**. This reproduces the model's
  historical behavior and is required by downstream lookups in this
  function, but is not obviously the physically correct choice for a
  domain that does not extend to full ice density -- a non-dense bottom
  layer could instead be allowed to let water leave the domain as
  drainage. This is an open physics question, not resolved here, and
  interacts with **Ponding**: water could pond on a sub-ice-density
  "impermeable" bottom layer.

References
==========

- Coléou, C. and Lesaffre, B. (1998). Irreducible water saturation in
  snow: experimental results in a cold laboratory. *Annals of Glaciology*,
  26, 64-68.
- Langen, P.L. et al. (2017). Liquid Water Flow and Retention on the
  Greenland Ice Sheet in the Regional Climate Model HIRHAM5: Local and
  Large-Scale Impacts. *Frontiers in Earth Science*, 4, 110.
- Wever, N. et al. (2014). Verification of the multi-layer SNOWPACK model
  with different water transport schemes. *The Cryosphere*, 8, 2687-2705.
- Yamaguchi, S. et al. (2010). Water retention curve of snow with
  different grain sizes. *Cold Regions Science and Technology*, 64(2),
  87-93.
- Zuo, Z. and Oerlemans, J. (1996). Modelling albedo and specific balance
  of the Greenland ice sheet: calculations for the Søndre
  Strømfjord transect. *Journal of Glaciology*, 42(141), 305-317.
