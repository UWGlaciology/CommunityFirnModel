Meltwater checks and strict mode
--------------------------------

The **bucket** and **darcy** meltwater schemes (``bucket`` and ``darcyscheme`` in ``melt.py``) do checks at the end of each call:

- A **budget** check. The water at the start of the call (liquid water in the column, plus melt, plus rain) must agree with the water at the end (liquid water, plus refreezing, plus runoff).
- A **cold-layer** check. No layer below 0 °C can hold liquid water.

The **prefsnowpack** and **resingledomain** schemes do not do these checks. Before model year 1980, the **darcy**, **prefsnowpack** and **resingledomain** options use the **bucket** scheme, which does the checks.

The checks
==========

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Scheme
     - Check
     - What it checks
   * - ``bucket``
     - ``cold_1``
     - After percolation and the first refreezing step: no layer below 0 °C holds liquid water.
   * - ``bucket``
     - ``budget_1``
     - After runoff: the water budget closes to 1e-3 m w.e.
   * - ``bucket``
     - ``cold_2``
     - After the temperature correction and the second refreezing step: no layer below 0 °C holds liquid water.
   * - ``bucket``
     - ``budget_2``
     - At the end of the call: the water budget closes to 1e-5 m w.e.
   * - ``darcy``
     - ``budget``
     - At the end of the call (all sub-steps): the water budget closes to 1e-12 m w.e.
   * - ``darcy``
     - ``cold``
     - At the end of the call: no layer below 0 °C holds liquid water.

Before each cold-layer check, the scheme removes round-off liquid water from cold layers. If each cold layer holds less than 1e-9 m w.e., the scheme sets the liquid water of all cold layers to zero. It records the removed amount as ``cold_lwc_removed``. If one or more cold layers hold more than this, the scheme removes nothing and the check fails.

The result of a failed check
============================

The **melt_strict** key in the .json configuration (see :doc:`../running/json`) sets what occurs when a check fails:

- **melt_strict** false (default). The scheme prints a message and the run continues.
- **melt_strict** true. The scheme raises ``melt.MeltCheckError`` and the run stops. The message gives the scheme, the check, the time step, and the recorded values, for example::

	MeltCheckError: bucket: check 'budget_1' failed at step 7: {'ok': False, 'residual': 0.002, 'threshold': 0.001}

Strict mode is useful for ensemble runs, because a run with a wrong water budget then does not complete.

The diagnostics record
======================

Each call of ``bucket`` or ``darcyscheme`` writes a new dictionary ``melt_diagnostics`` on the model object. It holds the results of that call only. It is not a time series.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Key
     - Value
   * - ``scheme``
     - ``'bucket'`` or ``'darcy'``.
   * - ``step``
     - The time-step index ``iii``.
   * - ``water_in``
     - The water at the start of the call: liquid water, plus melt, plus rain [m w.e.].
   * - ``lwc_end``
     - The liquid water in the column at the end of the call [m w.e.].
   * - ``refrozen``
     - The water that refroze in the call [m w.e.].
   * - ``runoff``
     - The runoff in the call [m w.e.].
   * - ``checks``
     - One entry for each check, with the check name as key. Each entry has ``ok`` (True or False). A budget check also has ``residual`` and ``threshold`` [m w.e.]. A cold-layer check also has ``n_wet_cold_layers`` and ``cold_lwc_removed`` [m w.e.].

For example, after a call of ``bucket``::

	>>> model.melt_diagnostics['checks']['budget_2']
	{'ok': True, 'residual': 1.7e-17, 'threshold': 1e-05}
