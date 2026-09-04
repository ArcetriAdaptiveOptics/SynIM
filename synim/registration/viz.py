"""
Schematic (altitude) visualization of a `System`'s geometry.

`plot_altitude_schematic` is deliberately a 2D SIDE-VIEW SLICE (like
Fig. 2a of the SPIE paper), not a full 3D rendering: the pupil is drawn
at altitude 0, each DM as a horizontal segment at its own altitude, and
each WFS's guide star as a chief ray from the pupil (z=0) up to the guide
star (finite altitude for an LGS, "infinity" - an arrow off the top of
the plot - for an NGS), using the physical, un-rescaled ray position
``theta_rad * z`` (theta being the guide star's angular position
projected onto `slice_axis`). This directly shows *why* a DM at a given
altitude only sees part of an LGS's cone, without trying to reproduce the
internal (cone-effect-rescaled) pixel-shift numbers of
`gs_parallax_transform` geometrically - those stay exact in
`model.py`/`analysis.py`; this module is for visual intuition, not for
extracting numbers from the plot.

The schematic itself only carries short name tags (colour-coded); the
actual mis-registration numbers go in a separate table
(`plot_mis_registration_table`), since packing them as inline plot text
turns unreadable as soon as a system has more than one or two elements.

`plot_dm_footprint`/`plot_dm_footprints` add the complementary TOP-DOWN
view at a given DM's altitude (pupil + every guide star's footprint
circle there), adapted from PASSATA's `compute_mcao_geom.pro` and Fig. 2a
of the SPIE paper to the `System`'s actual (possibly non-symmetric,
mis-registered) guide star constellation.

`plot_system_overview` puts the schematic, the table and (optionally,
if a pupil diameter is given) the footprint panels together in one figure.

`plot_mode_bar`/`plot_mode_gs_quiver` visualize one SVD mode of a
sensitivity matrix (`analysis.svd_of_jacobian`, one row of `Vt`, i.e. one
coefficient per `ParameterSpec`): a bar chart is a complete, unfiltered
view (unlike `analysis.describe_mode`'s small-coefficient cutoff, which
can hide contributions that are individually small but collectively
needed for an exact degeneracy), and the sky quiver is a geometric view
of the guide-star-position part of a mode, when it has one.
"""
import numpy as np

__all__ = [
    "plot_altitude_schematic", "plot_mis_registration_table",
    "plot_dm_footprint", "plot_dm_footprints", "plot_system_overview",
    "plot_mode_bar", "plot_mode_gs_quiver",
]

ARCSEC2RAD = np.pi / 180 / 3600


def _fmt(x, sig=3):
    """Format a float with `.{sig}g`, snapping -0 to 0 for readability."""
    x = float(x) + 0.0  # snap -0.0 to 0.0
    return f"{x:.{sig}g}"


def _element_colors(system):
    """A stable name -> matplotlib color mapping, shared by the schematic
    and the table so the same DM/WFS is recognizable in both."""
    import matplotlib.pyplot as plt
    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    names = list(system.dms.keys()) + list(system.wfss.keys())
    return {name: cycle[i % len(cycle)] for i, name in enumerate(names)}


def plot_altitude_schematic(system, ax=None, slice_axis=0, half_width=1.0,
                             height_margin=1.3, show_labels=False, legend=True,
                             fontsize=9, ngs_ray_length=None, colors=None):
    """
    Draw the altitude schematic of `system` on `ax` (a new figure/axes is
    created if `ax` is None) - see `plot_mis_registration_table` for the
    actual mis-registration numbers.

    Parameters
    ----------
    system : registration.model.System
    ax : matplotlib.axes.Axes, optional
    slice_axis : int
        Which guide star position component (0=x, 1=y) to project onto
        the horizontal axis - this is a 2D slice, so guide stars that
        differ only in the other component will appear to overlap (and,
        with several guide stars at symmetric angles - e.g. an evenly
        spaced ring - end up literally on top of each other).
    half_width : float
        Half-width of the pupil/DM reference bars, in the same unit as
        the plotted guide star ray positions (arbitrary/cosmetic:
        `System`/`DM` carry no telescope diameter).
    height_margin : float
        The plot's vertical extent is `height_margin` times the highest
        DM/finite-GS altitude.
    show_labels : bool
        If True, tag each DM/WFS with its short name next to the element
        (inline text) - gets crowded with more than a handful of
        elements, especially combined with the `slice_axis` overlap above;
        `legend` (the default) is the readable alternative for those cases.
    legend : bool
        If True (default), add a legend mapping colour -> DM/WFS name
        instead of (or alongside) inline labels.
    ngs_ray_length : float, optional
        How far up to draw an NGS's (infinite-height) chief ray; defaults
        to the plot's top margin.
    colors : dict, optional
        name -> matplotlib color (default: `_element_colors(system)`, also
        used by `plot_mis_registration_table` to keep colors consistent
        across both panels).

    Returns
    -------
    ax : the matplotlib Axes used.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 6))
    if colors is None:
        colors = _element_colors(system)

    dms = list(system.dms.values())
    finite_heights = [dm.height for dm in dms]
    finite_heights += [wfs.guide_star.actual_height for wfs in system.wfss.values()
                        if not np.isinf(wfs.guide_star.actual_height)]
    max_height = max(finite_heights) if finite_heights else 0.0
    if max_height <= 0.0:
        max_height = 1.0  # degenerate case (e.g. a single pupil DM + NGS only)
    top = max_height * height_margin
    if ngs_ray_length is None:
        ngs_ray_length = top

    # Pupil, at altitude 0.
    ax.plot([-half_width, half_width], [0.0, 0.0], color="0.2", lw=2, zorder=1)
    ax.text(-half_width * 1.05, 0.0, "pupil", va="center", ha="right", fontsize=fontsize)

    # Each DM as a horizontal bar, with a small marker for its own shift.
    for dm in dms:
        color = colors[dm.name]
        ax.plot([-half_width, half_width], [dm.height, dm.height], color=color, lw=2,
                 zorder=1, label=dm.name)
        ax.plot(dm.shift[slice_axis], dm.height, marker="^", color=color, markersize=8, zorder=3)
        if show_labels:
            ax.annotate(dm.name, (half_width, dm.height), xytext=(4, 4),
                         textcoords="offset points", fontsize=fontsize, color=color,
                         fontweight="bold")

    # Each WFS's guide star chief ray, from the pupil up to the GS.
    for wfs in system.wfss.values():
        color = colors[wfs.name]
        gs = wfs.guide_star
        theta_rad = gs.actual_position[slice_axis] * ARCSEC2RAD
        gs_height = gs.actual_height

        if np.isinf(gs_height):
            z_top = ngs_ray_length
            ax.annotate("", xy=(theta_rad * z_top, z_top), xytext=(0.0, 0.0),
                        arrowprops=dict(arrowstyle="->", color=color, lw=1.5), zorder=2)
            # `annotate`'s arrow is not legend-tracked: an invisible proxy
            # line carries the label instead.
            ax.plot([], [], color=color, lw=1.5, label=wfs.name)
            label_pos = (theta_rad * z_top, z_top)
        else:
            ax.plot([0.0, theta_rad * gs_height], [0.0, gs_height], color=color, lw=1.5,
                     zorder=2, label=wfs.name)
            ax.plot(theta_rad * gs_height, gs_height, marker="*", color=color,
                     markersize=14, zorder=3)
            label_pos = (theta_rad * gs_height, gs_height)

        # Mark where this GS's chief ray crosses each DM.
        for dm in dms:
            ax.plot(theta_rad * dm.height, dm.height, marker="o", color=color,
                     markersize=5, fillstyle="none", zorder=3)

        if show_labels:
            ax.annotate(wfs.name, label_pos, xytext=(4, 4), textcoords="offset points",
                         fontsize=fontsize, color=color, fontweight="bold")

    ax.set_xlabel(f"position along slice_axis={slice_axis} (chief-ray projection)")
    ax.set_ylabel("altitude [m]")
    ax.set_ylim(-0.05 * top, top * 1.05)
    ax.margins(x=0.15)
    ax.set_title("Altitude schematic")
    if legend:
        ax.legend(fontsize=fontsize - 1, loc="best", framealpha=0.9, ncol=2)
    return ax


def plot_mis_registration_table(system, ax=None, colors=None, fontsize=9):
    """
    Render the current mis-registration parameters of every DM and WFS/GS
    in `system` as a compact table (two stacked blocks: DM, then WFS),
    color-coded to match `plot_altitude_schematic`.

    Parameters
    ----------
    system : registration.model.System
    ax : matplotlib.axes.Axes, optional
        A new figure/axes is created if None; the axes' frame/ticks are
        hidden (this is a table, not a plot).
    colors : dict, optional
        name -> matplotlib color (default: `_element_colors(system)`).

    Returns
    -------
    ax : the matplotlib Axes used.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    if colors is None:
        colors = _element_colors(system)
    ax.axis("off")

    dm_header = ["DM", "height [m]", "shift x", "shift y", "rot [deg]", "mag", "anam45", "anam90"]
    dm_rows = [[dm.name, f"{dm.height:.0f}", _fmt(dm.shift[0]), _fmt(dm.shift[1]),
                _fmt(dm.rotation), _fmt(dm.magnification, 4), _fmt(dm.anamorphosis_45, 4),
                _fmt(dm.anamorphosis_90, 4)]
               for dm in system.dms.values()]

    wfs_header = ["WFS", "guide star", "shift x", "shift y", "rot [deg]", "mag", "anam45", "anam90"]
    wfs_rows = []
    gs_header = ["guide star", "pos [\"]", "pos err [\"]", "height [m]", "height err [m]"]
    gs_rows = []
    for wfs in system.wfss.values():
        gs = wfs.guide_star
        wfs_rows.append([wfs.name, gs.name, _fmt(wfs.shift[0]), _fmt(wfs.shift[1]),
                          _fmt(wfs.rotation), _fmt(wfs.magnification, 4),
                          _fmt(wfs.anamorphosis_45, 4), _fmt(wfs.anamorphosis_90, 4)])
        height_str = "inf" if np.isinf(gs.height) else _fmt(gs.height, 6)
        gs_rows.append([gs.name, f"({_fmt(gs.position[0])}, {_fmt(gs.position[1])})",
                         f"({_fmt(gs.position_shift[0])}, {_fmt(gs.position_shift[1])})",
                         height_str, _fmt(gs.height_shift)])

    def _row_colors(rows, name_col, names_to_colors):
        return [[names_to_colors.get(row[name_col], "black")] + ["black"] * (len(row) - 1)
                for row in rows]

    blocks = [("Deformable mirrors", dm_header, dm_rows, 0),
              ("Wavefront sensors", wfs_header, wfs_rows, 0),
              ("Guide stars", gs_header, gs_rows, 0)]

    y = 1.0
    row_h = 0.045
    for title, header, rows, name_col in blocks:
        ax.text(0.0, y, title, fontsize=fontsize + 1, fontweight="bold",
                transform=ax.transAxes, va="top")
        y -= row_h * 1.3
        if not rows:
            continue
        table = ax.table(cellText=rows, colLabels=header, cellLoc="center",
                          bbox=[0.0, y - row_h * len(rows), 1.0, row_h * len(rows)],
                          transform=ax.transAxes)
        table.auto_set_font_size(False)
        table.set_fontsize(fontsize)
        table.auto_set_column_width(col=list(range(len(header))))
        for (r, c), cell in table.get_celld().items():
            if r == 0:
                cell.set_text_props(fontweight="bold")
            elif c == name_col:
                cell.set_text_props(color=colors.get(rows[r - 1][name_col], "black"),
                                     fontweight="bold")
        y -= row_h * len(rows) + row_h

    return ax


def plot_dm_footprint(system, dm_name, pupil_diameter, ax=None, colors=None,
                       technical_fov_arcsec=None, fontsize=9):
    """
    Top-down view, AT ONE DM'S ALTITUDE, of the pupil and of every WFS's
    guide-star footprint on that DM - the same picture as
    `compute_mcao_geom.pro`'s `display` keyword and Fig. 2a of the SPIE
    paper, generalized from an idealized evenly-spaced guide-star ring to
    the `System`'s actual (possibly mis-registered) guide star positions.

    For guide star `gs` and this DM at height `h`:
      - footprint center = chief-ray offset ``actual_position * h``
        (same physical ray used by `plot_altitude_schematic`);
      - footprint diameter = ``pupil_diameter`` for an NGS (infinite
        height: the cone never converges, so the full pupil is seen),
        or ``pupil_diameter * (H - h) / H`` for a finite-height source at
        height H > h (the cone-effect shrink of `compute_mcao_geom.pro`'s
        `lgs_patch_diam`).

    The DM's own shift mis-registration (`DM.shift`, in `pixel_pitch`
    units) is shown as a triangle marker, converted to the same physical
    (meter) unit as everything else via `system.pixel_pitch`.

    Parameters
    ----------
    system : registration.model.System
    dm_name : str
    pupil_diameter : float
        Telescope pupil diameter, in meters (or any consistent physical
        unit - `System`/`DM` carry no telescope diameter of their own).
    ax : matplotlib.axes.Axes, optional
    colors : dict, optional
        name -> matplotlib color (default: `_element_colors(system)`).
    technical_fov_arcsec : float, optional
        If given, also draws the (dashed) metapupil circle needed to
        cover this full technical FoV at this DM's height (independent of
        any specific guide star) - `compute_mcao_geom.pro`'s `meta_diam`.

    Returns
    -------
    ax : the matplotlib Axes used.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))
    if colors is None:
        colors = _element_colors(system)

    dm = system.dms[dm_name]

    if technical_fov_arcsec is not None:
        meta_diam = pupil_diameter + 2.0 * dm.height * np.tan(
            0.5 * technical_fov_arcsec * ARCSEC2RAD)
        ax.add_patch(plt.Circle((0, 0), meta_diam / 2.0, fill=False,
                                 color="0.5", ls=":", lw=1.5))

    ax.add_patch(plt.Circle((0, 0), pupil_diameter / 2.0, fill=False, color="0.2", ls="--"))

    for wfs in system.wfss.values():
        color = colors[wfs.name]
        gs = wfs.guide_star
        center = gs.actual_position * ARCSEC2RAD * dm.height
        if np.isinf(gs.actual_height):
            footprint_diam = pupil_diameter
        else:
            footprint_diam = pupil_diameter * max(gs.actual_height - dm.height, 0.0) / gs.actual_height

        ax.add_patch(plt.Circle(center, footprint_diam / 2.0, fill=False, color=color, lw=2))
        ax.plot(*center, marker="+", color=color, markersize=10, label=wfs.name)

    dm_shift_m = np.array(dm.shift) * system.pixel_pitch
    ax.plot(*dm_shift_m, marker="^", color="0.2", markersize=8)

    # A legend rather than per-marker text: with several guide stars close
    # together (small angular separation vs. altitude) inline labels just
    # overlap and become unreadable.
    if system.wfss:
        ax.legend(fontsize=fontsize - 1, loc="best", framealpha=0.9)

    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(f"{dm_name} footprint (h={dm.height:.0f} m)", fontsize=fontsize + 1)
    ax.relim()
    ax.autoscale_view()
    ax.margins(0.15)
    return ax


def plot_dm_footprints(system, pupil_diameter, dm_names=None, colors=None,
                        technical_fov_arcsec=None, fontsize=9, figsize=None, axes=None):
    """
    `plot_dm_footprint` for several (default: all) DMs of `system`, laid
    out in a row of subplots sharing the same colour mapping.

    Returns (fig_or_None, axes) - `fig` is None if `axes` was provided.
    """
    import matplotlib.pyplot as plt

    if dm_names is None:
        dm_names = list(system.dms.keys())
    if colors is None:
        colors = _element_colors(system)

    fig = None
    if axes is None:
        fig, axes = plt.subplots(1, max(len(dm_names), 1),
                                  figsize=figsize or (4.5 * max(len(dm_names), 1), 4.5))
        if len(dm_names) == 1:
            axes = [axes]

    for ax, dm_name in zip(axes, dm_names):
        plot_dm_footprint(system, dm_name, pupil_diameter, ax=ax, colors=colors,
                           technical_fov_arcsec=technical_fov_arcsec, fontsize=fontsize)

    return fig, axes


def plot_system_overview(system, slice_axis=0, pupil_diameter=None, figsize=None,
                          fontsize=9, schematic_kwargs=None, table_kwargs=None,
                          footprint_kwargs=None):
    """
    Convenience: the altitude schematic and the mis-registration table
    side by side on top, and (if `pupil_diameter` is given) one top-down
    DM-footprint panel per DM (`plot_dm_footprints`) in a row underneath -
    all sharing the same colour mapping.

    Returns (fig, axes) where `axes` is a dict with keys "schematic",
    "table" and (if `pupil_diameter` is given) "footprints" (a list, one
    Axes per DM, in `system.dms` order).
    """
    import matplotlib.pyplot as plt

    colors = _element_colors(system)
    n_dms = len(system.dms)

    if pupil_diameter is None:
        if figsize is None:
            figsize = (13, 6)
        fig, (ax_schem, ax_table) = plt.subplots(
            1, 2, figsize=figsize, gridspec_kw={"width_ratios": [1.4, 1.0]})
        axes = {"schematic": ax_schem, "table": ax_table}
    else:
        # Two independent nested grids (top: schematic+table, bottom: one
        # column per DM) so the table's width never shrinks with n_dms.
        if figsize is None:
            figsize = (max(13, 3.2 * max(n_dms, 1)), 11)
        fig = plt.figure(figsize=figsize)
        outer = fig.add_gridspec(2, 1, height_ratios=[1.3, 1.0])
        top = outer[0].subgridspec(1, 2, width_ratios=[1.4, 1.0])
        bottom = outer[1].subgridspec(1, max(n_dms, 1))

        ax_schem = fig.add_subplot(top[0, 0])
        ax_table = fig.add_subplot(top[0, 1])
        footprint_axes = [fig.add_subplot(bottom[0, i]) for i in range(n_dms)]
        plot_dm_footprints(system, pupil_diameter, colors=colors, fontsize=fontsize,
                            axes=footprint_axes, **(footprint_kwargs or {}))
        axes = {"schematic": ax_schem, "table": ax_table, "footprints": footprint_axes}

    plot_altitude_schematic(system, ax=ax_schem, slice_axis=slice_axis,
                             fontsize=fontsize, colors=colors, **(schematic_kwargs or {}))
    plot_mis_registration_table(system, ax=ax_table, colors=colors, fontsize=fontsize,
                                 **(table_kwargs or {}))
    fig.tight_layout()
    return fig, axes


_KIND_COLORS = {"wfs": "tab:blue", "dm": "tab:orange", "gs": "tab:green"}


def plot_mode_bar(specs, coefficients, ax=None, top_n=30, fontsize=8, title=None):
    """
    Horizontal bar chart of one mode's coefficients (e.g. one row of
    `Vt` from `analysis.svd_of_jacobian`, one value per entry of `specs`),
    sorted by absolute value and colour-coded by element kind (WFS/DM/GS).

    A complete, unfiltered view: `analysis.describe_mode`'s small-
    coefficient cutoff is meant for a quick read of the dominant terms,
    but it can hide contributions that are individually small yet
    collectively necessary for an exact degeneracy - see the module
    docstring.

    Parameters
    ----------
    specs : list of reconstruction.ParameterSpec
    coefficients : array-like, same length as `specs`
    ax : matplotlib.axes.Axes, optional
    top_n : int
        Show only the `top_n` largest-magnitude coefficients (default
        30); pass `None` to show all of `specs`.
    title : str, optional

    Returns
    -------
    ax : the matplotlib Axes used.
    """
    import matplotlib.pyplot as plt

    coefficients = np.asarray(coefficients, dtype=float)
    n_show = len(specs) if top_n is None else min(top_n, len(specs))
    order = np.argsort(-np.abs(coefficients))[:n_show]

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 0.28 * n_show + 1))

    y = np.arange(n_show)
    bar_colors = [_KIND_COLORS.get(specs[i].kind, "0.5") for i in order]
    ax.barh(y, coefficients[order], color=bar_colors)
    ax.set_yticks(y)
    ax.set_yticklabels([specs[i].get_label() for i in order], fontsize=fontsize)
    ax.invert_yaxis()
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_xlabel("coefficient")
    ax.set_title(title or f"Mode coefficients (top {n_show} of {len(specs)})")

    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in _KIND_COLORS.values()]
    ax.legend(handles, _KIND_COLORS.keys(), loc="lower right", fontsize=fontsize)
    return ax


def plot_mode_gs_quiver(system, specs, coefficients, ax=None, fontsize=9, title=None):
    """
    Sky-position quiver for the `("gs", name, "position_shift", 0/1)`
    entries of one mode: one arrow per guide star, at its actual sky
    position, in the (coefficient_x, coefficient_y) direction of this
    mode - a geometric view of the guide-star-position part of a
    degenerate/poorly-observed combination (see `plot_mode_bar` for the
    other element kinds and for a complete, non-geometric accounting).

    Guide stars with no `position_shift` entry in `specs` are skipped.

    Parameters
    ----------
    system : registration.model.System
    specs : list of reconstruction.ParameterSpec
    coefficients : array-like, same length as `specs`
    ax : matplotlib.axes.Axes, optional
    title : str, optional

    Returns
    -------
    ax : the matplotlib Axes used.
    """
    import matplotlib.pyplot as plt

    coefficients = np.asarray(coefficients, dtype=float)
    coeff_by_gs = {}
    for spec, c in zip(specs, coefficients):
        if spec.kind == "gs" and spec.field == "position_shift":
            coeff_by_gs.setdefault(spec.name, [0.0, 0.0])[spec.component] = c

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    xs, ys, us, vs, names = [], [], [], [], []
    for wfs in system.wfss.values():
        gs = wfs.guide_star
        if gs.name not in coeff_by_gs:
            continue
        cx, cy = coeff_by_gs[gs.name]
        xs.append(gs.position[0])
        ys.append(gs.position[1])
        us.append(cx)
        vs.append(cy)
        names.append(gs.name)

    if not xs:
        raise ValueError("None of `specs` is a ('gs', ..., 'position_shift', ...) entry.")

    ax.scatter(xs, ys, color="0.3", zorder=3)
    ax.quiver(xs, ys, us, vs, angles="xy", scale_units="xy", color="tab:green",
              width=0.008, zorder=2)
    for x, y, name in zip(xs, ys, names):
        ax.annotate(name, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=fontsize)
    ax.set_xlabel('x ["]')
    ax.set_ylabel('y ["]')
    ax.set_aspect("equal")
    ax.margins(0.3)
    ax.set_title(title or "Guide star position_shift coefficients (this mode)")
    return ax
