#!/usr/bin/env python3

import dataclasses
import functools
import io
import pdb
import pickle
from typing import Callable, Sequence

import matplotlib
import matplotlib.axes
import matplotlib.collections
import matplotlib.patches
import matplotlib.patheffects
import matplotlib.pyplot
import matplotlib.legend
import matplotlib.legend_handler
import matplotlib.lines
import mpllayout
import numpy
import pandas
import scipy
import scipy.stats
import sklearn.model_selection
import sklearn.metrics

import pylib

matplotlib.rcParams["font.size"] = 8


@dataclasses.dataclass
class TrajDataCfg(object):
	@dataclasses.dataclass
	class TrajTermCfg(object):
		key: tuple
		display_name: str

	@dataclasses.dataclass
	class RefX(object):
		data: pandas.Series
		match_col: str
		display_name: str

	key: str
	file: str
	traj_term_cfgs: Sequence[TrajTermCfg]
	ref_x: RefX | None = None

	@functools.cached_property
	def data(self) -> dict:
		with open(self.file, "rb") as fp:
			ret = pickle.load(fp)
		return ret


@dataclasses.dataclass
class DHCompCfg(object):
	key: str
	d_file: str
	h_file: str
	d_label: str = "Disease"
	h_label: str = "Healthy"
	plot_relative: bool = False

	@functools.cached_property
	def d_res(self) -> dict[str]:
		with open(self.d_file, "rb") as fp:
			ret = pickle.load(fp)
		return ret

	@functools.cached_property
	def h_res(self) -> dict[str]:
		with open(self.h_file, "rb") as fp:
			ret = pickle.load(fp)
		return ret


@dataclasses.dataclass
class DHCompStatCfg(object):
	key: str
	umeth: Callable
	display_name: str


class MultiPatchDummy(matplotlib.patches.Patch):
	def __init__(self, *ka, patches: Sequence[matplotlib.patches.Patch], **kw):
		self.patches = patches
		return super().__init__(*ka, **kw)


class MultiPatchHandler(matplotlib.legend_handler.HandlerPatch):
	def create_artists(self, legend: matplotlib.legend.Legend,
		orig_handle: MultiPatchDummy, xdescent: float, ydescent: float,
		width: float, height: float, fontsize: float,
		transform: matplotlib.transforms.Transform
	) -> list:
		patches = orig_handle.patches
		n_patches = len(patches)

		gap_frac = 0.20
		patch_w = width / (1 * n_patches + gap_frac * (n_patches - 1))
		gap_w = patch_w * gap_frac

		artists = list()
		for i, patch in enumerate(patches):
			x = (patch_w + gap_w) * i
			p = matplotlib.patches.Rectangle((x, 0),
				patch_w, height,
				linewidth=patch.get_linewidth(),
				edgecolor=patch.get_edgecolor(),
				facecolor=patch.get_facecolor(),
			)
			artists.append(p)
			l = matplotlib.lines.Line2D([x, x + patch_w], [0, height],
				color=patch.get_edgecolor(),
				linewidth=patch.get_linewidth(),
			)
			artists.append(l)

		return artists


class ClockMarkerHandler(matplotlib.legend_handler.HandlerPathCollection):
	def __init__(self, *ka, clock_cfgs: pylib.ClockCfgLib, **kw):
		self.clock_cfgs = clock_cfgs
		return super().__init__(*ka, **kw)

	def create_artists(self, legend, orig_handle: matplotlib.collections.PathCollection,
		xdescent: float, ydescent: float, width: float, height: float,
		fontsize: float, transform
	) -> list:
		artists = super().create_artists(legend, orig_handle,
			xdescent, ydescent, width, height, fontsize, transform)
		# overlay text
		s = self.clock_cfgs[orig_handle.get_label()].marker_char
		fontsize = 5 + 2 / len(s)
		t = matplotlib.text.Text(width / 2, height / 3, s,
			fontsize=fontsize, color="#000000",
			transform=transform, zorder=artists[0].zorder + 1,
			path_effects=[
				matplotlib.patheffects.Stroke(linewidth=1.5, foreground="#ffffff"),
			],
			horizontalalignment="center", verticalalignment="center",
		)
		artists.append(t)
		t = matplotlib.text.Text(width / 2, height / 3, s,
			fontsize=fontsize, color="#000000",
			transform=transform, zorder=artists[0].zorder + 1,
			horizontalalignment="center", verticalalignment="center",
		)
		artists.append(t)
		return artists


def setup_layout(
	traj_data_cfgs: Sequence[TrajDataCfg],
	dh_comp_cfgs: Sequence[DHCompCfg],
	clocks: Sequence[str],
) -> dict:
	lc = mpllayout.LayoutCreator(
		left_margin=0.7,
		right_margin=0.4,
		bottom_margin=0.7,
		top_margin=0.5,
	)

	# disease/health comp
	n_dh_comps = len(dh_comp_cfgs)
	n_clocks = len(clocks)
	n_col = 2
	n_row = (n_clocks + n_col - 1) // n_col

	dh_col_w = 4.0
	dh_h = 5.2
	dh_gap_w_rela = 0.1

	dh_w = dh_col_w / (n_col + dh_gap_w_rela * (n_col - 1))
	dh_gap_w = dh_w * dh_gap_w_rela
	dh_gap_h = 0.3
	dh_col_gap_w = 1.0

	for i, comp_cfg in enumerate(dh_comp_cfgs):
		for j in range(n_col):
			dh = lc.add_frame(f"dh_{comp_cfg.key}_{j}")
			dh.set_anchor("bottomleft", offsets=(
				(dh_col_w + dh_col_gap_w) * i + (dh_w + dh_gap_w) * j,
				0.0,
			))
			dh.set_size(dh_w, dh_h)

	# performance comp
	perf_w = 2.8
	perf_h = 2.2
	perf_gap_w = 1.5

	perf_offset_w = dh_col_w * n_dh_comps + dh_col_gap_w * (n_dh_comps - 1) \
		+ perf_gap_w
	perf = lc.add_frame("perf")
	perf.set_anchor("bottomleft", offsets=(perf_offset_w, 0.0))
	perf.set_size(perf_w, perf_h)

	# tl entropy comp layout
	te_w = perf_w
	te_scan_h = 1.2
	te_abs_h = 0.6
	te_gap_h = 0.6

	te_scan = lc.add_frame("te_scan")
	te_scan.set_anchor("bottomleft", perf, "topleft", offsets=(0.0, te_gap_h))
	te_scan.set_size(te_w, te_scan_h)

	te_abs = lc.add_frame("te_abs")
	te_abs.set_anchor("bottomleft", te_scan, "topleft", offsets=(0.0, te_gap_h))
	te_abs.set_size(te_w, te_abs_h)

	# im+tl trajectory
	traj_sec_gap_h = 0.4
	traj_offset_h = dh_h + dh_gap_h + traj_sec_gap_h
	traj_w = 3.2
	traj_h = traj_w
	traj_gap_w = 0.7

	for i, cfg in enumerate(traj_data_cfgs):
		traj = lc.add_frame(f"traj_{cfg.key}")
		traj.set_anchor("bottomleft",
			offsets=((traj_w + traj_gap_w) * i, traj_offset_h))
		traj.set_size(traj_w, traj_h)

	layout = lc.create_figure_layout()
	layout["n_row"] = n_row
	layout["n_col"] = n_col

	return layout


def _plot_traj_single(axes: matplotlib.axes.Axes, data_cfg: TrajDataCfg,
	*, dataset_cfgs: pylib.DatasetCfgLib,
	clocks: Sequence[str], clock_cfgs: pylib.clock_cfg.ClockCfgLib,
	add_legend: bool = False,
	stat_fp: io.TextIOBase | None = None,
) -> None:
	res_dict = data_cfg.data["res"]

	# extract data
	true_vals = list()
	for cfg in data_cfg.traj_term_cfgs:
		_res = res_dict[cfg.key]
		if data_cfg.ref_x is not None:
			ref_x = data_cfg.ref_x
			true_vals.append(ref_x.data.loc[_res[ref_x.match_col]])
			ref_label = ref_x.display_name
		else:
			true_vals.append(_res["age"])
			ref_label = "CA"

	# calculate and plot per clock
	delta_mae = list()
	delta_r = list()

	clock_handles = list()
	for c in clocks:
		pred_vals = list()
		for cfg in data_cfg.traj_term_cfgs:
			_res: pylib.AgePredRes = res_dict[cfg.key]
			pred_vals.append(_res[c])

		# calculate stats for plot
		# x is mae
		x0 = sklearn.metrics.mean_absolute_error(true_vals[0], pred_vals[0])
		x1 = sklearn.metrics.mean_absolute_error(true_vals[1], pred_vals[1])
		# y is corr r
		y0, *_ = scipy.stats.pearsonr(true_vals[0], pred_vals[0])
		y1, *_ = scipy.stats.pearsonr(true_vals[1], pred_vals[1])

		# record for min/max calculations
		delta_mae.append(x0 - x1)  # mae decrease
		delta_r.append(y1 - y0)  # r increase

		# plot scatter
		clock_cfg: pylib.ClockCfg = clock_cfgs[c]
		p = axes.scatter([x0], [y0], marker="o", s=100, linewidths=0.5,
			edgecolors="#000000", facecolors=clock_cfg.color,
			label=clock_cfg.key,
			zorder=15,
		)
		clock_handles.append(p)
		p = axes.scatter([x1], [y1], marker="^", s=100, linewidths=0.5,
			edgecolors="#000000", facecolors=clock_cfg.color,
			label=clock_cfg.key,
			zorder=12,
		)
		# add clock text only for primary marker (legend)
		mch = clock_cfg.marker_char
		fontsize = 5 + 2 / len(mch)
		axes.text(x0, y0, mch, fontsize=fontsize, zorder=20,
			path_effects=[
				matplotlib.patheffects.Stroke(linewidth=1.5,
				foreground="#ffffff"),
			],
			horizontalalignment="center", verticalalignment="center",
		)
		axes.text(x0, y0, mch, fontsize=fontsize, zorder=25,
			horizontalalignment="center", verticalalignment="center",
		)
		# add annotation
		anno_xy = numpy.asarray([(x0 + x1) / 2, (y0 + y1) / 2])
		anno_fullvec = numpy.asarray([x1 - x0, y1 - y0])
		anno_vec = 0.1 * anno_fullvec / numpy.linalg.norm(anno_fullvec)
		anno_xytext = anno_xy - anno_vec
		axes.annotate("", xy=anno_xy, xytext=anno_xytext, zorder=10,
			xycoords=axes.transData, textcoords=axes.transData,
			arrowprops=dict(arrowstyle="fancy", color=clock_cfgs[c].color,
				shrinkA=6, shrinkB=6, linewidth=0.6),
		)
		axes.plot([x0, x1], [y0, y1], linestyle="-", linewidth=0.5,
			color=clock_cfgs[c].color, zorder=10,
		)

	# always add imput legend
	imput_handles = [
		axes.scatter([], [], marker="o", s=80, linewidths=0.5,
			edgecolors="#000000", facecolors="#e0e0e0",
			label=data_cfg.traj_term_cfgs[0].display_name,
		),
		axes.scatter([], [], marker="^", s=80, linewidths=0.5,
			edgecolors="#000000", facecolors="#e0e0e0",
			label=data_cfg.traj_term_cfgs[1].display_name,
		),
	]
	legend = axes.legend(handles=imput_handles, ncol=2,
		loc=1, bbox_to_anchor=[1.00, 1.00], frameon=True, fontsize=8,
		handlelength=0.75, title="Method", title_fontsize=8, borderpad=0.5,
	)
	axes.add_artist(legend)

	if add_legend:
		# clock legend
		handler_map = matplotlib.legend.Legend.get_default_handler_map().copy()
		handler_map[matplotlib.collections.PathCollection] = ClockMarkerHandler(
			clock_cfgs=clock_cfgs
		)
		legend = axes.legend(handles=clock_handles, handler_map=handler_map,
			loc=2, bbox_to_anchor=[1.02, 1.10], frameon=False, fontsize=8,
			handlelength=0.75, title="Clocks (across a-b)", title_fontsize=10,
		)

	# misc
	axes.grid(linestyle="-", linewidth=0.5, color="#e0e0e0", zorder=0)
	axes.set_xlabel(f"MAE (BA~{ref_label}, years)", fontsize=10)
	axes.set_ylabel(f"Pearson's R (BA~{ref_label})", fontsize=10)
	dataset_cfg = dataset_cfgs[data_cfg.key]
	axes.set_title(f"Training/testing: {dataset_cfg.display_name} (N={dataset_cfg.n_samples})\n5-fold cross-validation",
		fontsize=10)

	# report min/max
	if stat_fp:
		for s, vals in zip(["MAE", "R"], [delta_mae, delta_r]):
			print(f"delta {s} for traj {data_cfg.key}", file=stat_fp)
			print(f"  min: {numpy.min(vals):.4f}", file=stat_fp)
			print(f"  max: {numpy.max(vals):.4f}", file=stat_fp)
			print(f"  mean: {numpy.mean(vals):.4f}", file=stat_fp)
			print(f"  median: {numpy.median(vals):.4f}", file=stat_fp)

	return


def _plot_traj_by_dataset(layout: dict[str], data_cfgs: Sequence[TrajDataCfg], *,
	dataset_cfgs: pylib.DatasetCfgLib = None, clocks: Sequence[str] = None,
	clock_cfgs: pylib.clock_cfg.ClockCfgLib = None,
	stat_output: str | None = None,
	skip: bool = False,
) -> None:
	if skip:
		pylib.logger.info("skipping plotting best strat dataset")
		return

	# load data
	if clocks is None:
		clocks = pylib.clock.load_clocks(exclude={"altumage"})
	if dataset_cfgs is None:
		dataset_cfgs = pylib.DatasetCfgLib.from_json()
	if clock_cfgs is None:
		clock_cfgs = pylib.clock_cfg.ClockCfgLib.from_json()

	if stat_output is not None:
		stat_fp = open(stat_output, "w")
	else:
		stat_fp = None

	for i, data_cfg in enumerate(data_cfgs):
		_plot_traj_single(layout[f"traj_{data_cfg.key}"], data_cfg,
			dataset_cfgs=dataset_cfgs,
			clocks=clocks, clock_cfgs=clock_cfgs,
			add_legend=(i == len(data_cfgs) - 1),
			stat_fp=stat_fp,
		)

	if stat_fp is not None:
		stat_fp.close()

	return


def _calc_delta_age(pred_res: dict[str], clocks: Sequence[str],
) -> tuple[pandas.DataFrame, pandas.DataFrame]:
	age: pandas.Series = pred_res["base_pred"]["age"]
	mask = age.notna()  # drop na values

	da_base = pred_res["base_pred"].loc[mask, clocks] - \
		age.loc[mask].values.reshape(-1, 1)
	da_cross = pred_res["cross_pred"].loc[mask, clocks] - \
		age.loc[mask].values.reshape(-1, 1)
	return da_base, da_cross


def _calc_dist_jsd(v1: numpy.ndarray, v2: numpy.ndarray) -> float:
	# kde to calculate density
	vmin = min(v1.min(), v2.min())
	vmax = max(v1.max(), v2.max())
	pos = numpy.linspace(vmin, vmax, 1000)
	kde1 = scipy.stats.gaussian_kde(v1)
	kde2 = scipy.stats.gaussian_kde(v2)
	p1 = kde1(pos)
	p2 = kde2(pos)
	# normalize to prob distribution
	p1 /= p1.sum()
	p2 /= p2.sum()
	# calculate jsd
	jsd = scipy.spatial.distance.jensenshannon(p1, p2)
	return jsd


def _plot_base_cross_comp(axes: matplotlib.axes.Axes, *,
	base_vals: numpy.ndarray, cross_vals: numpy.ndarray,
	pos: float,
	val_offset: float = 0, pos_shift: float = 0.02, width: float = 0.36,
	color: str = "#808080", zorder_base: int = 10,
) -> None:
	for d, p, s, ec, fc in zip(
		[base_vals, cross_vals],
		[pos - pos_shift, pos + pos_shift],
		["low", "high"],
		[color, color],
		["#ffffff", color],
	):
		violin = axes.violinplot(d - val_offset, [p], orientation="horizontal", side=s,
			widths=width, showmeans=False, showmedians=True, showextrema=False,
		)
		for b in violin["bodies"]:
			b.set(linewidth=0.5, edgecolor=ec, facecolor=fc + "40", alpha=None,
				zorder=zorder_base,
		)
		for part in ["cbars", "cmins", "cmaxes", "cmedians"]:
			if part in violin:
				violin[part].set(linewidth=0.5, color=ec, zorder=zorder_base + 1)

	return


def _plot_dh_comp_single_column(axes: matplotlib.axes.Axes, *,
	h_da_base: pandas.DataFrame,
	h_da_cross: pandas.DataFrame,
	d_da_base: pandas.DataFrame,
	d_da_cross: pandas.DataFrame,
	clocks: Sequence[str],
	clock_order_values: Sequence[float],
	clock_cfgs: pylib.ClockCfgLib,
	h_color: str,
	d_color: str,
	n_row: int,
	plot_relative: bool = False,
) -> None:

	# plot by clock
	for i, c in enumerate(clocks):
		# val_offset = numpy.median(h_da_base[c])
		if plot_relative:
			val_offset = numpy.median(h_da_base[c])
			xlabel = "Delta age shift (years)"
		else:
			val_offset = 0.0
			xlabel = "Delta age (years)"

		# healthy
		_plot_base_cross_comp(axes,
			base_vals=h_da_base[c].values,
			cross_vals=h_da_cross[c].values,
			val_offset=val_offset,
			pos=n_row - i - 0.80,
			color=h_color,
			zorder_base=10,
		)
		# diseased
		_plot_base_cross_comp(axes,
			base_vals=d_da_base[c].values,
			cross_vals=d_da_cross[c].values,
			val_offset=val_offset,
			pos=n_row - i - 0.45,
			color=d_color,
			zorder_base=10,
		)

		# add clock text
		clock_cfg: pylib.ClockCfg = clock_cfgs[c]
		axes.text(0.02, (n_row - i) / n_row - 0.005, clock_cfg.key,
			fontsize=8, color="#202020", zorder=20,
			horizontalalignment="left", verticalalignment="top",
			transform=axes.transAxes,
		)
		# continue
		axes.text(0.98, (n_row - i) / n_row - 0.005, "ΔJSD={:.3f}".format(clock_order_values[i]),
			fontsize=6, color="#202020", zorder=20,
			horizontalalignment="right", verticalalignment="top",
			transform=axes.transAxes,
		)

	# grid
	axes.grid(axis="x", linestyle="-", linewidth=0.5, color="#e0e0e0", zorder=0)
	for i in range(1, n_row):
		axes.axhline(i, linestyle="-", linewidth=0.5, color="#c0c0c0", zorder=0)
	# axes.axvline(0, linestyle="-", linewidth=1.0, color="#808080", zorder=1)

	# misc
	axes.tick_params(
		left=False, labelleft=False,
		right=False, labelright=False,
		bottom=True, labelbottom=True,
		top=False, labeltop=False,
	)
	axes.set_ylim(0, n_row)
	axes.set_xlabel(xlabel, fontsize=10)

	return


def _plot_dh_comp_single(layout: dict, comp_cfg: DHCompCfg, *,
	dataset_cfgs: pylib.DatasetCfgLib, clocks: Sequence[str],
	clock_cfgs: pylib.clock_cfg.ClockCfgLib,
	dist_meth: Callable = _calc_dist_jsd,
	show_legend: bool = False,
) -> None:
	h_color = "#4076ff"
	d_color = "#ff4040"
	n_row = layout["n_row"]
	n_col = layout["n_col"]

	# calc delta age
	d_da_base, d_da_cross = _calc_delta_age(comp_cfg.d_res, clocks)
	h_da_base, h_da_cross = _calc_delta_age(comp_cfg.h_res, clocks)

	# sort clocks by delta dist
	dist_diff = pandas.Series(index=clocks, dtype=float)
	for c in clocks:
		dist_base = dist_meth(
			d_da_base[c].values,
			h_da_base[c].values,
		)
		dist_cross = dist_meth(
			d_da_cross[c].values,
			h_da_cross[c].values,
		)
		dist_diff[c] = dist_cross - dist_base
	dist_diff.sort_values(ascending=False, inplace=True)

	# plot
	layout[f"dh_{comp_cfg.key}_0"].sharex(layout[f"dh_{comp_cfg.key}_1"])

	# plot by column
	for col in range(n_col):
		axes = layout[f"dh_{comp_cfg.key}_{col}"]
		clock_ordering = dist_diff[col * n_row:(col + 1) * n_row]
		_plot_dh_comp_single_column(axes,
			h_da_base=h_da_base,
			h_da_cross=h_da_cross,
			d_da_base=d_da_base,
			d_da_cross=d_da_cross,
			clocks=clock_ordering.index.tolist(),
			clock_order_values=clock_ordering.values.tolist(),
			clock_cfgs=clock_cfgs,
			h_color=h_color,
			d_color=d_color,
			n_row=n_row,
			plot_relative=comp_cfg.plot_relative,
		)

	# legend
	if show_legend:
		handles = [
			matplotlib.patches.Patch(linewidth=0.5, edgecolor=d_color,
				facecolor=d_color + "40", label=f"{comp_cfg.d_label}\n(IM+TL)"),
			matplotlib.patches.Patch(linewidth=0.5, edgecolor=d_color,
				facecolor="#ffffff", label=f"{comp_cfg.d_label}\n(Direct)"),
			matplotlib.patches.Patch(linewidth=0.5, edgecolor=h_color,
				facecolor=h_color + "40", label=f"{comp_cfg.h_label}\n(IM+TL)"),
			matplotlib.patches.Patch(linewidth=0.5, edgecolor=h_color,
				facecolor="#ffffff", label=f"{comp_cfg.h_label}\n(Direct)"),
		]
		handler_map = matplotlib.legend.Legend.get_default_handler_map().copy()
		handler_map[MultiPatchDummy] = MultiPatchHandler()
		legend = axes.legend(handles=handles, handler_map=handler_map,
			loc=2, bbox_to_anchor=[1.02, 1.02], frameon=False, fontsize=8,
			handlelength=0.75,
		)

	# add dataset text
	axes = layout[f"dh_{comp_cfg.key}_1"]
	lines: list[str] = list()
	# training
	lines.append("Training:")
	for train_key in comp_cfg.d_res["train_list"]:
		train_cfg: pylib.DatasetCfg = dataset_cfgs[train_key]
		lines.append(f"  {train_cfg.display_name}")
	# testing
	lines.append("")
	lines.append("Testing:")
	test_cfg: pylib.DatasetCfg = dataset_cfgs[comp_cfg.key]
	lines.append(f"  {test_cfg.display_name}")
	axes.text(1.02, 0.00, "\n".join(lines),
		fontsize=6, color="#202020",
		horizontalalignment="left", verticalalignment="bottom",
		transform=axes.transAxes,
	)

	return


def _plot_dh_comp(layout: dict[str], comp_cfgs: Sequence[DHCompCfg],
	dataset_cfgs: pylib.DatasetCfgLib = None, clocks: Sequence[str] = None,
	clock_cfgs: pylib.clock_cfg.ClockCfgLib = None,
	dist_meth: Callable = _calc_dist_jsd,
	skip: bool = False,
) -> None:
	if skip:
		pylib.logger.info("skipping plotting dh comp")
		return

	# config
	if clocks is None:
		clocks = pylib.clock.load_clocks(exclude={"altumage"})
	if clock_cfgs is None:
		clock_cfgs = pylib.clock_cfg.ClockCfgLib.from_json()
	if dataset_cfgs is None:
		dataset_cfgs = pylib.DatasetCfgLib.from_json()

	for comp_cfg in comp_cfgs:
		_plot_dh_comp_single(layout, comp_cfg,
			dataset_cfgs=dataset_cfgs,
			clocks=clocks,
			clock_cfgs=clock_cfgs,
			dist_meth=dist_meth,
			show_legend=True,
		)

	return


def _plot_box_single(axes: matplotlib.axes.Axes, values: numpy.ndarray, *,
	pos: float, width: float = 0.5, edgecolor: str = "#000000",
	facecolor: str = None, zorder_ref: int = 10,
) -> None:
	if facecolor is None:
		facecolor = edgecolor + "40"

	# box
	box = axes.boxplot(values, positions=[pos],
		widths=width, orientation="horizontal", patch_artist=True,
		showfliers=False, capwidths=width * 0.6,  # whis=(0, 100),
		boxprops=dict(linewidth=0.5, edgecolor=edgecolor, facecolor=facecolor),
		medianprops=dict(linewidth=0.5, color=edgecolor),
		whiskerprops=dict(linewidth=0.5, color=edgecolor),
		capprops=dict(linewidth=0.5, color=edgecolor),
		zorder=zorder_ref,
	)
	# scatter
	axes.scatter(values, [pos] * len(values), marker="o", s=8, linewidths=0.5,
		edgecolors=edgecolor, facecolors="none", zorder=zorder_ref + 10,
	)

	return


def _plot_te_abs(layout: dict[str], comp_cfg: DHCompCfg, *,
	dataset_cfgs: pylib.DatasetCfgLib, clocks: Sequence[str],
	dist_meth: Callable = _calc_dist_jsd,
) -> None:
	clock_cfgs = pylib.clock_cfg.ClockCfgLib.from_json()

	# calc delta age
	d_da_base, d_da_cross = _calc_delta_age(comp_cfg.d_res, clocks)
	h_da_base, h_da_cross = _calc_delta_age(comp_cfg.h_res, clocks)

	# calculate dist diff
	dist_base = pandas.Series(index=clocks, dtype=float)
	dist_cross = pandas.Series(index=clocks, dtype=float)
	for c in clocks:
		dist_base[c] = dist_meth(
			d_da_base[c].values,
			h_da_base[c].values,
		)
		dist_cross[c] = dist_meth(
			d_da_cross[c].values,
			h_da_cross[c].values,
		)
	dist_diff = dist_cross - dist_base

	# plot dist box
	te_abs = layout["te_abs"]
	_plot_box_single(te_abs, dist_base.values, pos=0.25, width=0.4,
		edgecolor="#808080", zorder_ref=10,
	)
	_plot_box_single(te_abs, dist_cross.values, pos=0.75, width=0.4,
		edgecolor="#ffa000", zorder_ref=10,
	)

	# plot lines
	for c in clocks:
		te_abs.plot(
			[dist_base[c], dist_cross[c]],
			[0.25, 0.75],
			linestyle="-",
			linewidth=0.5,
			color="#80808040",
			zorder=5,
		)

	# misc
	te_abs.set_ylim(0.0, 1.0)
	te_abs.set_xlabel("JSD (ALS vs. Normal delta ages, AU)",
		fontsize=10,
	)
	te_abs.set_yticks((0.25, 0.75))
	te_abs.set_yticklabels(("Direct", "IM+TL"))

	return


def _calc_classif_res(d_da: numpy.ndarray, h_da: numpy.ndarray,
) -> tuple[float, float, float]:
	# create matrix
	x = numpy.vstack([d_da, h_da])
	y = numpy.asarray([1] * len(d_da) + [0] * len(h_da), dtype=int)
	pred_y = numpy.zeros_like(y, dtype=int)

	# model
	driver = sklearn.svm.SVC(kernel="linear", probability=True)

	# cross-validation
	# cv = sklearn.model_selection.LeaveOneOut()
	cv = sklearn.model_selection.StratifiedKFold(12, shuffle=True,
		random_state=42,
	)
	for train_idx, test_idx in cv.split(x, y):
		driver.fit(x[train_idx], y[train_idx])
		# pred_y[test_idx] = (driver.predict_proba(x[test_idx])[:, 1] > 0.5).astype(int)
		pred_y[test_idx] = driver.predict(x[test_idx])

	# calculate auc
	auc = sklearn.metrics.roc_auc_score(y, pred_y)
	if auc > 0.5:
		sens = sklearn.metrics.recall_score(y, pred_y, pos_label=1)
		spec = sklearn.metrics.recall_score(y, pred_y, pos_label=0)
	else:
		# in a binary classification, if AUC < 0.5, simply invert the prediction
		auc = 1.0 - auc
		sens = sklearn.metrics.recall_score(y, 1 - pred_y, pos_label=0)
		spec = sklearn.metrics.recall_score(y, 1 - pred_y, pos_label=1)

	return auc, sens, spec


def _plot_te_optim_top_clocks(layout: dict[str], comp_cfg: DHCompCfg, *,
	dataset_cfgs: pylib.DatasetCfgLib, clocks: Sequence[str],
	dist_meth: Callable = _calc_dist_jsd,
) -> Sequence[str]:
	# calc delta age
	d_da_base, d_da_cross = _calc_delta_age(comp_cfg.d_res, clocks)
	h_da_base, h_da_cross = _calc_delta_age(comp_cfg.h_res, clocks)

	# calculate dist diff
	dist_base = pandas.Series(index=clocks, dtype=float)
	dist_cross = pandas.Series(index=clocks, dtype=float)
	for c in clocks:
		dist_base[c] = dist_meth(
			d_da_base[c].values,
			h_da_base[c].values,
		)
		dist_cross[c] = dist_meth(
			d_da_cross[c].values,
			h_da_cross[c].values,
		)
	dist_diff = dist_cross - dist_base
	dist_diff.sort_values(ascending=False, inplace=True)

	# select top clocks
	res_list = list()
	for n in range(len(dist_diff)):
		top_clocks = dist_diff.index.tolist()[:n + 1]
		if dist_diff.values[n] < 0:
			break

		# calc auc
		res = _calc_classif_res(d_da_cross[top_clocks], h_da_cross[top_clocks])
		res_list.append((n, *res))

	# plot
	te_scan = layout["te_scan"]
	x = [v[0] + 1 for v in res_list]
	y = [v[1] for v in res_list]
	te_scan.bar(x, y, width=0.95, align="center", linewidth=0.5,
		edgecolor="#228DFF", facecolor="#228DFF40", zorder=10,
	)

	# misc
	te_scan.set_xlim(0.5, len(res_list) + 0.5)
	te_scan.set_ylim(0.0, None)
	te_scan.set_xlabel("Top 1-n (by ΔJSD) clocks used for classification",
		fontsize=10)
	te_scan.set_ylabel("AUC (ALS vs. Normal)", fontsize=10)
	te_scan.set_xticks(x)

	# find best
	sorted_res = sorted(res_list, key=lambda v: v[1], reverse=True)
	best_n = sorted_res[0][0]
	top_clocks = dist_diff.index.tolist()[:best_n + 1]

	return top_clocks


def _plot_te_perf_single(axes: matplotlib.axes.Axes,
	d_values: numpy.ndarray, h_values: numpy.ndarray, *,
	color: str = "#808080", zorder_ref: int = 10,
) -> None:

	auc, sens, spec = _calc_classif_res(d_values, h_values)

	# plot ROC curve
	axes.plot([0, 1 - spec, 1], [0, sens, 1], linestyle="-", linewidth=1.0,
		color=color, zorder=zorder_ref + 5,
	)
	axes.plot([0, 1], [0, 1], linestyle="--", linewidth=0.5, color="#c0c0c0",
		zorder=zorder_ref,
	)
	# add text
	pathfx = matplotlib.patheffects.Stroke(linewidth=1.5, foreground="#ffffff")
	axes.text(1 - spec, sens, f" AUC={auc:.2f}",
		fontsize=10, color="#ffffff", zorder=zorder_ref + 10,
		path_effects=[pathfx],
		horizontalalignment="left", verticalalignment="center",
	)
	axes.text(1 - spec, sens, f" AUC={auc:.3f}",
		fontsize=10, color=color, zorder=zorder_ref + 10,
		horizontalalignment="left", verticalalignment="center",
	)

	return


def _plot_te_perf(layout: dict[str], comp_cfg: DHCompCfg, *,
	dataset_cfgs: pylib.DatasetCfgLib, clocks: Sequence[str],
) -> None:
	# calc delta age
	d_da_base, d_da_cross = _calc_delta_age(comp_cfg.d_res, clocks)
	h_da_base, h_da_cross = _calc_delta_age(comp_cfg.h_res, clocks)

	perf = layout["perf"]
	_plot_te_perf_single(perf,
		d_da_base[clocks].values,
		h_da_base[clocks].values,
		color="#808080",
		zorder_ref=10,
	)
	_plot_te_perf_single(perf,
		d_da_cross[clocks].values,
		h_da_cross[clocks].values,
		color="#ffa000",
		zorder_ref=15,
	)

	# legends
	handles = [
		matplotlib.lines.Line2D([0], [0], linestyle="-", linewidth=1.0,
			color="#808080", label="Direct"),
		matplotlib.lines.Line2D([0], [0], linestyle="-", linewidth=1.0,
			color="#ffa000", label="IM+TL"),
	]
	legend = perf.legend(handles=handles, loc=4, bbox_to_anchor=[0.98, 0.02],
		frameon=True, fontsize=8, handlelength=1.2,
	)

	# misc
	perf.set_xlim(0.0, 1.0)
	perf.set_ylim(0.0, 1.0)
	perf.set_xlabel("1 - Specificity", fontsize=10)
	perf.set_ylabel("Sensitivity", fontsize=10)

	return


def _main():
	# configs
	srrsh_ref_x = TrajDataCfg.RefX(
		data=pandas.read_csv("srrsh_predicted_age_PhenoAge.txt", index_col=0,
			sep="\t")["Predicted_Age"],
		match_col="sample",
		display_name="PhenoAge",
	)

	traj_data_cfgs = [
		TrajDataCfg(
			key="BUCCAL_TWIST",
			file="BUCCAL_TWIST_filter_imput_transfer.age_pred.pkl",
			traj_term_cfgs=[
				TrajDataCfg.TrajTermCfg(key=(0, False, False), display_name="Direct"),
				TrajDataCfg.TrajTermCfg(key=(0, False, True), display_name="TL"),
			],
		),
		TrajDataCfg(
			key="RB_SRRSH",
			file="RB_SRRSH_filter_imput_transfer.age_pred.pkl",
			traj_term_cfgs=[
				TrajDataCfg.TrajTermCfg(key=(0, False, False), display_name="Direct"),
				TrajDataCfg.TrajTermCfg(key=(0, True, True), display_name="IM+TL"),
			],
			ref_x=srrsh_ref_x,
		),
	]
	dh_comp_cfgs = [
		DHCompCfg(
			key="GSE164600",
			d_file="cross_test.RB_GALAXY-RB_TWIST-RB_SRRSH.GSE164600_D.no_df.with_imput.age_pred.pkl",
			h_file="cross_test.RB_GALAXY-RB_TWIST-RB_SRRSH.GSE164600_N.no_df.with_imput.age_pred.pkl",
			d_label="ALS",
			h_label="Normal",
			plot_relative=True,
		),
	]

	# load data
	dataset_cfgs = pylib.DatasetCfgLib.from_json()
	clocks = pylib.clock.load_clocks(exclude={"altumage"})

	# plot
	layout = setup_layout(traj_data_cfgs, dh_comp_cfgs, clocks)

	figure = layout["figure"]

	_plot_traj_by_dataset(layout, traj_data_cfgs,
		stat_output="tl_improve_stats.txt",
		skip=False,
	)

	_plot_dh_comp(layout, dh_comp_cfgs,
		skip=False,
	)

	_plot_te_abs(layout, dh_comp_cfgs[0],
		dataset_cfgs=dataset_cfgs,
		clocks=clocks,
		dist_meth=_calc_dist_jsd,
	)

	top_clocks = _plot_te_optim_top_clocks(layout, dh_comp_cfgs[0],
		dataset_cfgs=dataset_cfgs,
		clocks=clocks,
		dist_meth=_calc_dist_jsd,
	)

	_plot_te_perf(layout, dh_comp_cfgs[0],
		dataset_cfgs=dataset_cfgs,
		clocks=top_clocks,
	)

	figure.savefig("fig_6.svg", dpi=600)
	figure.savefig("fig_6.png", dpi=600)
	figure.savefig("fig_6.pdf", dpi=600)
	matplotlib.pyplot.close(figure)
	return


if __name__ == "__main__":
	_main()
