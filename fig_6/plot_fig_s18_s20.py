#!/usr/bin/env python3

import dataclasses
import functools
import pdb
import pickle
from typing import Callable, Self, Sequence

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
import sklearn.metrics

import pylib

matplotlib.rcParams["font.size"] = 8


@dataclasses.dataclass
class SeriesCfg(object):
	key: tuple
	display_name: str
	marker: str = "o"
	edgecolor_override: str | None = None
	facecolor_override: str | None = None
	text_xy: tuple[float, float] = (0.5, 0.5)
	horizontalalignment: str = "center"
	verticalalignment: str = "center"
	zorder: int = 0


@dataclasses.dataclass
class TransfPlotCfg(object):
	@dataclasses.dataclass
	class RefX(object):
		data: pandas.Series
		match_col: str
		displayname: str

	key: str
	file: str
	plot_basename: str
	series_cfgs: Sequence[SeriesCfg]
	ref_x: RefX | None = None

	@functools.cached_property
	def data(self) -> dict:
		with open(self.file, "rb") as fp:
			ret = pickle.load(fp)
		return ret


def setup_layout(clocks: Sequence[str]) -> dict:
	lc = mpllayout.LayoutCreator(
		left_margin=0.7,
		right_margin=0.2,
		bottom_margin=0.7,
		top_margin=0.3,
	)

	n_clocks = len(clocks)
	n_row = int(numpy.ceil(numpy.sqrt(n_clocks)))
	n_col = (n_clocks + n_row - 1) // n_row

	clock_w = 1.5
	clock_h = clock_w
	clock_gap_w = 0.2
	clock_gap_h = 0.4

	axes_pos = dict()
	for i, clock in enumerate(clocks):
		ri = n_row - 1 - i // n_col
		ci = i % n_col
		axes = lc.add_frame(f"clock_{clock}")
		axes.set_anchor("bottomleft", offsets=(
			clock_gap_w + (clock_w + clock_gap_w) * ci,
			clock_gap_h + (clock_h + clock_gap_h) * ri,
		))
		axes.set_size(clock_w, clock_h)
		axes_pos[clock] = (ri, ci)

	layout = lc.create_figure_layout()
	layout["n_row"] = n_row
	layout["n_col"] = n_col
	layout["axes_pos"] = axes_pos

	return layout


def _plot_transf_comp_clock_single(axes: matplotlib.axes.Axes, clock: str, *,
	plot_cfg: TransfPlotCfg, series_cfgs: Sequence[SeriesCfg],
	clock_cfgs: pylib.ClockCfgLib,
	show_xlabel: bool = False,
	show_ylabel: bool = False,
) -> tuple[float, float]:

	res_dict = plot_cfg.data["res"]
	clock_cfg = clock_cfgs[clock]
	mae_list = [0] * len(series_cfgs)
	r_list = [0] * len(series_cfgs)

	# plot per series
	for i, cfg in enumerate(series_cfgs):
		res = res_dict[cfg.key]

		# get x and y
		if plot_cfg.ref_x is not None:
			x = plot_cfg.ref_x.data[res[plot_cfg.ref_x.match_col]].values
		else:
			x = res["age"].values
		y = res[clock].values

		# plot scatter
		edgecolor = "#000000" if cfg.edgecolor_override is None \
			else cfg.edgecolor_override
		facecolor = clock_cfg.color if cfg.facecolor_override is None \
			else cfg.edgecolor_override

		axes.scatter(x, y, marker=cfg.marker, s=10, linewidths=0.5,
			edgecolors=edgecolor, facecolors=facecolor, zorder=cfg.zorder,
		)

		# add reg line
		slope, intercept, r_val, p_val, serr = scipy.stats.linregress(x, y)
		reg_x = numpy.linspace(-1000, 1000, 100)
		reg_y = intercept + slope * reg_x
		line = matplotlib.lines.Line2D(reg_x, reg_y, linestyle="-",
			linewidth=1.0, color=facecolor, zorder=15)
		axes.add_artist(line)

		# add text
		mae = sklearn.metrics.mean_absolute_error(x, y)
		text_lines = [
			f"[{cfg.display_name}]",
			f"R={r_val:.2f}",
			f"$p$={p_val:.2e}",
			f"MAE={mae:.1f} yrs",
		]
		text = "\n".join(text_lines)
		axes.text(*cfg.text_xy, text,
			horizontalalignment=cfg.horizontalalignment,
			verticalalignment=cfg.verticalalignment,
			transform=axes.transAxes,
			zorder=30,
		)

		# record
		mae_list[i] = mae
		r_list[i] = r_val

	# add diag line
	diag_x = numpy.linspace(-1000, 1000, 100)
	diag_y = diag_x
	line = matplotlib.lines.Line2D(diag_x, diag_y, linestyle="--",
		linewidth=1.0, color="#000000", zorder=15)
	axes.add_artist(line)

	# grid
	axes.grid(linewidth=0.5, color="#c0c0c0", zorder=0)

	# misc
	axes.tick_params(
		left=True, labelleft=show_ylabel,
		right=False, labelright=False,
		bottom=True, labelbottom=show_xlabel,
		top=False, labeltop=False,
	)
	axes.set_title(clock_cfg.key)
	if show_xlabel:
		xlabel = ("CA" if plot_cfg.ref_x is None else plot_cfg.ref_x.displayname) \
			+ " (years)"
		axes.set_xlabel(xlabel)
	if show_ylabel:
		axes.set_ylabel("Predicted age (years)")

	d_mae = mae_list[1] - mae_list[0]
	d_r = r_list[1] - r_list[0]

	return d_mae, d_r


def _plot_transf_comp_dataset(plot_cfg: TransfPlotCfg) -> None:
	# load data
	series_cfgs = plot_cfg.series_cfgs
	clocks = pylib.clock.load_clocks(exclude={"altumage"})
	clock_cfgs = pylib.ClockCfgLib.from_json()
	dataset_cfgs = pylib.DatasetCfgLib.from_json()

	# plot
	layout = setup_layout(clocks)
	figure = layout["figure"]
	axes_pos = layout["axes_pos"]

	ref_axes = None

	d_mae_list = list()
	d_r_list = list()

	for c in clocks:
		ri, ci = axes_pos[c]
		axes = layout[f"clock_{c}"]
		if ref_axes is None:
			ref_axes = axes
		else:
			axes.sharex(ref_axes)
			axes.sharey(ref_axes)

		mae, r_val = _plot_transf_comp_clock_single(axes, c,
			plot_cfg=plot_cfg,
			series_cfgs=series_cfgs,
			clock_cfgs=clock_cfgs,
			show_xlabel=(ri == 0),
			show_ylabel=(ci == 0),
		)
		d_mae_list.append(mae)
		d_r_list.append(r_val)

	# legend
	handles = list()
	for cfg in series_cfgs:
		edgecolor = "#000000" if cfg.edgecolor_override is None \
			else cfg.edgecolor_override
		s = figure.gca().scatter([], [], marker=cfg.marker, s=10, linewidths=0.5,
			edgecolors=edgecolor, facecolors="#c0c0c0", label=cfg.display_name,
		)
		handles.append(s)
	# legend
	legend = figure.legend(handles=handles, loc=3, bbox_to_anchor=(0.05, 0.01),
		frameon=False, fontsize=8, ncol=len(series_cfgs), handlelength=0.75,
		title="Improvement strategies", title_fontsize=8,
	)

	# figure title
	dataset_cfg = dataset_cfgs[plot_cfg.key]
	title = f"Training/testing: {dataset_cfg.display_name} (N={dataset_cfg.n_samples}), " \
		"5-fold cross-validation"
	figure.text(0.60, 0.02, title, fontsize=12,
		horizontalalignment="center", verticalalignment="bottom"
	)

	figure.savefig(plot_cfg.plot_basename + ".svg", dpi=600)
	figure.savefig(plot_cfg.plot_basename + ".png", dpi=600)
	figure.savefig(plot_cfg.plot_basename + ".pdf", dpi=600)
	matplotlib.pyplot.close(figure)

	# report stats
	for name, vals in zip(["MAE", "R"], [d_mae_list, d_r_list]):
		# print title, min, max, mean and median
		arr = numpy.asarray(vals)
		print(f"[{plot_cfg.key}] {name} improvement:")
		print(f"  min: {arr.min():.3f}")
		print(f"  max: {arr.max():.3f}")
		print(f"  mean: {arr.mean():.3f}")
		print(f"  median: {numpy.median(arr):.3f}")

	return


def _main():
	# plot config
	srrsh_ref_x = TransfPlotCfg.RefX(
		data=pandas.read_csv("srrsh_predicted_age_PhenoAge.txt", index_col=0,
			sep="\t")["Predicted_Age"],
		match_col="sample",
		displayname="PhenoAge",
	)

	plot_cfgs = [
		TransfPlotCfg(
			key="BUCCAL_TWIST",
			file="BUCCAL_TWIST_filter_imput_transfer.age_pred.pkl",
			plot_basename="fig_s18",
			series_cfgs=[
				SeriesCfg(key=(0, False, False), display_name="Direct", marker="o",
					edgecolor_override="#c0c0c0", facecolor_override="#c0c0c0",
					text_xy=(0.02, 0.98), horizontalalignment="left", verticalalignment="top",
					zorder=10,
				),
				SeriesCfg(key=(0, False, True), display_name="TL",
					marker="^", text_xy=(0.98, 0.02),
					horizontalalignment="right", verticalalignment="bottom",
					zorder=20,
				),
			],
		),
		TransfPlotCfg(
			key="RB_SRRSH",
			file="RB_SRRSH_filter_imput_transfer.age_pred.pkl",
			plot_basename="fig_s20",
			series_cfgs=[
				SeriesCfg(key=(0, False, False), display_name="Direct", marker="o",
					edgecolor_override="#c0c0c0", facecolor_override="#c0c0c0",
					text_xy=(0.02, 0.98), horizontalalignment="left", verticalalignment="top",
					zorder=10,
				),
				SeriesCfg(key=(0, True, True), display_name="IM+TL",
					marker="^", text_xy=(0.98, 0.02),
					horizontalalignment="right", verticalalignment="bottom",
					zorder=20,
				),
			],
			ref_x=srrsh_ref_x,
		),
	]

	for plot_cfg in plot_cfgs:
		_plot_transf_comp_dataset(plot_cfg)

	return


if __name__ == "__main__":
	_main()
