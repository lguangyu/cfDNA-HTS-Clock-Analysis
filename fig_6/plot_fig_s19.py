#!/usr/bin/env python3

import argparse
import collections
import dataclasses
import pickle
import pdb
import types
from typing import Callable, NoReturn, Self, Sequence

import matplotlib
import matplotlib.patches
import matplotlib.pyplot
import matplotlib.axes
import mpllayout
import numpy
import pandas

import pylib


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("dataset", type=str, help="dataset name")
	ap.add_argument("--with-rd", action="store_true",
		help="set to RD metric [no]"
	)
	ap.add_argument("-p", "--plot-basename", type=str, required=True,
		metavar="basename",
		help="base name for output plot files")

	args = ap.parse_args()
	return args


class KeyCfgLib(set):
	@dataclasses.dataclass
	class KeyCfg(object):
		depth_filtering: int | None = None
		imputation: bool = False
		transfer_learning: bool = False

		def __post_init__(self):
			if (self.depth_filtering is None) or (self.depth_filtering < 0):
				self.depth_filtering = 0
			return

		def to_tuple(self) -> tuple:
			return (self.depth_filtering, self.imputation, self.transfer_learning)

		def __hash__(self) -> int:
			return hash(self.to_tuple())

		@property
		def order_key(self) -> int:
			ret = -int(self.transfer_learning) * 100 \
				- int(self.imputation) * 10 \
				- bool(self.depth_filtering)
			return ret

	def add_new(self, **kw) -> KeyCfg:
		new = self.KeyCfg(**kw)
		self.add(new)
		return new

	@property
	def key_field_values(self) -> dict[str, set]:
		ret = collections.defaultdict(set)
		for cfg in self:
			for field in dataclasses.fields(cfg):
				name = field.name
				ret[name].add(getattr(cfg, name))
		return ret


class KeyFieldList(list["KeyFieldList.KeyFieldCfg"]):
	@dataclasses.dataclass
	class KeyFieldCfg(object):
		key: str
		full_name: str
		short_name: str
		color: str
		# return if the field is 'active' given a field value
		is_active_in: Callable[[Self, KeyCfgLib.KeyCfg], bool] = \
			lambda self, v: bool(getattr(v, self.key))
		# return the representative label given a field value
		get_active_label: Callable[[Self, KeyCfgLib.KeyCfg], str | None] = \
			lambda self, v: self.short_name

		def __post_init__(self):
			self.is_active_in = types.MethodType(self.is_active_in, self)
			self.get_active_label = types.MethodType(self.get_active_label, self)
			return

		def get_label(self, v: KeyCfgLib.KeyCfg) -> str | None:
			return self.get_active_label(v) if self.is_active_in(v) else None

	def get_field_labels(self, cfg: KeyCfgLib.KeyCfg) -> list[str]:
		ret = list()
		for field in self:
			if (v := field.get_label(cfg)) is not None:
				ret.append(v)
		return ret


@dataclasses.dataclass
class MetricCfg(object):
	key: str
	display_name: str
	unit: str = None

	@property
	def label(self) -> str:
		if self.unit is None:
			return self.display_name
		else:
			return f"{self.display_name} ({self.unit})"


def _load_results(dataset: str,
) -> tuple[KeyCfgLib, list[str], dict[KeyCfgLib.KeyCfg, pandas.DataFrame]]:
	fname = f"output/figure_6.{dataset}_filter_imput_transfer.age_pred.pkl"
	with open(fname, "rb") as fp:
		raw = pickle.load(fp)

	key_fields = raw["key_fields"]

	key_cfgs = KeyCfgLib()
	res = dict()

	for k, df in raw["res"].items():
		cfg = key_cfgs.add_new(**dict(zip(key_fields, k)))
		res[cfg] = df

	return key_cfgs, raw["clocks"], res


def setup_layout(metric_cfgs: list[MetricCfg],
	n_key_cfgs: int, n_key_fields: int,
) -> dict:
	lc = mpllayout.LayoutCreator(
		left_margin=0.8,
		right_margin=0.2,
		bottom_margin=2.0,
		top_margin=0.4,
	)

	nrows = len(metric_cfgs)

	key_field_cfg_sp = 0.3
	key_field_elem_sp = 0.3

	# horizontal layout
	key_field_w = key_field_elem_sp * n_key_fields
	key_field_h = (key_field_cfg_sp * n_key_cfgs)
	key_field_gap = 0.1

	box_w = 3.5
	box_h = key_field_h

	col_gap_w = 0.8
	col_w = key_field_w + key_field_gap + box_w

	for i, cfg in enumerate(metric_cfgs):
		key_field = lc.add_frame(f"key_field_{cfg.key}")
		key_field.set_size(key_field_w, key_field_h)
		key_field.set_anchor("bottomleft", offsets=((col_w + col_gap_w) * i, 0))

		box = lc.add_frame(f"box_{cfg.key}")
		box.set_size(box_w, box_h)
		box.set_anchor("bottomleft", key_field, "bottomright",
		    offsets=(key_field_gap, 0))

	layout = lc.create_figure_layout()
	layout["nrows"] = nrows

	return layout


def _calc_stats(data: dict[tuple, pandas.DataFrame], clocks: Sequence[str],
	*, repl_group: pylib.dataset.ReplGroup = None,
) -> dict[tuple, pandas.DataFrame]:
	ret = dict()
	for k, df in data.items():
		age_pred_res = pylib.age_pred_res.AgePredRes(df)
		ret[k] = age_pred_res.calc_stat(clocks, repl_group=repl_group)
	return ret


def _sort_stat_keys(stats: dict[tuple, pandas.DataFrame], metric: MetricCfg, *,
	reverse: bool = False,
) -> list[KeyCfgLib.KeyCfg]:

	ret = sorted(stats.keys(), key=lambda k: k.order_key,
		# key=lambda k: numpy.median(stats[k][metric.key]),
		reverse=reverse)

	return ret


def _plot_stats(box_axes: matplotlib.axes.Axes,
	key_field_axes: matplotlib.axes.Axes, *,
	stats: dict[tuple, pandas.DataFrame],
	metric: MetricCfg,
	key_order=list[KeyCfgLib.KeyCfg],
	key_field_cfgs: KeyFieldList,
	highlight_clock_cfgs: pylib.ClockCfgLib,
	box_ec: str = "#000000",
	box_fc: str = "#c0c0c0",
) -> NoReturn:

	n_keys = len(key_order)
	n_key_fields = len(key_field_cfgs)

	# plot key field patches
	for i, key in enumerate(key_order):
		for j, field in enumerate(key_field_cfgs):
			if field.is_active_in(key):
				p = matplotlib.patches.Circle((j + 0.5, i + 0.5), 0.3,
					linewidth=0.5, edgecolor="#000000", facecolor=field.color,
					zorder=10
				)
				key_field_axes.add_patch(p)

	# add horizontal arrows
	for i in range(n_keys):
		a = matplotlib.patches.FancyArrowPatch(
			(-0.2, i + 0.5), (n_key_fields + 0.4, i + 0.5),
			arrowstyle='-|>', mutation_scale=10, linewidth=2.0, color="#c0c0c0",
			clip_on=False, zorder=5,
		)
		key_field_axes.add_patch(a)

	# misc
	key_field_axes.tick_params(
		left=False, labelleft=True,
		right=False, labelright=False,
		bottom=True, labelbottom=True,
		top=False, labeltop=False,
	)
	for sp in key_field_axes.spines.values():
		sp.set_visible(False)
	key_field_axes.set(
		xlim=(0, len(key_field_cfgs)),
		ylim=(0, n_keys),
	)
	key_field_axes.set_xticks(numpy.arange(n_key_fields) + 0.5, fontsize=8,
		rotation=90,
		labels=[f.full_name for f in key_field_cfgs],
	)
	key_field_axes.set_yticks(numpy.arange(n_keys) + 0.5, fontsize=8, labels=[
		(("+").join(key_field_cfgs.get_field_labels(k)) or "Direct") for k in key_order
	])

	# plot boxes
	for i, key in enumerate(key_order):
		vals = stats[key][metric.key].values

		box = box_axes.boxplot(vals,
			orientation="horizontal", positions=[i + 0.5],
			widths=0.8, showfliers=False, patch_artist=True, whis=(0, 100),
			boxprops=dict(linewidth=0.5, edgecolor=box_ec, facecolor=box_fc),
			medianprops=dict(linewidth=1.0, color=box_ec),
			whiskerprops=dict(linewidth=0.5, color=box_ec),
			capprops=dict(linewidth=0.5, color=box_ec),
			capwidths=0.7,
		)

	# add clock scatters
	clocks = sorted(stats[key_order[0]].index)
	y = numpy.arange(n_keys) + 0.5

	handles = list()
	for c in clocks:
		x = [stats[key].loc[c, metric.key] for key in key_order]

		if c in highlight_clock_cfgs:
			cfg = highlight_clock_cfgs[c]
			marker = cfg.marker
			edgecolor = "#000000"
			facecolor = cfg.color + "a0"
			zorder = 12
			label = cfg.name
		else:
			marker = "."
			edgecolor = "none"
			facecolor = box_ec
			zorder = 10
			label = None

		h = box_axes.scatter(x, y, marker=marker, s=24, linewidth=0.5,
			edgecolors=edgecolor, facecolors=facecolor,
			zorder=zorder, label=label,
		)

		if c in highlight_clock_cfgs:
			handles.append(h)

	# clock legend
	box_axes.figure.legend(handles=handles, loc=4, bbox_to_anchor=(1.00, 0.00),
		frameon=True, ncol=3, fontsize=8, handlelength=0.75, title_fontsize=10,
		title="Clock markers",
	)

	# misc
	box_axes.tick_params(
		left=False, labelleft=False,
		right=False, labelright=False,
		bottom=True, labelbottom=True,
		top=False, labeltop=False,
	)
	box_axes.set_xlim(0, None)
	box_axes.set_ylim(0, n_keys)
	box_axes.set_yticks(numpy.arange(n_keys) + 0.5)
	box_axes.set_xlabel(metric.label)

	return


def _plot_main(dataset: str, plot_bname: str, *,
	highlight_clocks: set[str] = None, with_rd: bool = False,
):
	# config
	metric_cfgs = [MetricCfg("mae", "MAE", unit="years")]
	if with_rd:
		metric_cfgs.append(MetricCfg("repl_mad", "RD", unit="years"))

	key_field_cfgs = KeyFieldList([
		KeyFieldList.KeyFieldCfg(key="depth_filtering", short_name="DF",
			full_name="Depth Filtering", color="#c4fbff",
			# get_active_label=lambda self, v: f"{self.short_name}$_{{{getattr(v, self.key)}}}$"
		),
		KeyFieldList.KeyFieldCfg(key="imputation", short_name="IM",
			full_name="0/1 Imputation", color="#daffc3",
		),
		KeyFieldList.KeyFieldCfg(key="transfer_learning", short_name="TL",
			full_name="Transfer Learning", color="#ecb3ff",
		),
	])

	# load results
	key_cfgs, clocks, results = _load_results(dataset)
	clock_cfgs = pylib.ClockCfgLib.from_json(select=highlight_clocks)
	repl_group = pylib.dataset.ReplGroup.load_dataset(dataset) if with_rd else None
	dataset_cfgs = pylib.DatasetCfgLib.from_json()

	# prep data
	stats = _calc_stats(results, clocks, repl_group=repl_group)

	# plot
	layout = setup_layout(metric_cfgs, len(key_cfgs), len(key_field_cfgs))
	figure = layout["figure"]

	for i, metric in enumerate(metric_cfgs):
		key_order = _sort_stat_keys(stats, metric, reverse=False)

		_plot_stats(
			layout[f"box_{metric.key}"],
			layout[f"key_field_{metric.key}"],
			stats=stats,
			metric=metric,
			key_order=key_order,
			key_field_cfgs=key_field_cfgs,
			highlight_clock_cfgs=clock_cfgs,
			box_ec="#a0a0a0",
			box_fc="#e0e0e0",
		)

	# figure title
	figure.suptitle(
		f"Training/testing: {dataset_cfgs[dataset].display_name} (N={dataset_cfgs[dataset].n_samples}), 5-fold cross-validation",
		fontsize=12,
	)

	figure.savefig(plot_bname + ".png", dpi=600)
	figure.savefig(plot_bname + ".svg", dpi=600)
	figure.savefig(plot_bname + ".pdf", dpi=600)
	matplotlib.pyplot.close(figure)
	return


def _main():
	args = get_args()
	_plot_main(
		dataset=args.dataset,
		plot_bname=args.plot_basename,
		with_rd=args.with_rd,
	)

	return


if __name__ == "__main__":
	_main()
