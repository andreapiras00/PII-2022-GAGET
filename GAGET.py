import json
import sys
import os
import threading
import networkx as nx
# import gfapy
import re
import cairo
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.cm as cm
import graph_tool.all as gt
import numpy as np
from copy import copy, deepcopy
from optparse import OptionParser, OptionGroup, Option
from collections import defaultdict
from operator import itemgetter
from multiprocessing import Pool
from functools import partial
from scripts.config import *
from scripts.graphaligner_runner import *
from scripts.utils import parse_assembler_output
from scripts.alignments_parser import *
from scripts.graph_analysis import best_alignments_selection
from gi.repository import Gtk, GdkPixbuf, Pango
from timeit import default_timer as timer
from matplotlib.backends.backend_gtk3agg import (FigureCanvasGTK3Agg as FigureCanvas)
from matplotlib.figure import Figure


timer_list = [] 
start_total = timer()

class GAGETOption(Option):
	''' Class used to parse options from terminal '''
	def check_file(option, opt, fpath):
		if fpath and os.path.isdir(fpath):
			print("ERROR! You specify a folder instead of a file: %s" % (fpath))
			sys.exit(2)
		if not fpath or not os.path.isfile(fpath):
			print("ERROR! File not found: %s" % (fpath))
			sys.exit(2)
		return fpath

	def check_dir(option, opt, dirpath):
		if not dirpath or not os.path.isdir(dirpath):
			print("ERROR! Folder not found: %s" % (dirpath))
			sys.exit(2)
		return dirpath

	TYPES = Option.TYPES + ('file', 'dir')
	TYPE_CHECKER = copy(Option.TYPE_CHECKER)
	TYPE_CHECKER['file'] = check_file
	TYPE_CHECKER['dir'] = check_dir


class GAGET_Window(Gtk.Window):
	''' Class GAGET_Window represents the GUI of the software as well as the 
		logic behind used to interact with the graph. '''

	def __init__(self, g, geometry, gr_stats=[], ge_stats=[], paths={}, al_info={}, v_ID={}, vprops=None, eprops=None, vorder=None,
				 eorder=None, nodesfirst=False, update_layout=False, **kwargs):
		''' Method __init__ initializes all the different sections of the GUI using the info passed as arguments. '''

		self.g = g
		self.paths = paths
		self.al_info = al_info
		self.geometry = geometry
		self.vprops = vprops
		self.eprops = eprops
		self.vorder = vorder
		self.eorder = eorder
		self.nodesfirst = nodesfirst
		self.update_layout = update_layout
		self.kwargs = kwargs
		self.ref_length = {}
		self.scaling_factor = {}
		self.refs = []
		self.is_fullscreen = False
		# values used to represent the genome coverage of each node as an array of 0s and 1s.
		# for each ref, its length and the consequent scaling factor are stored
		for ref in ge_stats.keys():
			self.ref_length[ref] = min(int(ge_stats[ref][0][1]), 10000)
			self.scaling_factor[ref] = self.ref_length[ref]/int(ge_stats[ref][0][1])
			self.refs.append(ref)
		
		# initialize the window with the app logo and the default sizes
		Gtk.Window.__init__(self, title="Quality Assessment for Genome Assembly Graphs")
		self.set_border_width(10)
		icon = GdkPixbuf.Pixbuf.new_from_file('logo.png')
		self.set_icon(icon)
		self.set_default_size(geometry[0], geometry[1])
		# self.fullscreen()

		# initialize root grid
		self.grid = Gtk.Grid()
		self.grid.set_column_homogeneous(True)
		self.grid.set_row_homogeneous(True)
		self.grid.set_column_spacing(10)
		self.grid.set_row_spacing(10)
		self.add(self.grid)
		
		# initialize graph stats section as a ListStore
		self.graph_stats = Gtk.ListStore(str, str)
		for stat in gr_stats:
			self.graph_stats.append(list(stat))
		# initialize TreeView
		self.treeview_gr = Gtk.TreeView.new_with_model(self.graph_stats)
		renderer = Gtk.CellRendererText(weight=430)
		for i, column_title in enumerate(["Metric", "Value"]):
			column = Gtk.TreeViewColumn(column_title, cell_renderer=renderer, text=i)
			column.set_alignment(0.5)
			self.treeview_gr.append_column(column)
			renderer = Gtk.CellRendererText(weight=370)
		# attach TreeView to scrollable window
		self.scrollable_treelist_gr = Gtk.ScrolledWindow()
		self.scrollable_treelist_gr.set_vexpand(True)
		self.scrollable_treelist_gr.add(self.treeview_gr)
		# create graph stats box and attach title and scrollable window
		self.gr_st_box = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 10)
		self.gr_st_box.set_homogeneous(False)
		label = Gtk.Label()
		label.set_markup("<b>Graph Statistics</b>")
		label.set_justify(Gtk.Justification.LEFT)
		self.gr_st_box.pack_start(label, False, False, 0)
		self.gr_st_box.pack_start(self.scrollable_treelist_gr, True, True, 0)

		# initialize genome stats section as a dictionary of ListStore
		self.genome_stats = {}
		for ref, stats in ge_stats.items():
			self.genome_stats[ref] = Gtk.ListStore(str, str)
			for stat in stats:
				self.genome_stats[ref].append(list(stat))
		# initialize TreeView with first (default) reference 
		self.treeview_ge = Gtk.TreeView.new_with_model(self.genome_stats[self.refs[0]])
		renderer = Gtk.CellRendererText(weight=430)
		for i, column_title in enumerate(["Metric", "Value"]):
			column = Gtk.TreeViewColumn(column_title, cell_renderer=renderer, text=i)
			column.set_alignment(0.5)
			self.treeview_ge.append_column(column)
			renderer = Gtk.CellRendererText(weight=370)
		# attach TreeView to scrollable window
		self.scrollable_treelist_ge = Gtk.ScrolledWindow()
		self.scrollable_treelist_ge.set_vexpand(True)
		self.scrollable_treelist_ge.set_hexpand(True)
		self.scrollable_treelist_ge.add(self.treeview_ge)
		# initialize button to switch between references as ComboBoxText
		self.refcombo = Gtk.ComboBoxText()
		self.refcombo.set_entry_text_column(0)
		self.refcombo.connect("changed", self.on_reference_changed)
		for ref in ge_stats.keys():
			self.refcombo.append_text(ref)
		# default reference is the first one
		self.refcombo.set_active(0)
		# create genome stats box and attach title, scrollable window and 
		# references button only if there are more than 1 reference
		self.ge_st_box = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 10)
		self.ge_st_box.set_homogeneous(False)
		label = Gtk.Label()
		label.set_markup("<b>Genome Statistics</b>")
		label.set_justify(Gtk.Justification.LEFT)
		if(len(self.refs) >= 2):
			box = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 10)
			box.pack_start(label, False, False, 0)
			box.pack_end(self.refcombo, True, True, 0)
			self.ge_st_box.pack_start(box, False, False, 0)
		else:
			self.ge_st_box.pack_start(label, False, False, 0)
		self.ge_st_box.pack_start(self.scrollable_treelist_ge, True, True, 0)

		# create filter box
		self.fil_box = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 10)
		label = Gtk.Label()
		label.set_markup("<b>Filters:</b>")
		self.prop_list = ["length", "coverage"]
		# initialize filters list as ListStore
		self.filters = Gtk.ListStore(bool, str, int)
		self.filters.append([False, "Length ≥", 0])
		self.filters.append([False, "% aligned bases ≥", 0])
		self.filters.append([False, "Aligned path", 1])
		# initialize TreeView
		self.filt_tree = Gtk.TreeView.new_with_model(self.filters)
		self.renderer_check = Gtk.CellRendererToggle()
		self.renderer_check.connect("toggled", self.on_cell_toggled)
		renderer_text = Gtk.CellRendererText()
		renderer_editable = Gtk.CellRendererText()
		renderer_editable.set_property("editable", True)
		renderer_editable.connect("edited", self.text_edited)
		column = Gtk.TreeViewColumn("", cell_renderer=self.renderer_check, active=0)
		self.filt_tree.append_column(column)
		column = Gtk.TreeViewColumn("Property", cell_renderer=renderer_text, text=1)
		column.set_alignment(0.5)
		self.filt_tree.append_column(column)
		column = Gtk.TreeViewColumn("Threshold", cell_renderer=renderer_editable, text=2)
		column.set_alignment(0.5)
		self.filt_tree.append_column(column)
		# attach TreeView to scrollable window
		self.scrollable_treelist_fil = Gtk.ScrolledWindow()
		self.scrollable_treelist_fil.set_vexpand(True)
		self.scrollable_treelist_fil.set_hexpand(True)
		self.scrollable_treelist_fil.add(self.filt_tree)
		# attach title and scrollable window to filters box
		self.fil_box.pack_start(label, False, False, 0)
		self.fil_box.pack_start(self.scrollable_treelist_fil, True, True, 0)
		
		# create paths list box
		self.path_box = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 10)
		label = Gtk.Label()
		label.set_markup("<b>Aligned paths:</b>")
		# initialize button to deselect paths
		self.reset_paths = Gtk.Button.new_with_label("Deselect Path")
		self.reset_paths.connect("clicked", self.deselect_path)
		# initialize paths list as ListStore
		self.aligns = Gtk.ListStore(int,int)
		for path in paths:
			self.aligns.append(list((path,len(paths[path]))))
		# initialize TreeView with paths list
		self.treeview_paths = Gtk.TreeView.new_with_model(self.aligns)
		renderer = Gtk.CellRendererText(weight=430)
		for i, column_title in enumerate(["Path", "# nodes"]):
			column = Gtk.TreeViewColumn(column_title, cell_renderer=renderer, text=i)
			column.set_alignment(0.5)
			column.set_sort_column_id(i)
			self.treeview_paths.append_column(column)
			renderer = Gtk.CellRendererText(weight=370)
		# connect path selection to method	
		self.select = self.treeview_paths.get_selection()
		self.select.connect("changed", self.path_selection_changed)
		# attach TreeView to scrollable window
		self.scrollable_treelist_paths = Gtk.ScrolledWindow()
		self.scrollable_treelist_paths.set_vexpand(True)
		self.scrollable_treelist_paths.set_hexpand(True)
		self.scrollable_treelist_paths.add(self.treeview_paths)
		# attach title, scrollable window and reset button to paths list box
		self.path_box.pack_start(label, False, False, 0)
		self.path_box.pack_start(self.reset_paths, False, False, 0)
		self.path_box.pack_start(self.scrollable_treelist_paths, True, True, 0)

		# initialize graph's info for graph widget
		self.v_size = gt.prop_to_size(g.vertex_properties["length"], mi=10, ma=50, power=1)
		self.v_col = g.new_vertex_property("vector<float>")
		# initialize colormap for coverage
		cov_min = 0.000000001
		cov_max = max([g.vertex_properties["coverage"][v] for v in g.vertices()])
		colormap = cm.viridis
		norm = mcol.Normalize(vmin=cov_min, vmax=cov_max)
		m = cm.ScalarMappable(norm = norm, cmap = colormap)
		# assign a color to each vertix based on its coverage value
		for v in g.vertices():
			self.v_col[v] = list(m.to_rgba(g.vertex_properties["coverage"][v]+0.000000001))
		# store original color of nodes	
		self.original_col = deepcopy(self.v_col)
		# initialize graph widget
		self.pos_sfdp = gt.sfdp_layout(self.g)
		self.graph = gt.GraphWidget(self.g, self.pos_sfdp, vprops, eprops, vorder, eorder,
								 nodesfirst, update_layout, vertex_size = self.v_size, vertex_fill_color = self.v_col, **kwargs)
		# store original porision of each node
		self.original_pos = deepcopy(self.graph.pos)
		# attach node selection to method
		self.graph.connect("button-press-event", self.print_info)
		# create graph box and attach all required buttons and graph widget
		self.gr_box = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 10)
		self.gr_box.set_homogeneous(False)
		label = Gtk.Label()
		label.set_markup("<b>Assembly Graph</b>")
		# initialize button to fit graph inside box
		fit_view = Gtk.Button.new_with_label("Fit View")
		fit_view.connect("clicked", self.fit_v)
		# initialize button to reset zoom
		reset_view = Gtk.Button.new_with_label("Reset View")
		reset_view.connect("clicked", self.reset_v)
		# initialize button to reset the layout
		reset_graph = Gtk.Button.new_with_label("Reset Layout")
		reset_graph.connect("clicked", self.reset_layout)
		# initialize button to deselect a node
		self.reset_node = Gtk.Button.new_with_label("Deselect Node")
		self.reset_node.connect("clicked", self.deselect_node)
		# initialize figure for coverage's colormap
		gradient = np.linspace(0, 1, 256)
		gradient = np.vstack((gradient, gradient))
		fig = Figure()
		ax = fig.add_subplot()
		fig.subplots_adjust(top=1, bottom=0, left=0, right=1)
		ax.imshow(gradient, aspect='auto', cmap=colormap)
		ax.set_axis_off()
		canvas = FigureCanvas(fig)
		canvas.set_size_request(300, 30)
		label1 = Gtk.Label()
		label1.set_markup("<i>Coverage</i>")
		label2 = Gtk.Label()
		label2.set_text("0.0")
		label3 = Gtk.Label()
		label3.set_text("{:.1f}".format(cov_max))
		hbox = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 0)
		hbox.pack_end(label3, False, False, 6)
		hbox.pack_end(canvas, False, False, 0)
		hbox.pack_end(label2, False, False, 6)
		cov_legend = Gtk.Frame()
		cov_legend.add(hbox)
		# attach all buttons and coverage to utils box
		utils_box = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
		utils_box.pack_start(fit_view, False, False, 0)
		utils_box.pack_start(reset_graph, False, False, 0)
		utils_box.pack_start(reset_view, False, False, 0)
		utils_box.pack_start(self.reset_node, False, False, 0)
		utils_box.pack_end(cov_legend, False, False, 0)
		utils_box.pack_end(label1, False, False, 0)
		# attach title, utils box and graph to graph box
		self.gr_box.pack_start(label, False, False, 0)
		self.gr_box.pack_start(utils_box, False, False, 0)
		self.gr_box.pack_end(self.graph, True, True, 0)

		# create reference box and attach name and ref figure
		self.ref_box = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 1)
		self.ref_name = Gtk.Label()
		self.ref_name.set_markup('<b>Reference Genome:</b>' + '\t' + self.refs[0])
		self.ref_box.pack_start(self.ref_name, False, False, 0)
		# plot default reference to its default value
		self.initial_ref_plot()

		# initialize TextView for nodes and paths info
		self.textview = Gtk.TextView()
		self.textview.set_editable(False)
		self.textview.set_cursor_visible(False)
		self.textview.set_wrap_mode(Gtk.WrapMode.CHAR)
		self.textview.set_left_margin(10)
		self.textview.set_top_margin(10)
		self.textview.set_right_margin(10)
		self.textbuffer = self.textview.get_buffer()
		# initialize info window and attach TextView to it
		info_window = Gtk.ScrolledWindow()
		info_window.set_hexpand(True)
		info_window.set_vexpand(True)
		info_window.add(self.textview)
		
		# attach all boxes and window to root grid
		self.grid.attach(self.gr_st_box, 0, 0, 6, 14)
		self.grid.attach(self.ge_st_box, 0, 14, 6, 10)
		self.grid.attach(self.gr_box, 6, 0, 24, 14)
		self.grid.attach(self.ref_box, 6, 14, 24, 1)
		self.grid.attach(info_window, 6, 15, 24, 9)
		self.grid.attach(self.fil_box, 30, 0, 6, 8)
		self.grid.attach(self.path_box, 30, 8, 6, 16)

		self.tag_bold = self.textbuffer.create_tag("bold", weight=Pango.Weight.BOLD, scale=1.1)
		self.tag_normal = self.textbuffer.create_tag("nenti", scale=1.1)

		self.connect('key-press-event', self.fullscreen_toggler)


	def on_reference_changed(self, combo):
		''' on_reference_changed method updates the genome info based on the reference selected from the combobox '''
		ref = combo.get_active_text()
		if ref is not None:
			self.scrollable_treelist_ge.remove(self.treeview_ge)
			self.treeview_ge.set_model(self.genome_stats[ref])
			self.scrollable_treelist_ge.add(self.treeview_ge)

	# function used to initialize the plot of the reference
	def initial_ref_plot(self):
		''' initial_ref_plot method initializes the FigureCanvas for the default reference using an 
			array of 0s '''
		self.fig = Figure(figsize=(20, 5))
		self.ax = self.fig.add_subplot()
		self.fig.subplots_adjust(top=1, bottom=0, left=0, right=1)
		# delete all unnecessary lines of the plot and reduce empty spaces
		self.ax.spines['top'].set_visible(False)
		self.ax.spines['right'].set_visible(False)
		self.ax.spines['left'].set_visible(False)
		self.ax.spines['bottom'].set_visible(False)
		self.ax.margins(x=0)
		self.ax.tick_params(labelleft=False, left=False, bottom=False, labelbottom=False)
		self.fig.tight_layout()
		self.canvas = FigureCanvas(self.fig)
		self.fig.canvas.set_size_request(1250, 25)
		self.ref_box.pack_start(self.canvas, True, True, 0)
		self.update_ref(np.zeros(self.ref_length[self.refs[0]], dtype=int))

	def update_ref(self, ref):
		''' update_ref method updates the reference figure based on the array passed as argument '''
		gradient = np.vstack((ref, ref))
		self.ax.clear()
		colormap = cm.Reds
		self.ax.imshow(gradient, aspect='auto', cmap=colormap)
		self.fig.canvas.draw()

	def print_info(self, widget, event):
		''' print_info method looks for the selected node, parses and adds its info to the info TextView, updates the reference
			based on the node's mapped positions and updates the graph's node colors '''
		# check that there's no path selected	
		selection = self.treeview_paths.get_selection()
		model, treeiter = selection.get_selected()
		# clear TextView
		self.textbuffer.set_text("")
		if treeiter is not None:
			return
		self.textbuffer.place_cursor(self.textbuffer.get_start_iter())
		# current_filter = self.g.get_vertex_filter()[0]	
		orange = [0.807843137254902, 0.3607843137254902, 0.0, 1.0]
		gray = [0.6, 0.6, 0.6, 1.0]
		vcolor = self.v_col
		# find the selected node and color it orange, all other nodes are colored gray
		for v in self.g.vertices():
			if self.graph.selected[v]:
				vcolor[v] = orange
				sel_node = v
			else:
				vcolor[v] = gray	
		# store selected node's info 			
		text = [("Node ID", " = {}\n\n".format(self.g.vp['ID'][sel_node]))]
		text.append(("Length", " = {}\n".format(self.g.vp['length'][sel_node])))
		text.append(("Unique length", " = {}\n\n".format(self.g.vp['u_length'][sel_node])))
		text.append(("Coverage", " = {:.4f}\n\n".format(self.g.vp['coverage'][sel_node])))
		text.append(("GC content", " = {:.2f}%\n\n".format(self.g.vp['GC'][sel_node]*100)))
		text.append(("MAPPING INFO", "\n"))
		# find selected node's reference and positions mapped
		posmppd = self.g.vp['pos_mapped'][sel_node]
		if posmppd:
			ref = posmppd[next(iter(posmppd))][0][4]
		else:
			ref = self.refs[0]
		# update selected reference	
		self.refcombo.set_active(self.refs.index(ref))	
		self.ref_name.set_markup('<b>Reference Genome:</b>' + '\t' + ref)	
		new_ref = np.zeros(self.ref_length[ref], dtype=int)
		# for each path in positions mapped, update reference and store path's info
		for path in posmppd:
			for el in posmppd[path]:
				#TODO: should this be done only for the node's coverage or for all paths?
				start = int(el[5]*self.scaling_factor[ref])
				end = int(el[6]*self.scaling_factor[ref]) 
				# update necessary caused by the scaling factor's rounding
				if start == end:
					if end != self.ref_length[ref]:
						end += 1
					else:
						start -= 1
				new_ref[start:end] = 1
				text.append(("\nPath", " = {}\n".format(path)))
				text.append(("Alignment length", " = {}\n".format(el[6]-el[5]+1)))
				text.append(("Start position (node)", " = {}\n".format(el[0])))
				text.append(("End position (node)", " = {}\n".format(el[1])))
				text.append(("Start position (ref)", " = {}\n".format(el[5])))
				text.append(("End position (ref)", " = {}\n".format(el[6])))
				text.append(("Reverse complement", " = {}\n".format(el[7])))
				text.append(("Score", " = {}\n".format(el[3])))
				text.append(("CIGAR", " = '{}'\n".format(el[2])))
		# disable paths list window and path deselection button		
		self.reset_paths.set_sensitive(False)
		self.scrollable_treelist_paths.set_sensitive(False)
		# redraw reference figure
		self.update_ref(new_ref)		
		text.append(("\nSequence\n", "{}\n".format(self.g.vp['sequence'][v])))
		# add stored info to TextView
		for line in text:
			insert = self.textbuffer.get_iter_at_mark(self.textbuffer.get_insert())
			self.textbuffer.insert_with_tags(insert, line[0], self.tag_bold)
			insert = self.textbuffer.get_iter_at_mark(self.textbuffer.get_insert())
			self.textbuffer.insert_with_tags(insert, line[1], self.tag_normal)
		# redraw graph	
		self.graph.regenerate_surface()
		self.graph.queue_draw()	
	
	# function used to parse the info from the selected path and write it on the info box
	def path_selection_changed(self, selection):
		''' path_selection_changed method parses info from the selected path and adds its nodes info to the info 
			TextView, updates the reference based on the nodes' mapped positions and updates the graph's nodes colors '''
		# update filter and clear TextView
		current_filter = self.g.get_vertex_filter()[0]
		self.g.set_vertex_filter(None)
		self.textbuffer.set_text("")
		self.textbuffer.place_cursor(self.textbuffer.get_start_iter())
		model, treeiter = selection.get_selected()
		vcolor = self.v_col
		# check that a path has been selected
		if treeiter is not None:
			path = model[treeiter][0]
			orange = [0.807843137254902, 0.3607843137254902, 0.0, 1.0]
			gray = [0.6, 0.6, 0.6, 1.0]
			# find nodes in the path and color them orange, all other nodes are colored gray
			for v in self.g.vertices():
				if self.g.vp['ID'][v] in self.paths[path]:
					vcolor[v] = orange
				else:
					vcolor[v] = gray
			# find path's reference
			alignment = self.al_info[path]
			text = [("Path", " = {}\n".format(path))]
			ref = alignment[0][5]
			# update selected reference	
			self.refcombo.set_active(self.refs.index(ref))
			self.ref_name.set_markup('<b>Reference Genome:</b>' + '\t' + ref)
			new_ref = np.zeros(self.ref_length[ref], dtype=int)
			# for each node in path, update reference and store its info
			for node in alignment:
				start = int(node[6]*self.scaling_factor[ref])
				end = int(node[7]*self.scaling_factor[ref]) 
				if start == end:
					if end != self.ref_length[ref]:
						end += 1
					else:
						start -= 1
				new_ref[start:end] = 1
				text.append(("\nNode ID", " = {}\n".format(node[0])))
				text.append(("Alignment length", " = {}\n".format(node[7]-node[6]+1)))
				text.append(("Start position (node)", " = {}\n".format(node[1])))
				text.append(("End position (node)", " = {}\n".format(node[2])))
				text.append(("Start position (ref)", " = {}\n".format(node[6])))
				text.append(("End position (ref)", " = {}\n".format(node[7])))
				text.append(("Reverse complement", " = {}\n".format(node[8])))
				text.append(("Score", " = {}\n".format(node[4])))
				text.append(("CIGAR", " = '{}'\n".format(node[3])))
			# add stored info to TextView	
			for line in text:
				insert = self.textbuffer.get_iter_at_mark(self.textbuffer.get_insert())
				self.textbuffer.insert_with_tags(insert, line[0], self.tag_bold)
				insert = self.textbuffer.get_iter_at_mark(self.textbuffer.get_insert())
				self.textbuffer.insert_with_tags(insert, line[1], self.tag_normal)
		# if no path has been selected set reference and nodes color to their default values		
		else:
			new_ref = np.zeros(self.ref_length[self.refs[0]], dtype=int)
			for v in self.g.vertices():
				vcolor[v] = self.original_col[v]	
		# disable node deselection button		
		self.reset_node.set_sensitive(False)	
		# redraw reference figure	
		self.update_ref(new_ref)
		# update filter and redraw graph
		self.g.set_vertex_filter(current_filter)
		self.graph.regenerate_surface()
		self.graph.queue_draw()

	
	def text_edited(self, widget, path, text):
		''' text_edited method updates the filter based on the number inputted and activates it if it's
			a valid number '''
		toggle = True
		print(path, self.filters[path][0], self.prop_list)
		# checks if the paths filter is being modified
		if (int(path) == len(self.prop_list)):
			# checks if the path number is present
			if (int(text) in self.paths):
				self.filters[path][2] = int(text)
			else:
				toggle = False
		else:
			# update threshold value
			self.filters[path][2] = int(text)
		if toggle:
			# apply filter
			self.on_cell_toggled(widget, path)
			if not self.filters[path][0]:
				self.on_cell_toggled(widget, path)


	def on_cell_toggled(self, widget, path):
		''' on_cell_toggled method applies the filter to the graph's nodes based on the threshold value '''
		self.filters[path][0] = not self.filters[path][0]
		if self.filters[path][0]:
			if int(path) < len(self.prop_list):
				prop = self.prop_list[int(path)]
				thres = self.filters[path][2]
				if prop == 'coverage':
					thres = float(thres)/100
				filter_new = self.g.get_vertex_filter()[0]
				if filter_new:
					for v in self.g.vertices():
						if self.g.vp[prop][v] < thres:
							filter_new[v] = filter_new[v] and False
						else:
							filter_new[v] = filter_new[v] and True
				else:
					filter_new = self.g.new_vertex_property("bool")
					for v in self.g.vertices():
						if self.g.vp[prop][v] < thres:
							filter_new[v] = False
						else:
							filter_new[v] = True
			else:
				al = self.filters[path][2]
				filter_new = self.g.get_vertex_filter()[0]
				if filter_new:
					for v in self.g.vertices():
						if self.g.vp['ID'][v] in self.paths[al]:
							filter_new[v] = filter_new[v] and True
						else:
							filter_new[v] = filter_new[v] and False
				else:
					filter_new = self.g.new_vertex_property("bool")						
					for v in self.g.vertices():
						if self.g.vp['ID'][v] in self.paths[al]:
							filter_new[v] = True
						else:
							filter_new[v] = False
		else:
			self.g.set_vertex_filter(None)
			filter_new = self.g.new_vertex_property("bool")
			for v in self.g.vertices():
				filter_new[v] = True
			for row in self.filters:
				if self.filters[row.path][0]:
					if int(str(row.path)) < len(self.prop_list):
						prop = self.prop_list[int(str(row.path))]
						thres = self.filters[row.path][2]
						if prop == 'coverage':
							thres = float(thres)/100
						for v in self.g.vertices():
							if self.g.vp[prop][v] < thres:
								filter_new[v] = filter_new[v] and False
							else:
								filter_new[v] = filter_new[v] and True
					else:
						al = self.filters[row.path][2]
						for v in self.g.vertices():
							if self.g.vp['ID'][v] in self.paths[al]:
								filter_new[v] = filter_new[v] and True
							else:
								filter_new[v] = filter_new[v] and False
		self.g.set_vertex_filter(filter_new)
		self.graph.regenerate_surface()
		self.graph.queue_draw()
	
	def reset_layout(self, button):
		''' reset_layout method resets the graph's layout to the default one '''
		self.graph.draw(0,cr=cairo.Context(self.graph.base))
		for v in self.g.vertices():
			self.graph.vertex_matrix.update_vertex(self.g.vertex(int(v)), self.original_pos[v].a)
		self.graph.regenerate_surface()
		self.graph.queue_draw()
	
	def fit_v(self, button):
		''' fit_v method updates the graph's view in order to fit it to the box space available '''
		self.graph.draw(0,cr=cairo.Context(self.graph.base))
		self.graph.fit_to_window()
		self.graph.regenerate_surface()
		self.graph.queue_draw()
	
	def reset_v(self, button):
		''' reset_v method resets the graph's view to the default one '''
		self.graph.draw(0,cr=cairo.Context(self.graph.base))
		self.g.set_vertex_filter(None)
		self.select.unselect_all()
		for row in self.filters:
			# remove filters
			self.filters[row.path][0] = False
		self.graph.fit_to_window()
		self.graph.regenerate_surface()
		self.graph.queue_draw()
	
	def deselect_path(self, button):
		''' deselect_path method deselects the current path (automatically calls path_selection_changed method)'''
		self.select.unselect_all()
		self.reset_node.set_sensitive(True)

	def deselect_node(self, button):
		''' deselect_node method deselects the current node and clears the info section '''
		vcolor = self.v_col
		self.textbuffer.set_text("")
		self.textbuffer.place_cursor(self.textbuffer.get_start_iter())
		for v in self.g.vertices():
			vcolor[v] = self.original_col[v]
		self.reset_paths.set_sensitive(True)
		self.scrollable_treelist_paths.set_sensitive(True)		
		self.graph.regenerate_surface()
		self.graph.queue_draw()		
		new_ref = np.zeros(self.ref_length[self.refs[0]], dtype=int)	
		self.update_ref(new_ref)

	def fullscreen_toggler(self, widget, event):
		''' fullscreen_toggler method sets to fullscreen (or viceversa) the window when the 'f10' key is pressed '''
		if event.keyval == 65479:
			if self.is_fullscreen:
				self.unfullscreen()
			else:
				self.fullscreen()
			self.is_fullscreen = not self.is_fullscreen	

	def __del__(self):
		self.graph.cleanup()


def main():
	description = ('Quality Assessment for Genome Assembly Graphs')
	print("\n{}\n".format(description))

	parser = OptionParser(description=description, option_class=GAGETOption)

	group = OptionGroup(parser, "Options")#, "Options that can be used in any mode")

	group.add_option('-o', dest='output_dir', help='Output directory [default: gaget_output]', default='gaget_output')
	group.add_option('-m', dest='min_edge_len', help='Lower threshold for edge length [default: %d]' % MIN_EDGE_LEN, default=MIN_EDGE_LEN)
	group.add_option('-t', dest='threads', help='Maximum number of threads [default: %d]' % DEFAULT_THREADS, default=DEFAULT_THREADS)
	parser.add_option_group(group)

	group = OptionGroup(parser, "Required Arguments")

	group.add_option('-g', dest='input_file', type="file", help='Assembly graph in .GFA') 
	group.add_option('-r', dest='reference', help='Path to the reference genome (FASTA)')
	parser.add_option_group(group)

	parser.set_usage('Usage: \n' +
					 '' + __file__ + ' [options] -g assembly_graph_file -r reference_file') # -a <assembler_name> [--fasta file_with_graph_edge_sequences]\n'
					 # '2) ' + __file__ + ' [options] -a <assembler_name> -i <assembler_output_dir>\tThis option is supported for %s assemblers only.' % ', '.join(SUPPORTED_ASSEMBLERS))

	opts, args = parser.parse_args()

	if (not opts.input_file):
		print('ERROR! You should specify an assembly graph file using the option -g\nUse --help to see the full usage information')
		sys.exit(1)

	if (not opts.reference):
		print('ERROR! You should specify a reference file using the option -r\nUse --help to see the full usage information')
		sys.exit(1)

	if not os.path.exists(opts.output_dir):
		os.makedirs(opts.output_dir)

	start_parse_g = timer()
	#dict_edges, contig_edges, edges_fpath 
	graph = parse_assembler_output(opts.input_file, opts.output_dir, opts.min_edge_len)	
	
	end_parse_g = timer()
	timer_list.append(["Graph parser", start_parse_g, end_parse_g])

	json_output_dirpath = join(opts.output_dir, "data")
	if not os.path.exists(json_output_dirpath):
		os.makedirs(json_output_dirpath)

	if opts.reference:
		reference = parse_reference(opts.reference)

	start_Graphaligner = timer()
	if opts.reference and opts.input_file:
		print("Aligning reference to the assembly graph\n")
		alignments_fpath = run_graph_alignment(opts.input_file, opts.reference, json_output_dirpath, opts.threads)
		print("Done\n")
	end_Graphaligner = timer()
	timer_list.append(["GraphAligner", start_Graphaligner, end_Graphaligner])

	start_parse_a = timer()
	alignments = parse_json(alignments_fpath)
	paths = defaultdict(list)
	NODES = defaultdict(list)
	
	try:
		i_al = range(1,len(alignments)+1)
		input_data = tuple(zip(list(i_al), alignments))
		pool = Pool()
		out_info = pool.map(partial(align_info2, graph = graph, ref = reference, n_align = len(alignments)), input_data)
		for path in out_info:
			nodes = [(path[0][j][0], j) for j in range(0,len(path[0]))]
			for j in range(0, len(path[0])):
				NODES[path[2]] = path[0]
			paths[path[2]] = [node[0] for node in nodes] 
			for node in nodes:
				graph.nodes[node[0]]["pos_mapped"][path[2]].append(path[0][node[1]][1:])
			reference[path[3]]["pos_mapped"][path[2]] = path[1]
		print("Finish parsing")
	except:
		i_al = range(1,len(alignments)+1)
		input_data = tuple(zip(list(i_al), alignments))
		out_info = [align_info2(inp, graph, reference, len(alignments)) for inp in input_data]
		for path in out_info:
			nodes = [(path[0][j][0], j) for j in range(0,len(path[0]))]
			paths[path[2]] = [node[0] for node in nodes] 
			for node in nodes:
				graph.nodes[node[0]]["pos_mapped"][path[2]].append(path[0][node[1]][1:])
			reference[path[3]]["pos_mapped"][path[2]] = path[1]
		print("Finish parsing")

	end_parse_a = timer()
	timer_list.append(["Alignments parser", start_parse_a, end_parse_a])

	start_stats = timer()
	Graph_Stats = []
	count_size = 0 		# >=500
	count_size2 = 0 	# >=1000
	count_indegree = 0 	
	count_outdegree = 0  
	longest = max([graph.nodes[node]["size"] for node in list(graph.nodes())])
	for node in list(graph.nodes()):
		if graph.nodes[node]["size"] >= 500:
			count_size += 1
			if graph.nodes[node]["size"] >= 1000:
				count_size2 += 1
		if graph.in_degree(node) >= 2:
			count_indegree += 1
		if graph.out_degree(node) >= 2:
			count_outdegree += 1
	Graph_Stats.append(("# of nodes", "{}".format(len(graph))))
	Graph_Stats.append(("# w/ len >= 500 bp", "{}".format(count_size)))
	Graph_Stats.append(("# w/ len >= 1000 bp", "{}".format(count_size2)))
	Graph_Stats.append(("# w/ in-degree >= 2", "{}".format(count_indegree)))
	Graph_Stats.append(("# w/ out-degree >= 2", "{}".format(count_outdegree)))
	Graph_Stats.append(("Longest node", "{}".format(longest)))

	print("\n# NODES = {}\n# with length >= 500bp = {}\n# with in-degree >= 2 = {}\n# with out-degree >= 2 = {}\nLongest node = {}\n".format(len(graph),count_size,count_indegree,count_outdegree, longest))
	#print("ID\tLEN\tU_LEN\tMAPD\tCOV\tPOS_MAPD")
	gra_mapped = 0
	Nbases_graph = 0
	GC_con = 0
	for node in list(graph.nodes()):
		node_gc = 0
		node_mapped = 0
		for base in graph.nodes[node]["seq"]:
			if base in ['G','g','C','c']:
				node_gc += 1
		node_gc = float(node_gc)/graph.nodes[node]["size"]
		for path in graph.nodes[node]["pos_mapped"]:
			for pos in graph.nodes[node]["pos_mapped"][path]:
				overlap_operations = re.split('(\d+)', pos[2].strip())
				for i in range(0, len(overlap_operations) - 1):
					if overlap_operations[i+1] == '=':
						node_mapped += int(overlap_operations[i])
		cov_node = float(node_mapped)/graph.nodes[node]["size"] #["unique_size"]
		graph.nodes[node]["coverage"] = cov_node
		graph.nodes[node]["GC"] = node_gc
		#print("{0}\t{1}\t{2}\t{3}\t{4:.2f}\t{5}".format(node, graph.nodes[node]["size"], graph.nodes[node]["unique_size"], node_mapped, cov_node, dict(graph.nodes[node]["pos_mapped"])))
		gra_mapped += node_mapped
		Nbases_graph += graph.nodes[node]["size"] #["unique_size"]
		GC_con += node_gc
	cov_graph = float(gra_mapped)/Nbases_graph
	GC_con = float(GC_con)/Nbases_graph
	print("GRAPH\tMAPD\tLEN\tCOV\tGC\n\t{0}\t{1}\t{2:.2f}\t{3:.2f}\n".format(gra_mapped, Nbases_graph, cov_graph, GC_con))

	Graph_Stats.append(("Total length", "{}".format(Nbases_graph)))
	Graph_Stats.append(("GC content [%]", "{:.2f}".format(GC_con*100)))

	print("# EDGES = {}\n".format(len(list(graph.edges()))))
	#for edge in list(graph.edges(data=True)):
	#	print(edge)
	print("# CYCLES = {}\n".format(len(list(nx.simple_cycles(graph)))))
	#for cycle in list(nx.simple_cycles(g)):
	#	print(cycle)

	Graph_Stats.append(("# of edges", "{}".format(len(list(graph.edges())))))
	Graph_Stats.append(("# of cycles", "{}".format(len(list(nx.simple_cycles(graph))))))
	
	count_outdegree = 0
	for node in list(graph.nodes()):
		count_outdegree += graph.out_degree(node)
	branching_factor = float(count_outdegree)/len(graph)
	print("Avg. Branching Factor = {:.2f}\n".format(branching_factor))

	Graph_Stats.append(("Avg. branching factor", "{:.2f}".format(branching_factor)))

	Nbases_ref = sum([reference[ref]['len'] for ref in reference])
	nodes = list(graph.nodes(data = 'size'))
	nodes = sorted(nodes, key = itemgetter(1))
	flag0 = False
	flag1 = False
	NG50 = nodes[-1][1]
	LG50 = len(nodes)
	for i in range(0, len(nodes)):
		tmp = sum([nodes[j][1] for j in range(i,len(nodes))])
		if not flag0 and tmp <= Nbases_graph/2:
			N50 = nodes[i-1][1]
			L50 = len(nodes)-(i-1)
			flag0 = True
			if flag0 and flag1:
				break
		if tmp <= Nbases_ref/2 and not flag1:
			NG50 = nodes[i-1][1]
			LG50 = len(nodes)-(i-1)
			flag1 = True
			if flag0 and flag1:
				break

	print("N50 = {}\nL50 = {}\nNG50 = {}\nLG50 = {}\n".format(N50, L50, NG50, LG50))

	Graph_Stats.append(("N50", "{}".format(N50)))
	Graph_Stats.append(("L50", "{}".format(L50)))
	Graph_Stats.append(("NG50", "{}".format(NG50)))
	Graph_Stats.append(("LG50", "{}".format(LG50)))

	unaligned_nodes = 0
	unaligned_len = 0
	for node in list(graph.nodes()):
		if not graph.nodes[node]["pos_mapped"]:
			unaligned_nodes += 1
			unaligned_len += graph.nodes[node]["size"]
	print("# fully unaligned nodes = {}/{} (cumulative length = {})\n".format(unaligned_nodes,len(graph),unaligned_len))		
	print("% of aligned nodes = {:.2f}%\n".format(float(len(graph)-unaligned_nodes)/len(graph)*100))
	
	Graph_Stats.append(("% of aligned bases", "{:.2f}".format(cov_graph*100)))
	Graph_Stats.append(("% of aligned nodes", "{:.2f}".format(float(len(graph)-unaligned_nodes)/len(graph)*100)))
	Graph_Stats.append(("# of fully unaligned nodes", "{}".format(unaligned_nodes)))
	Graph_Stats.append(("Fully unaligned length", "{}".format(unaligned_len)))

	Genome_Stats = {}
	

	print("\nREFERENCE\tLEN\tMAPD\tCOV\tGC")
	for ref in reference:
		SNPs = 0
		mismatches = 0
		indels_len = 0
		indels = 0
		short_INDELs = 0
		long_INDELs = 0
		GC_con = 0
		for base in reference[ref]["seq"]:
			if base in ['g','G','c','C']:			
				GC_con += 1	
		GC_con = float(GC_con)/reference[ref]["len"]
		pos_mapped = reference[ref]['pos_mapped']
		ref_mapped = 0
		for path in pos_mapped:
			for pos in pos_mapped[path]:
				overlap_operations = re.split('(\d+)', pos[2].strip())
				for i in range(0, len(overlap_operations) - 1):
					if overlap_operations[i+1] == '=':
						ref_mapped += int(overlap_operations[i])
					elif overlap_operations[i+1] == 'X':
						mismatches += int(overlap_operations[i])
						if int(overlap_operations[i]) == 1:
							SNPs += 1
					elif overlap_operations[i+1] == 'I' or overlap_operations[i+1] == 'D':
						indels_len += int(overlap_operations[i])
						indels += 1
						if int(overlap_operations[i]) > 5:
							long_INDELs += 1
						else:
							short_INDELs += 1

		cov_ref = float(ref_mapped)/reference[ref]['len']
		print("{0}\t{1}\t{2}\t{3:.2f}\t{4:.2f}".format(ref, reference[ref]['len'], ref_mapped, cov_ref, GC_con))
		#print("PATH\tPOS_MAPD")
		#for path in reference[ref]['pos_mapped']:
		#	for pos in reference[ref]['pos_mapped'][path]:
		#		print("{}\t{}".format(path,pos))
	
		print("\n#mismatches = {}\n#SNPs = {}\n\n#indels = {}\nIndels length = {}\n#short = {}\n#long = {}\n".format(mismatches, SNPs, indels, indels_len, short_INDELs, long_INDELs))
		Genome_Stats[ref] = []
		Genome_Stats[ref].append(("Reference length", str(reference[ref]['len'])))
		Genome_Stats[ref].append(("Reference GC [%]", "{:.2f}".format(GC_con*100)))
		Genome_Stats[ref].append(("Genome fraction [%]", "{:.2f}".format(cov_ref*100)))
		Genome_Stats[ref].append(("# of mismatches", str(mismatches)))
		Genome_Stats[ref].append(("# mismatches per 100kbp", "{:.2f}".format(float(mismatches)/reference[ref]['len']*100000)))
		Genome_Stats[ref].append(("# of SNPs", str(SNPs)))
		Genome_Stats[ref].append(("# of indels", str(indels)))
		Genome_Stats[ref].append(("# indels per 100kbp", "{:.2f}".format(float(indels)/reference[ref]['len']*100000)))
		Genome_Stats[ref].append(("# indels <= 5bp", str(short_INDELs)))
		Genome_Stats[ref].append(("# indels > 5bp", str(long_INDELs)))
		Genome_Stats[ref].append(("Indels length", str(indels_len)))

	end_stats = timer()
	timer_list.append(["Statistics", start_stats, end_stats])

	start_best = timer()
	Best_alignments = best_alignments_selection(graph, reference)
	end_best = timer()
	timer_list.append(["Best set selection", start_best, end_best])

	for ref in Best_alignments:
		print(ref)
		best_ref_set = Best_alignments[ref]['reference'] 
		print('Set score: ', best_ref_set[1])
		for al in best_ref_set[0]:
			print(al)
		print("\n")

	start_visualization = timer()
	g = gt.Graph()
	v_ID = g.new_vertex_property("string")
	v_GC = g.new_vertex_property("float")
	v_len = g.new_vertex_property("int")
	v_uni = g.new_vertex_property("int")
	v_cov = g.new_vertex_property("float")
	v_seq = g.new_vertex_property("string")
	v_pos = g.new_vertex_property("object")
	v_col = g.new_vertex_property("vector<float>")
	v_ID_plot = g.new_vertex_property("string")
	v_len_plot = g.new_vertex_property("string")
	v_cov_plot = g.new_vertex_property("string")

	v_by_ID = {}

	for node in list(graph.nodes()):
		v = g.add_vertex()
		v_ID[v] = node
		v_ID_plot[v] = "Node ID = " + node
		v_len[v] = graph.nodes[node]["size"]
		v_len_plot[v] = "Length = " + str(v_len[v])
		v_cov[v] = graph.nodes[node]["coverage"]
		v_cov_plot[v] = "Coverage = " + "{:.2f}".format(v_cov[v])
		v_col[v] = [0.0, 0.0, 0.0, 0.0]
		v_seq[v] = graph.nodes[node]["seq"]
		v_pos[v] = graph.nodes[node]["pos_mapped"]
		v_uni[v] = graph.nodes[node]["unique_size"]
		v_GC[v] = graph.nodes[node]["GC"]
		v_by_ID[node] = v

	g.vertex_properties["ID"] = v_ID
	g.vertex_properties["GC"] = v_GC
	g.vertex_properties["length"] = v_len
	g.vertex_properties["u_length"] = v_uni
	g.vertex_properties["coverage"] = v_cov
	g.vertex_properties["sequence"] = v_seq
	g.vertex_properties["pos_mapped"] = v_pos

	e_FO = g.new_edge_property("string")
	e_TO = g.new_edge_property("string")
	e_over = g.new_edge_property("string")
	e_wid = g.new_edge_property("float")

	for edge in list(graph.edges()):
		v1 = v_by_ID[edge[0]]
		v2 = v_by_ID[edge[1]]
		e = g.add_edge(v1, v2)
		e_FO[e] = graph.edges[edge[0],edge[1]]["from_ori"]
		e_TO[e] = graph.edges[edge[0],edge[1]]["to_ori"]
		e_over[e] = graph.edges[edge[0],edge[1]]["overlap"]
		e_wid[e] = 1.0

	g.edge_properties["from_ori"] = e_FO
	g.edge_properties["to_ori"] = e_TO
	g.edge_properties["overlap"] = e_over

	win = GAGET_Window(g, geometry=(1500, 1000), gr_stats=Graph_Stats, ge_stats=Genome_Stats, paths=paths, al_info=NODES, v_ID=v_by_ID,
	 	edge_pen_width=e_wid, display_props=[v_ID_plot], display_props_size=12.) #display_props=[v_ID_plot, v_len_plot, v_cov_plot], display_props_size=14.
	
	win.connect("destroy", Gtk.main_quit)
	win.maximize()
	win.show_all()
	end_visualization = timer()
	timer_list.append(["UI generation", start_visualization, end_visualization])
	end_total = timer()
	timer_list.append(["Total time", start_total, end_total])

	for task in timer_list:
		print("{} = {} s\n".format(task[0], task[2]-task[1]))
	Gtk.main()

if __name__ == '__main__':
	main()