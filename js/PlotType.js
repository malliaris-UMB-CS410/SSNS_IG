///////////////////////////////////////////////////////////////////////////////////////////////
////////  PlotType  ///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

//
// * flot plot options are often static, then possibly copied and "injected" with data in "flot_obj['data']" property
//
class PlotType {

    static registered_PTs = ["XT", "HX", "HM", "PP"];
    static square_plot_width;   // may eventually tweak method of assigning numbers, so keep width/height separate for "square" plots
    static square_plot_height;  // may eventually tweak method of assigning numbers, so keep width/height separate for "square" plots
    static non_square_plot_width;
    static non_square_plot_height;
    static set_static_member_vals(spw, sph, nspw, nsph) {
	PlotType.square_plot_width = spw;
	PlotType.square_plot_height = sph;
	PlotType.non_square_plot_width = nspw;
	PlotType.non_square_plot_height = nsph;
    }

    constructor() {

	if (!this.plot) throw new Error("Derived PlotType must define plot()");
	if (!this.get_html_targ_id_str) throw new Error("Derived PlotType must define get_html_targ_id_str()");
	if (!this.get_plot_width) throw new Error("Derived PlotType must define get_plot_width()");
	if (!this.get_plot_height) throw new Error("Derived PlotType must define get_plot_height()");
	if (!this.get_ext_x_axis_lbl_str) throw new Error("Derived PlotType must define get_ext_x_axis_lbl_str()");
	if (!this.get_ext_y_axis_lbl_str) throw new Error("Derived PlotType must define get_ext_y_axis_lbl_str()");
    }

    apply_plot_dims_to_targ_css() {
	$("#" + this.get_html_targ_id_str()).css("width", this.get_plot_width());
	$("#" + this.get_html_targ_id_str()).css("height", this.get_plot_height());
    }

    switch_plot_and_update_ext_axes() {

	// update plot itself; this will sometimes lead to unnecessary hide/show, but no harm done
	$(".plot_target:not(:hidden)").hide();
	$("#" + this.get_html_targ_id_str()).show();
	if (this.get_html_targ_id_str() == "plot_flot") {
	    this.apply_plot_dims_to_targ_css();
	}

	// ... and label element dimensions, which depend on whether plot is square
	let pw = this.get_plot_width();  // pw = plot width; for convenience
	let ph = this.get_plot_height();  // ph = plot height; for convenience
	$("#plot_x_axis_lbl").css("width", pw);
	$("#plot_y_axis_lbl").css("width", ph);  // y axis div **width** assigned **height h_to_use** since it is rotated 90 degrees
	let navbar_height = parseInt($("#navbar_div_container").css("height"));
	let y_lbl_xfrm_X = ph + navbar_height - 8;  // after rotation, xfrm_X is up-down movement; -8 ~corrects for flot x ticklabel height
	let y_lbl_xfrm_str = "rotate(-90deg) translateX(-" + y_lbl_xfrm_X + "px)";  // NOTE: translation is in negative direction
	$("#plot_y_axis_lbl").css("transform", y_lbl_xfrm_str);  // apply the transform that we've measured out

	// update labels...
	CU.sk("plot_x_axis_lbl", this.get_ext_x_axis_lbl_str());
	CU.sk("plot_y_axis_lbl", this.get_ext_y_axis_lbl_str());
    }

    // take flot line plot data [ [x0, y0], [x1, y1], ... ], and return histogram data; x coordinates all shifted
    // left by 1/2 and an extra pair appended with x right-shifted and same y as the former last pair
    // *** ASSUMES x VALUES ARE CONSECUTIVE ***
    get_hist_data_flot(line_plot_data) {
	let hist_data = [];
	for (let i = 0; i < line_plot_data.length; i++) {
	    let new_pair = [ line_plot_data[i][0] - 0.5, line_plot_data[i][1] ];
	    hist_data.push(new_pair);
	}
	let extra_pair = [ hist_data.at(-1)[0] + 1.0, hist_data.at(-1)[1] ];
	hist_data.push(extra_pair);
	return hist_data;
    }

    set_xlim_flot(flot_gen_opts, left, right) {  // set x limits for flot_gen_opts obj passed in
	if ( ! ("xaxis" in flot_gen_opts)) {
	    flot_gen_opts["xaxis"] = {};
	}
	flot_gen_opts["xaxis"]["autoScale"] = "none";
	flot_gen_opts["xaxis"]["min"] = left;
	flot_gen_opts["xaxis"]["max"] = right;
    }

    set_xticks_intgr_flot(flot_gen_opts, L, R) {  // set x ticks for flot_gen_opts obj passed in to be some subset of the integers
	if ( ! ("xaxis" in flot_gen_opts)) {
	    flot_gen_opts["xaxis"] = {};
	}
	flot_gen_opts["xaxis"]["ticks"] = Array.from(Array(R - L + 1), (e,i) => i + L);  // like Python/NumPy's range/arange methods
    }

    set_ylim_flot(flot_gen_opts, bottom, top) {  // set y limits for flot_gen_opts obj passed in
	if ( ! ("yaxis" in flot_gen_opts)) {
	    flot_gen_opts["yaxis"] = {};
	}
	flot_gen_opts["yaxis"]["autoScale"] = "none";
	flot_gen_opts["yaxis"]["min"] = bottom;
	flot_gen_opts["yaxis"]["max"] = top;
    }

    set_yticks_intgr_flot(flot_gen_opts, B, T) {  // set y ticks for flot_gen_opts obj passed in to be some subset of the integers
	if ( ! ("yaxis" in flot_gen_opts)) {
	    flot_gen_opts["yaxis"] = {};
	}
	flot_gen_opts["yaxis"]["ticks"] = Array.from(Array(T - B + 1), (e,i) => i + B);  // like Python/NumPy's range/arange methods
    }

    add_vert_line_flot(flot_gen_opts, x_loc, lw, clr_str) {  // add vertical line to flot_gen_opts obj passed in
	if ( ! ("grid" in flot_gen_opts)) {
	    flot_gen_opts["grid"] = {};
	}
	flot_gen_opts["grid"]["markings"] = [ { xaxis: { from: x_loc, to: x_loc }, lineWidth: lw, color: clr_str } ];
    }    

    remove_vert_line_flot(flot_gen_opts) {  // remove vertical line (and all other grid entries) from flot_gen_opts obj passed in
	if ("grid" in flot_gen_opts) {
	    flot_gen_opts["grid"]["markings"] = {};
	}
    }
}
