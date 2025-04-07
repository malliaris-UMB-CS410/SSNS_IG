
///////////////////////////////////////////////////////////////////////////////////////////////
////////  PlotTypePP_*  ///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

// PlotTypePP = plot type Phase Portrait (loosely interpreted)
class PlotTypePP extends PlotType {

    constructor() {
	
	super();

	if (!this.get_flot_data_series) throw new Error("Derived PlotTypePP must define get_flot_data_series()");
	if (!this.get_flot_gen_opts) throw new Error("Derived PlotTypePP must define get_flot_gen_opts()");
    }

    get_html_targ_id_str() {
	return "plot_flot";  // PP for any ST can be done with flot library on this shared div
    }

    get_plot_width() {
	return PlotType.square_plot_width;
    }

    get_plot_height() {
	return PlotType.square_plot_height;
    }

    plot(t) {
	$.plot($("#" + this.get_html_targ_id_str()), this.get_flot_data_series(t), this.get_flot_gen_opts());
    }
}

// a crude version of the LM bifurcation diagram can be outlined by varying r (not quite a phase portrait -- I know)
class PlotTypePP_LM extends PlotTypePP {

    constructor(trj) {

	super();
	this.trj = trj;

	// flot plot options -- actual data written into property "flot_data_obj['data']" in plot() method
	this.flot_data_opts_reg = {  // for regular points
	    color: "blue",
	    points: {
		fillColor: "blue",//"rgba(40, 40, 200, 0.2)",
		show: true,
		radius: 2,
	    },
	};

	this.flot_data_opts_curr_t = {  // additional outline circle to indicate current t
	    color: "#000000",
	    points: {
		show: true,
		radius: 6,
		fill: 0.3,
		fillColor: null,
	    },
	};

	this.flot_gen_opts = {
	    xaxis: {  // really the r axis
		autoScale: "none",
		min: 0,
		max: 4,
	    },
	    yaxis: {  // really the x axis
		autoScale: "none",
		min: 0,
		max: 1,
	    }
	};
    }

    get_ext_x_axis_lbl_str() {
	return "r";
    }

    get_ext_y_axis_lbl_str() {
	return "x";
    }

    get_plot_data_reg(t) {  // returns what flot documentation call "rawdata" which has format [ [x0, y0], [x1, y1], ... ]

	let max_num_data_pts = 1000;  // eventually put this in UI?
	let num_data_pts_up_to_t = t - this.trj.t_0 + 1;  // how much past trajectory is there to possibly plot?
	let num_data_pts = Math.min(max_num_data_pts, num_data_pts_up_to_t);  // determine how much we will actually plot
	let t_i = t - num_data_pts + 1;
	let data = [];
	for (let s = t_i; s <= t; s++) {
	    let curr_r_val = this.trj.segs[this.trj.get_si(s)].p.r;
	    data.push( [ curr_r_val, this.trj.get_x(s).x ] );
	}
	return data;
    }

    get_plot_data_curr_t(t) {  // returns what flot documentation call "rawdata" which has format [ [x0, y0], [x1, y1], ... ]

	let curr_r_val = this.trj.segs[this.trj.get_si(t)].p.r;
	return [ [ curr_r_val, this.trj.get_x(t).x ] ];
    }

    get_flot_data_series(t) {
	let data_series = [];
	this.flot_data_opts_reg["data"] = this.get_plot_data_reg(t);
	data_series.push(this.flot_data_opts_reg);
	this.flot_data_opts_curr_t["data"] = this.get_plot_data_curr_t(t);
	data_series.push(this.flot_data_opts_curr_t);
	return data_series;
    }

    get_flot_gen_opts() {
	return this.flot_gen_opts;
    }
}

class PlotTypePP_GM extends PlotTypePP {

    constructor(trj) {

	super();
	this.trj = trj;

	// flot plot options -- actual data written into property "flot_data_obj['data']" in plot() method
	this.flot_data_opts_reg = {  // for regular points
	    color: "#c7894e",
	    points: {
		fillColor: "#c7894e",
		show: true,
		radius: 2,
	    },
	};

	this.flot_data_opts_curr_t = {  // additional outline circle to indicate current t
	    color: "#000000",
	    points: {
		show: true,
		radius: 6,
		fill: 0.3,
		fillColor: null,
	    },
	};

	this.flot_gen_opts = {
	    xaxis: {
		autoScale: "none",
		min: -4,
		max: 9,
		transform: function (v) { return -v; },  // make the gingerbread-man easier to recognize by flipping axis
	    },
	    yaxis: {
		autoScale: "none",
		min: -4,
		max: 9,
		transform: function (v) { return -v; },  // make the gingerbread-man easier to recognize by flipping axis
	    }
	};
    }

    get_ext_x_axis_lbl_str() {
	return "\\xleftarrow{\\hspace*{1.0cm}} x";
    }

    get_ext_y_axis_lbl_str() {
	return "\\xleftarrow{\\hspace*{1.0cm}} y";
    }

    get_plot_data_reg(t) {  // returns what flot documentation call "rawdata" which has format [ [x0, y0], [x1, y1], ... ]

	let max_num_data_pts = 2000;  // eventually put this in UI?
	let num_data_pts_up_to_t = t - this.trj.t_0 + 1;  // how much past trajectory is there to possibly plot?
	let num_data_pts = Math.min(max_num_data_pts, num_data_pts_up_to_t);  // determine how much we will actually plot
	let t_i = t - num_data_pts + 1;
	let data = [];
	for (let s = t_i; s <= t; s++) {
	    data.push( [ this.trj.get_x(s).x, this.trj.get_x(s).y ] );
	}
	return data;
    }

    get_plot_data_curr_t(t) {  // returns what flot documentation call "rawdata" which has format [ [x0, y0], [x1, y1], ... ]

	return [ [ this.trj.get_x(t).x, this.trj.get_x(t).y ] ];
    }

    get_flot_data_series(t) {
	let data_series = [];
	this.flot_data_opts_reg["data"] = this.get_plot_data_reg(t);
	data_series.push(this.flot_data_opts_reg);
	this.flot_data_opts_curr_t["data"] = this.get_plot_data_curr_t(t);
	data_series.push(this.flot_data_opts_curr_t);
	return data_series;
    }

    get_flot_gen_opts() {
	return this.flot_gen_opts;
    }
}
