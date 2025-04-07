
///////////////////////////////////////////////////////////////////////////////////////////////
////////  PlotTypeHM_*  ///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

// NOTE: heatmaps are rendered on an html canvas and maintain state, so their get_html_targ_id_str() must be ST-dependent... explain this a bit more
class PlotTypeHM extends PlotType {

    constructor() {

	super();

	if (!this.get_color_from_spin_val) throw new Error("Derived PlotTypeHM must define get_color_from_spin_val()");
	this.tile_dim = undefined;  // set in determine_tile_canv_dims()
	this.canv_dim = undefined;  // set in determine_tile_canv_dims()
    }

    get_plot_width() {
	return this.canv_dim;
    }

    get_plot_height() {
	return this.canv_dim;
    }

    determine_tile_canv_dims(N) {
	let rough_tile_dim = Math.floor(PlotType.square_plot_width / N);  // "target" a plot dim that is just under PlotType.square_plot_width
	this.tile_dim = parseInt(Math.max(rough_tile_dim, 1.0));  // require at least 1 pixel per tile (which may put plot dim above target!)
	this.canv_dim = this.tile_dim * N;
    }

    setup_canvas_etc() {
	$("#" + this.get_html_targ_id_str()).attr("width", this.canv_dim);   // NOTE: HM canvas plot_targets use **attr* not *css* for h/w
	$("#" + this.get_html_targ_id_str()).attr("height", this.canv_dim);  // NOTE: HM canvas plot_targets use **attr* not *css* for h/w
	this.cc = document.getElementById(this.get_html_targ_id_str()).getContext("2d", { alpha: false });  // cc = canvas context, for plotting calls
	let init_sa = this.trj.get_x(this.trj.t_0).spins;  // spin array associated with initial time step t_0 used to initialize canvas
	this.draw_entire_canv_from_spin_arr(init_sa);
	this.last_t_displayed = this.trj.t_0;
    }
    
    set_pixel(x, y, v) {

	let xc = y * this.tile_dim;  // xc = x on canvas; NOTE THE TRANSPOSE OPERATION y --> xc and x --> yc
	let yc = x * this.tile_dim;  // yc = y on canvas; NOTE THE TRANSPOSE OPERATION y --> xc and x --> yc
	this.cc.fillStyle = this.get_color_from_spin_val(v);
	this.cc.fillRect(xc, yc, this.tile_dim, this.tile_dim);
    }

    draw_entire_canv_from_spin_arr(sa) {

	for (let x = 0; x < sa.shape[0]; x++) {
	    for (let y = 0; y < sa.shape[1]; y++) {
		this.set_pixel(x, y, sa.get(x, y));
	    }
	}
    }

    update_canvas(t) {

	let moving_fwd = (this.last_t_displayed < t);
	if (moving_fwd) {
	    for (let s = this.last_t_displayed; s < t; s++) {
		let nt = this.trj.get_x(s).next_trans;  // nt = next transition object (for coords at time s)
		if (nt.move_occurred) {
		    this.set_pixel(nt.x, nt.y, nt.new_val);
		}
	    }
	} else {  // moving backward (or at t == this.last_t_displayed, in which case for loop will execute 0 times)
	    for (let s = this.last_t_displayed; s > t; s--) {
		let pt = this.trj.get_x(s).prev_trans;  // pt = previous transition object (for coords at time s)
		if (pt.move_occurred) {
		    this.set_pixel(pt.x, pt.y, pt.old_val);
		}
	    }
	}
	this.last_t_displayed = t;  // update for next call
    }

    get_ext_x_axis_lbl_str() {
	return "\\mathrm{spin} \\hspace{0.4em} x \\hspace{0.4em} \\mathrm{index}";
    }

    get_ext_y_axis_lbl_str() {
	return "\\mathrm{spin} \\hspace{0.4em} y \\hspace{0.4em} \\mathrm{index}";
    }

    plot(t) {
	this.update_canvas(t);
    }
}

class PlotTypeHM_IS extends PlotTypeHM {

    constructor(trj) {

	super();

	this.trj = trj;
	this.determine_tile_canv_dims(this.trj.mc.N);
	this.setup_canvas_etc();
    }

    get_html_targ_id_str() {
	return "plot_HM_IS";  // heatmaps are rendered on an html canvas and maintain state, so must be ST-dependent
    }

    get_color_from_spin_val(v) {  // translates a numeric spin value into an html color string

	return ((v == 0) ? "hsl(29, 85%, 44%)" : "hsl(52, 100%, 51%)");  // dark orange hsl(29, 85%, 44%), bright yellow hsl(52, 100%, 51%)
    }
}

class PlotTypeHM_XY extends PlotTypeHM {

    constructor(trj) {

	super();

	this.trj = trj;
	this.determine_tile_canv_dims(this.trj.mc.N);
	this.setup_canvas_etc();
    }

    get_html_targ_id_str() {
	return "plot_HM_XY";  // heatmaps are rendered on an html canvas and maintain state, so must be ST-dependent
    }

    get_color_from_spin_val(v) {  // translates a numeric spin value into an html color string

	let hue_val = v * 360.0;  // convert fraction of revolution --> degrees
	let hsl_str = "hsl(" + hue_val + ", 55%, 50%)";
	return hsl_str;
    }
}
