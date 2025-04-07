
///////////////////////////////////////////////////////////////////////////////////////////////
////////  SSNS          ///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////////  Javascript code for Simple Stochastic and Nonlinear Simulator (SSNS)  ///////////////
///////////////////////////////////////////////////////////////////////////////////////////////

//
// UserInterface -- class that represents all interaction -- graphical, parameter input, etc. -- with the user
//
// * 
//
class UserInterface {

    constructor(sim) {

	this.sim = sim;

	this.bsvs = this.get_bsvs();  // bsvs = Bootstrap viewport size
	console.log("INFO:\tBootstrap viewport size =", this.bsvs);

	// UI controls
	this.t_jump = new UINI_int(this, "UI_CTRL_JT", false);
	this.Delta_t = new UINI_int(this, "UI_CTRL_JD", false);
	this.delay = new UINI_int(this, "UI_CTRL_delay", false);
	PlotTypeXT.window_size = new UINI_int(this, "UI_CTRL_window_size", false);
	this.t_0 = new UINI_int(this, "UI_CTRL_t_0", false);
	this.t_max = new UINI_int(this, "UI_CTRL_t_max", false);
	this.t_stop = new UINI_int(this, "UI_CTRL_t_stop", false);
	ModelCalc_Stoch.rng_seed = new UINI_int(this, "UI_CTRL_rng_seed", false);
	Coords_SP.num_IEM = new UINI_int(this, "UI_CTRL_SP_NI", false);
	Coords_SP.num_GEM = new UINI_int(this, "UI_CTRL_SP_NG", false);

	// system-type parameters
	Params_RW.l = new UINI_float(this, "UI_P_SP_RW_l", true);
	Params_RW.r = new UINI_float(this, "UI_P_SP_RW_r", true);
	Coords_RW.N = new UINI_int(this, "UI_P_SP_RW_N", false);
	Coords_RW.x_0 = new UINI_int(this, "UI_P_SP_RW_x_0", false);
	Params_MN.mu = new UINI_float(this, "UI_P_SP_MN_mu", true);
	Params_MN.s = new UINI_float(this, "UI_P_SP_MN_s", true);
	Coords_MN.N = new UINI_int(this, "UI_P_SP_MN_N", false);
	Coords_MN.x_0 = new UINI_int(this, "UI_P_SP_MN_x_0", false);
	Params_CH.alpha = new UINI_float(this, "UI_P_SP_CH_alpha", true);
	Params_CH.beta = new UINI_float(this, "UI_P_SP_CH_beta", true);
	Coords_CH.x_0 = new UINI_int(this, "UI_P_SP_CH_x_0", false);
	Params_IS.T = new UINI_float(this, "UI_P_SM_IS_T", true);
	Coords_IS.N = new UINI_int(this, "UI_P_SM_IS_N", false);
	Params_XY.T = new UINI_float(this, "UI_P_SM_XY_T", true);
	Coords_XY.N = new UINI_int(this, "UI_P_SM_XY_N", false);
	Params_LM.r = new UINI_float(this, "UI_P_ND_LM_r", true);
	Coords_LM.x_0 = new UINI_float(this, "UI_P_ND_LM_x_0", false);
	Coords_GM.x_0 = new UINI_float(this, "UI_P_ND_GM_x_0", false);
	Coords_GM.y_0 = new UINI_float(this, "UI_P_ND_GM_y_0", false);

	// other UI stuff
	this.dropdown_show_vals = {};  // a bool for each dropdown indicating whether is is to be shown or hidden
	this.dropdown_show_vals["UI_parameters"] = true;  // #UI_parameters dropdown default setting
	this.dropdown_show_vals["UI_controls"] = true;  // #UI_controls dropdown default setting
	this.dropdown_show_vals["UI_system_type"] = true;  // #UI_system_type dropdown default setting
	this.push_dd_state_to_UI();  // apply default settings
	this.update_UI_dd_btn_sizes();
	$("#navbar_collapse").on("show.bs.collapse", function ()   {  $("#plot_y_axis_lbl").hide();  });
	$("#navbar_collapse").on("hidden.bs.collapse", function () {  $("#plot_y_axis_lbl").show();  });
	this.hv = new HelpViewer(this.sim, this.bsvs);
    }

    get_bsvs() {  // bsvs = Bootstrap viewport size
	let curr_bsvs_indicator = $(".bsvs_indicator:not(:hidden)");
	return $(curr_bsvs_indicator).attr('data-bsvs');
    }

    process_resize() {
	this.bsvs = this.get_bsvs();
	this.update_UI_dd_btn_sizes();
	console.log("INFO:\tBootstrap viewport size =", this.bsvs);
    }

    push_dd_state_to_UI() {

	for (let [dd_name, dd_bool] of Object.entries(this.dropdown_show_vals)) {
	    if (dd_bool) {
		$("#" + dd_name).collapse("show");
		$("#" + dd_name + "_dd_caret").removeClass("fa-caret-right");
		$("#" + dd_name + "_dd_caret").addClass("fa-caret-down");
	    } else {
		$("#" + dd_name).collapse("hide");
		$("#" + dd_name + "_dd_caret").removeClass("fa-caret-down");
		$("#" + dd_name + "_dd_caret").addClass("fa-caret-right");
	    }
	}
    }

    update_UI_dd_btn_sizes() {

	if (this.bsvs <= 2)  {  // for xs and sm screens...

	    for (let dd_name of Object.keys(this.dropdown_show_vals)) {  // eliminate the first line on the dropdown buttons to conserve space
		let span_id_str = "#" + dd_name + "_dd_btn_txt";
		$(span_id_str).html("");
	    }
	    $("#UI_help_btn_qm").attr("width", "20px");
	    $("#UI_help_btn_qm").attr("height", "30px");

	} else {  // for md, lg, and xl screens...

	    for (let dd_name of Object.keys(this.dropdown_show_vals)) {  // keep the first line in the dropdown buttons
		let span_id_str = "#" + dd_name + "_dd_btn_txt";
		let dd_short_name = $(span_id_str).attr("data-ddsn");  // short name is stored in data-ddsn attribute
		let html_str = "<samp>" + dd_short_name + "</samp><br/>";  // <br/> creates 2 lines on button; <samp> gives monospace look
		$(span_id_str).html(html_str);
	    }
	    $("#UI_help_btn_qm").attr("width", "30px");
	    $("#UI_help_btn_qm").attr("height", "45px");
	}
    }
		
    update_dropdowns(dd_btn_clicked) {  // arg is jQuery element

	let dd_name = dd_btn_clicked.attr("data-target").slice(1);  // slice(1) removes leading "#"
	this.dropdown_show_vals[dd_name] = ! this.dropdown_show_vals[dd_name];  // no matter what, we toggle dropdown clicked
	if (this.bsvs <= 2)  {  // also, for xs and sm screens, only one dropdown should be shown at a time to conserve space, so...
	    for (let dd_name_i of Object.keys(this.dropdown_show_vals)) {  // cycle through all dropdowns
		if (dd_name_i != dd_name) {
		    this.dropdown_show_vals[dd_name_i] = false;  // and close all other than the one clicked
		}
	    }
	}
	this.push_dd_state_to_UI();
    }
    
    update_ST_area_settings() {

	let ST_area_selected = CU.grb("UI_ST_area");
	for (let ST_area of Simulator.registered_ST_areas.values()) {
	    let id_str = "#UI_ST_area_" + ST_area;
	    let form_check_id_str = "#UI_ST_" + ST_area + "_form_check";
	    if (ST_area == ST_area_selected) {
		$(id_str).show();
		$(form_check_id_str).removeClass("SA_unselected");
		$(form_check_id_str).addClass("SA_selected");
	    } else {
		$(id_str).hide();
		$(form_check_id_str).removeClass("SA_selected");
		$(form_check_id_str).addClass("SA_unselected");
	    }   
	}
	console.log("INFO:\tshowing " + ST_area_selected + " system types...");
    }

    update_SYS_dd_ST_background(new_ST) {

	for (let ST_val of Simulator.registered_STs.values()) {  // not terribly efficient to cycle through all ST's, but it's quick and easy
	    let ST_area = Simulator.lookup_area_from_ST[ST_val];
	    let form_check_id_str = "#UI_ST_" + ST_area + "_" + ST_val + "_form_check";
	    if (ST_val == new_ST) {
		$(form_check_id_str).removeClass("ST_unselected");
		$(form_check_id_str).addClass("ST_selected");
	    } else {
		$(form_check_id_str).removeClass("ST_selected");
		$(form_check_id_str).addClass("ST_unselected");
	    }   
	}
    }

    show_ST_params(ST_val_to_show) {  // changing of ST value causes updating of UI parameters section

	for (let ST_val of Simulator.registered_STs.values()) {
	    let id_str = "#UI_P_" + Simulator.lookup_area_from_ST[ST_val] + "_" + ST_val;  // reconstruct, e.g., #UI_P_SM_IS
	    if (ST_val == ST_val_to_show) {
		$(id_str).show();
	    } else {
		$(id_str).hide();
	    }   
	}
    }

    update_plot_type_buttons() {

	for (let PT_val of PlotType.registered_PTs.values()) {  // iterate over ["XT", "HX", "HM", ...]

	    let ui_btn_str = "#UI_PT_" + PT_val + "_btn";  // id for the Bootstrap button

	    // start with a blank slate by clearing all attributes/modifiers for each button (not an error if absent)
	    $(ui_btn_str).attr("onclick", "");
    	    $(ui_btn_str).removeClass("btn-ssns btn-ssns-emphasis");
    	    $(ui_btn_str).removeAttr("disabled");
	    $(ui_btn_str).css("pointer-events", "none");

	    // put each button into 1 of 3 categories -- disabled, being used, and available -- by adding attributes/modifiers
	    if ( ! Object.keys(this.sim.pm.plots[this.sim.ST]).includes(PT_val)) {  // if PT is not available for this ST, disable
    		$(ui_btn_str).addClass("btn-ssns");
    		$(ui_btn_str).attr("disabled", "true");
	    } else if (PT_val == this.sim.PT) {  // if this is the PT currently being used, make it "chosen", but no onclick needed
    		$(ui_btn_str).addClass("btn-ssns-emphasis");
	    } else {  // if this PT is available to switch to, add onclick and enable clicking
    		$(ui_btn_str).addClass("btn-ssns");
		$(ui_btn_str).attr("onclick", "window.sim.switch_PT_within_ST('" + PT_val + "'); ");
		$(ui_btn_str).css("pointer-events", "auto");
	    }
	}
    }

    close_HV_load_ST(ST_to_load) {

	let area_of_ST_to_load = Simulator.lookup_area_from_ST[ST_to_load];
	let ST_to_load_id = "UI_ST_" + area_of_ST_to_load + "_" + ST_to_load;
	let area_of_ST_to_load_id = "UI_ST_" + area_of_ST_to_load;

	$("#md_container").modal("hide");  // close modal
	CU.srb(area_of_ST_to_load_id);  // "manually" set ST area radio button
	this.update_ST_area_settings();  // then call update method that would normally be called from onclick="..."
	CU.srb(ST_to_load_id);  // "manually" set ST radio button
	this.sim.update_ST();  // then call update method that would normally be called from onclick="..."

	this.dropdown_show_vals["UI_parameters"] = false;  // for clarity, make it so that only SYS dropdown is showing
	this.dropdown_show_vals["UI_controls"] = false;
	this.dropdown_show_vals["UI_system_type"] = true;
	this.push_dd_state_to_UI();
    }

    handle_ws_update(str_val) {  // trigger replot if window size changes
	let curr_ws = PlotTypeXT.window_size.v;
	PlotTypeXT.window_size.try_update(str_val); 
	let new_ws = PlotTypeXT.window_size.v;
	if (new_ws != curr_ws) {
	    this.sim.update_plot_and_UI();
	}
    }

    handle_t_0_update(str_val) {  // preserve correct order: t_0 < t_max; duration limit is enforced in Trajectory get_max_num_t_steps(), etc.
	if (this.t_0.yields_valid_value(str_val)) {
	    let req_t_0 = this.t_0.parse_from_str(str_val);
	    let max_t_0 = this.t_max.v - 1;  // enforce t_0 <= t_max - 1
	    let new_t_0 = Math.min(req_t_0, max_t_0);
	    this.t_0.sv(new_t_0);
	} else {
	    console.log("INFO:\thandle_t_0_update() failed to parse a valid number... restoring previous valid value to UI...");
	    this.t_0.sv(this.t_0.v);  // if proposed new value is invalid, push current valid value back to UI
	}
    }

    handle_t_max_update(str_val) {  // preserve correct order: t_0 < t_max; duration limit is enforced in Trajectory get_max_num_t_steps(), etc.
	if (this.t_max.yields_valid_value(str_val)) {
	    let req_t_max = this.t_max.parse_from_str(str_val);
	    let min_t_max = this.t_0.v + 1;  // enforce t_max >= t_0 + 1
	    let new_t_max = Math.max(req_t_max, min_t_max);
	    this.t_max.sv(new_t_max);
	} else {
	    console.log("INFO:\thandle_t_max_update() failed to parse a valid number... restoring previous valid value to UI...");
	    this.t_max.sv(this.t_max.v);  // if proposed new value is invalid, push current valid value back to UI
	}
    }

    handle_RW_l_update(str_val) {  // enforce 0 <= l <= 1 and 0 <= l + r = 1 - s <= 1
	if (Params_RW.l.yields_valid_value(str_val)) {
	    let req_l = Params_RW.l.parse_from_str(str_val);
	    let max_l = 1.0 - Params_RW.r.v;
	    let new_l = CU.gmv(0.0, req_l, max_l);  // gmv = get middle value, i.e., the one that's between the other two numerically
	    Params_RW.l.sv(new_l);
	    if (Params_RW.l.indicate_params_changed) {
		this.indicate_new_param_vals_ready_to_pull_UI_to_traj();
	    }
	} else {
	    console.log("INFO:\thandle_RW_l_update() failed to parse a valid number... restoring previous valid value to UI...");
	    Params_RW.l.sv(Params_RW.l.v);  // if proposed new value is invalid, push current valid value back to UI
	}
    }

    handle_RW_r_update(str_val) {  // enforce 0 <= r <= 1 and 0 <= r + l = 1 - s <= 1
	if (Params_RW.r.yields_valid_value(str_val)) {
	    let req_r = Params_RW.r.parse_from_str(str_val);
	    let max_r = 1.0 - Params_RW.l.v;
	    let new_r = CU.gmv(0.0, req_r, max_r);  // gmv = get middle value, i.e., the one that's between the other two numerically
	    Params_RW.r.sv(new_r);
	    if (Params_RW.r.indicate_params_changed) {
		this.indicate_new_param_vals_ready_to_pull_UI_to_traj();
	    }
	} else {
	    console.log("INFO:\thandle_RW_r_update() failed to parse a valid number... restoring previous valid value to UI...");
	    Params_RW.r.sv(Params_RW.r.v);  // if proposed new value is invalid, push current valid value back to UI
	}
    }

    // general and helper methods
    indicate_new_param_vals_ready_to_pull_UI_to_traj() {
	this.sim.rs.params_changed = true;  // and signal the internal machinery via RunState
	$("#PRMS_icon_and_indicator").css("border-width", "1px");  // helpful to see visual indication of parameter change if not running continuously
    }

    indicate_params_pushed_traj_to_UI() {
	this.sim.rs.params_changed = false;  // i.e., there hasn't been a **user-initiated** change of parameters since we just pushed values from app
	$("#PRMS_icon_and_indicator").css("border-width", "0px");  // visual indication from above has served its purpose
    }

    update_ST_PT_dep_settings() {

	// plot stuff
	this.sim.pm.plots[this.sim.ST][this.sim.PT].switch_plot_and_update_ext_axes();

	// interface stuff
	CU.sk("UI_INDCTR_ST", "\\texttt{" + this.sim.ST + "}");
	this.update_SYS_dd_ST_background(this.sim.ST);
	this.update_plot_type_buttons();
	this.show_ST_params(this.sim.ST);  // show/hide ST-specific parameters
	if (this.sim.trjs[this.sim.ST].mc.model_is_stoch()) {  // show/hide stochastic system machinery
	    $(".UI_CTRL_STOCH").show();
	    if (this.sim.ST_area == "SP") {  // show/hide stochastic process ensemble parameters (which are subset of stochastic system machinery)
		$(".UI_CTRL_SP").show();
	    } else {
		$(".UI_CTRL_SP").hide();
	    }
	} else {
	    $(".UI_CTRL_STOCH").hide();
	}
	if (this.sim.PT == "XT") {  // window size input only applies to PlotTypeXT
	    $("#UI_window_size_input_group").show();
	} else {
	    $("#UI_window_size_input_group").hide();
	}
    }
}
