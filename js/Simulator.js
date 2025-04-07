
///////////////////////////////////////////////////////////////////////////////////////////////
////////  SSNS          ///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////////  Javascript code for Simple Stochastic and Nonlinear Simulator (SSNS)  ///////////////
///////////////////////////////////////////////////////////////////////////////////////////////

//
// Simulator -- outermost container class of which there will be only a single instance
//
class Simulator {

    static registered_ST_areas = [ "SP", "SM", "ND", "FD" ];
    static registered_STs = [];  // filled in populate_ST_info_from_html()
    static lookup_area_from_ST = {};  // filled in populate_ST_info_from_html()

    constructor() {
	this.populate_ST_info_from_html();
	this.rs = new RunState(this);
	this.ui = new UserInterface(this);
	this.trjs = {};
	this.instantiate_all_Trajectories();
	this.pm = new PlottingMachinery(this, this.ui.bsvs);

	this.ST_area = undefined;  // "declaration" only; default is indicated in HTML and loaded in update_ST() call below
	this.ST = undefined;  // "declaration" only; default is indicated in HTML and loaded in update_ST() call below
	this.PT = undefined;  // "declaration" only; value is set in update_ST() call below
	this.ui.update_ST_area_settings();
	this.update_ST();  // contains initial call to this.update_plot_and_UI()
	this.timeout_id = undefined;  // "declaration" only; used with setTimeout() and clearTimeout()


	////////////////////////////  TESTING AREA -- EMPTY WHEN DONE TESTING!  ///////////////////////////////
	//let tempMCS = new ModelCalc_Stoch();
	//tempMCS.test_rngs();
    }

    populate_ST_info_from_html() {
	for (let ST_area of Simulator.registered_ST_areas.values()) {
	    let id_prefix_str = "UI_ST_" + ST_area + "_";
	    let enabled_STs = $("input[id*=" + id_prefix_str + "]:not(:disabled)");
	    let disabled_STs = $("input[id*=" + id_prefix_str + "]:disabled");
	    enabled_STs.each( function( index, input_element ){
		let ST_val = input_element.value;
		Simulator.registered_STs.push(ST_val);
		Simulator.lookup_area_from_ST[ST_val] = ST_area;  // a reverse lookup used, e.g., in ui.show_ST_params()
	    });
	    disabled_STs.each( function( index, input_element ){
		let ST_val = input_element.value;
		let div_id_str = id_prefix_str + ST_val + "_form_check";
		$("#" + div_id_str).addClass("ST_unselected_disabled");
	    });
	}
    }

    instantiate_Trajectory_ST(ST_str) {

	switch(ST_str) {

	case "RW":
	    this.trjs["RW"] = new Trajectory_RW(this);
	    break;
	case "MN":
	    this.trjs["MN"] = new Trajectory_MN(this);
	    break;
	case "CH":
	    this.trjs["CH"] = new Trajectory_CH(this);
	    break;
	case "IS":
	    this.trjs["IS"] = new Trajectory_IS(this);
	    break;
	case "XY":
	    this.trjs["XY"] = new Trajectory_XY(this);
	    break;
	case "LM":
	    this.trjs["LM"] = new Trajectory_LM(this);
	    break;
	case "GM":
	    this.trjs["GM"] = new Trajectory_GM(this);
	    break;
	case "PF":
	    this.trjs["PF"] = new Trajectory_LM(this);  // DUMMY FOR TESTING.... FIX WHEN READY
	    break;
	default:
	    throw new Error("ERROR 298291: Invalid ST_str!  Exiting...");
	}
    }

    instantiate_all_Trajectories() {

	for (let ST_val of Simulator.registered_STs.values()) {
	    this.instantiate_Trajectory_ST(ST_val);
	}
    }

    process_cmd(cmd_str) {

	this.rs.set_NR();  // exit from any running mode; will be reinstated for fresh RB/RF/RC commands below...
	clearTimeout(this.timeout_id);  // clears a queued setTimeout() from RB/RF/RC if it exists

	switch(cmd_str) {
	case "ST":  // ST == System Type, i.e., which model/simulation we are exploring
	    this.update_ST();
	    break;
	case "JB":  // JB == Jump to Beginning of Recorded Trajectory
	    let t_JB = this.trjs[this.ST].segs[0].t_first;
	    this.trjs[this.ST].jump_to(t_JB);
	    this.update_plot_and_UI();
	    break;
	case "SB":  // SB == Step Backward 1 step in time
	    this.trjs[this.ST].step_backward();
	    this.update_plot_and_UI();
	    break;
	case "SF":  // SF == Step Forward 1 step in time
	    this.trjs[this.ST].step_forward();
	    this.update_plot_and_UI();
	    break;
	case "SR":  // SR == Step (forward 1 step in time and) Record
	    this.trjs[this.ST].step_forward_and_record();
	    this.update_plot_and_UI();
	    break;
	case "JE":  // JE == Jump to End of Recorded Trajectory
	    let t_JE = this.trjs[this.ST].segs.at(-1).t_last;
	    this.trjs[this.ST].jump_to(t_JE);
	    this.update_plot_and_UI();
	    break;
	case "RB":  // RB == Run Backward
	    this.rs.set_RB();
	    this.run_backward();  // update_plot_and_UI() called within...
	    break;
	case "PS":  // PS == Pause
	    clearTimeout(this.timeout_id);  // clears a queued setTimeout() from RB/RF/RC if it exists
	    break;  // no update_plot_and_UI() call needed
	case "RF":  // RF == Run Forward
	    this.rs.set_RF();
	    this.run_forward();  // update_plot_and_UI() called within...
	    break;
	case "RC":  // RC == Record
	    this.rs.set_RC();
	    this.run_record();  // update_plot_and_UI() called within...
	    break;
	case "JT":  // JT == Jump to specific Time
	    let t_JT = this.ui.t_jump.v;
	    this.trjs[this.ST].jump_to(t_JT);  // range check/correction takes place within jump_to()
	    this.update_plot_and_UI();
	    break;
	case "JD":  // JD == Jump by specific Delta t
	    let Dt_by_which_to_jump = this.ui.Delta_t.v;
	    let t_JD = this.trjs[this.ST].t + Dt_by_which_to_jump;
	    this.trjs[this.ST].jump_to(t_JD);  // range check/correction takes place within jump_to()
	    this.update_plot_and_UI();
	    break;
	case "RT":  // RT == Reload This ST's Trajectory Object only
	    this.instantiate_Trajectory_ST(this.ST);
	    this.pm.instantiate_PlotTypes_ST(this.ST);
	    this.update_plot_and_UI();
	    break;
	case "RA":  // RA == Reload All ST's Trajectory Objects
	    this.instantiate_all_Trajectories();
	    this.pm.instantiate_all_PlotTypes();
	    this.update_plot_and_UI();
	    break;
	case "CK":
	    //this.trjs[this.ST].check_integrity();
	    //this.trjs[this.ST].output_curr_vals();
	    break;
	default:
	    console.log("ERROR 912352: bad code");
	}
	
    }

    update_ST() {

	this.ST = CU.grb("UI_ST");  // "ST" cmd only indicated a change -- this actually stores new ST value
	this.ST_area = Simulator.lookup_area_from_ST[this.ST];  // store ST_area value (e.g., "ND" for Nonlinear Dynamics) as well
	this.PT = this.pm.active_PT_by_ST[this.ST];  // PT value in use when previously in ST is reinstated
	this.ui.update_ST_PT_dep_settings();
	this.update_plot_and_UI();
	console.log("INFO:\tnow using system", this.ST, "...");
    }

    switch_PT_within_ST(new_PT) {  // a user-click-initiated switch of plot type without a change in system type (ST)

	this.PT = new_PT;  // store new PT value
	this.pm.active_PT_by_ST[this.ST] = this.PT;  // store value to reinstate when ST switches away and then back
	this.ui.update_ST_PT_dep_settings();
	if (!this.rs.is_RN()) {  // if we are not running, show new PT immediately (if we are running, new PT is shown
	    this.update_plot_and_UI();  // next time update_plot_and_UI() is called, which is usually ~instantaneous)
	}
	console.log("INFO:\tnow using plot type", this.PT, "...");
    }

    run_backward() {
	let trj = this.trjs[this.ST];  // for convenience
	// if already running backward AND attempt to step_backward returns success, update things...
	if (this.rs.is_RB() && trj.step_backward()) {
	    this.update_plot_and_UI();
	    // if, in addition, we're not at an enabled t_stop **after** step, keep running backward...
	    if (trj.at_enabled_t_stop()) {
		this.rs.set_NR();
	    } else {
		this.timeout_id = setTimeout(this.run_backward.bind(this), this.ui.delay.v);
	    }
	} else {
	    this.rs.set_NR()
	}
    }
    
    run_forward() {
	let trj = this.trjs[this.ST];  // for convenience
	// if already running forward AND attempt to step_forward returns success, update things...
	if (this.rs.is_RF() && trj.step_forward()) {
	    this.update_plot_and_UI();
	    // if, in addition, we're not at an enabled t_stop **after** step, keep running forward...
	    if (trj.at_enabled_t_stop()) {
		this.rs.set_NR();
	    } else {
		this.timeout_id = setTimeout(this.run_forward.bind(this), this.ui.delay.v);
	    }
	} else {
	    this.rs.set_NR()
	}
    }

    run_record() {  // note that failure in step_forward_and_record() occurs if t == t_max and will result in user notification
	let trj = this.trjs[this.ST];  // for convenience
	// if already recording forward AND attempt to record another step returns success, update things...
	if (this.rs.is_RC() && trj.step_forward_and_record()) {
	    this.update_plot_and_UI();
	    // if, in addition, we're not at an enabled t_stop **after** step, keep recording forward...
	    if (trj.at_enabled_t_stop()) {
		this.rs.set_NR();
	    } else {
		this.timeout_id = setTimeout(this.run_record.bind(this), this.ui.delay.v);
	    }
	} else {
	    this.rs.set_NR()
	}
    }

    update_plot_and_UI() {

	let current_t = this.trjs[this.ST].t;
	CU.sk("UI_INDCTR_t", "t = " + current_t);
	this.pm.plots[this.ST][this.PT].plot(current_t);  // COMMENT OUT THIS LINE TO TEMPORARILY DISABLE PLOTTING
    }
}
