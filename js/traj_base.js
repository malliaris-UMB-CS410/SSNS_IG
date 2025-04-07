///////////////////////////////////////////////////////////////////////////////////////////////
////////  SSNS          ///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////////  Javascript code for Simple Stochastic and Nonlinear Simulator (SSNS)  ///////////////
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
////////  ABSTRACT BASE CLASSES  //////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////////  these indicate requirements/relationships and are never instantiated  ///////////////
///////////////////////////////////////////////////////////////////////////////////////////////

class ModelCalc {

    constructor() {
	if (!this.model_is_stoch) throw new Error("Derived ModelCalc must define model_is_stoch()");
    }
}

class Params {

    constructor() {
	if (!this.push_vals_to_UI) throw new Error("Derived Params must define push_vals_to_UI()");
	if (!this.get_info_str) throw new Error("Derived Params must define get_info_str()");
    }
}

// UPDATE THIS TO INCLUDE Coords_SP case which is 2,4 arguments instead of 1,3
//constructor() {  // Coords_*() will be called with either:
//constructor(mc) {  // 1 argument, mc, in which case initial condition will be used to create coordinate data
//constructor(mc, p, c_prev) {  // or 3 arguments, in which case data will be calculated by mc using p and c_prev as inputs
//      // SAMPLE LINES:
//	    super();
//	    let mc = arguments[0];  (first arg in both cases)
//	    if (arguments.length == 1) {           -----> coordinate values from initial condition
//	    } else if (arguments.length == 3) {    -----> coordinate values update equation
class Coords {

    constructor() {  // see discussion of # args at definition of abstract Coords()

	if ((arguments.length == 2) || (arguments.length == 4)) {

	    this.constructing_init_cond = (arguments.length == 2);

	    this.mc = arguments[0];
	    this.extra_args = this.constructing_init_cond ? arguments[1] : arguments[3];
	    if ( ! this.constructing_init_cond) {
		this.p = arguments[1];
		this.c_prev = arguments[2];
	    }
	    
	} else throw new Error("ERROR 490122: Invalid number of arguments to Coords(): " + arguments.length + ".  (Should be 2 or 4.)  Exiting...");
    }
}


// TrajSeg = Trajectory Segment, the "monomer" from which a "polymeric" Trajectory class is built
//
// * of the 3 pieces of information -- t_first, t_last, and number of t steps -- we could get by with only storing 1...
// * we choose to store both t_first and t_last (even though there's redundant info) to make t navigation easier
// * number of t steps is always available via the length of various arrays, or from # = t_last - t_first + 1
//
class TrajSeg {

    constructor(mc, t_first, p, c0, gc_nv) {
	this.t_first = t_first;
	this.t_last = t_first;  // since, initially, there's only one time step in the segment
	this.mc = mc;
	this.p = p;
	this.cv = [ c0 ];  // cv = coordinate values array, started with initial coord object c0
	this.gc_nv = gc_nv;  // routine to use for generating further coord objects
    }

    output_info() {
	let output_str = "";
	output_str += "\t\t" + this.p.get_info_str() + ";\t\t" + this.cv.length + " coord obj(s);\t\t";  // raw data
	output_str += "# t,x steps match = " + this.t_x_lengths_match();  // and a check
	console.log(output_str);
    }

    get_num_t_steps() {  // NOTE: this is Delta(t) + 1
	return (this.t_last - this.t_first + 1);
    }

    get_t_arr() {
	return incrspace(this.t_first, this.t_last + 1, 1);
    }

    append_step() {
	let c_new = this.gc_nv(this.mc, this.p, this.cv.at(-1));
	this.cv.push(c_new);
	this.t_last++;	    
    }

    truncate_seg(new_t) {  // "shorten" a TrajSeg, saving some already-assembled data (Coords, rng_states, etc.), in preparation for recording
	this.t_last = new_t;  // (will affect get_num_t_steps() in next step)
	this.cv.splice(this.get_num_t_steps());  // splice(n) retains [0, n - 1] and discards the tail
    }

    t_x_lengths_match() {  // # t steps (calculated via t_first/last) should always match # x steps, i.e., length of coord vals array
	return (this.get_num_t_steps() == this.cv.length);
    }
}


//
// Trajectory -- a class that encapsulates all trajectory information
//
// * composed of one or more TrajSeg's, stored in sequence in an array
// * each TrajSeg contains a Params object, rng states, and relevant variable values, i.e., trajectory info
// * "rewinding" and changing parameter values causes creation/editing of relevant TrajSeg and "re-recording"
//
class Trajectory {

    static DEFAULT_MAX_NUM_T_STEPS = 10001;  // i.e., t_max == 10000 if t_0 = 0

    // unlike ModelCalc, Params, Coords, and TrajSeg, Trajectory does have a non-trivial constructor since there is a lot to set up...
    constructor(sim) {

	// UPDATE THESE LINES TO GENERALLY DESCRIBE REQUIRED METHODS BELOW Trajectory_RW uses ModelCalc_RW as its ModelCalc class
	// Trajectory_RW uses Params_RW as its Params class, and, most commonly, we pull parameter values from the user interface
	// gp = get Params object from UI (see discussion at definition)
	// Trajectory_RW uses Coords_SP as its Coords class; one method returns IC version, the other, new value version
	if (!this.gmc) throw new Error("Derived Trajectory must define gmc()");
	if (!this.gp) throw new Error("Derived Trajectory must define gp()");
	if (!this.gc_ic) throw new Error("Derived Trajectory must define gc_ic()");
	if (!this.gc_nv) throw new Error("Derived Trajectory must define gc_nv()");
	if (!this.get_max_num_t_steps) throw new Error("Derived Trajectory must define get_max_num_t_steps()");  // ST-specific to sensibly limit memory usage

	this.sim = sim;

	// store/determine time boundaries of trajectory
	this.t_0 = this.sim.ui.t_0.v;
	let requested_t_max = this.sim.ui.t_max.v;  // Javascript ensures in "real time" that t_max >= t_0 + 1, so no need to check for that
	let largest_t_max = this.t_0 + this.get_max_num_t_steps() - 1;  // largest ("max") t_max that could be allowed based on sensible memory limits
	this.t_max = Math.min(requested_t_max, largest_t_max);

	// nascent Trajectory includes a single TrajSeg which includes a single time step
	this.segs = new Array();
	this.t = this.t_0;
	this.t_edge = this.t;
	this.csi = 0;  // csi = current segment index

	// create ModelCalc and initial TrajSeg, which incorporates initial condition and default parameter values
	this.mc = this.gmc();
	let init_traj_seg = new TrajSeg(this.mc, this.t_0, this.gp(), this.gc_ic(this.mc), this.gc_nv);
	this.segs.push(init_traj_seg);
    }

    add_seg_at_edge() {
	let c_curr = this.get_x(this.t_edge);
	let p_new = this.gp();
	let t_new = this.t_edge + 1;
	let c_new = this.gc_nv(this.mc, p_new, c_curr)
	let new_traj_seg = new TrajSeg(this.mc, t_new, p_new, c_new, this.gc_nv);  // (last arg stores the function gc_nv() for future use)
	this.segs.push(new_traj_seg);
    }
    
    push_curr_Params_to_UI() {
	this.segs[this.csi].p.push_vals_to_UI();
	this.sim.ui.indicate_params_pushed_traj_to_UI();
    }

    check_integrity() {
	console.log("CHECKING Trajectory OBJECT INTEGRITY...")
	console.log("INFO:\tobject properties:\tcsi, t_0, t_edge, t_max =", this.csi, ",", this.t_0, ",", this.t_edge, ",", this.t_max);
	let csi_t_consis = ( (this.segs[this.csi].t_first <= this.t) && (this.t <= this.segs[this.csi].t_last) );
	let t_last_edge_consis = (this.segs.at(-1).t_last == this.t_edge);
	console.log("\tcsi_t_consis = ", csi_t_consis, "\t\tt_last_edge_consis = ", t_last_edge_consis, "\t\t# segs = ", this.segs.length);
	if (!csi_t_consis) throw new Error("ERROR 823826: mismatch between values of t and csi!  Exiting...");
	if (!t_last_edge_consis) throw new Error("ERROR 823827: mismatch between values of t_last of last seg and t_edge!  Exiting...");
	for (let i = 0; i < this.segs.length; i++) {
	    console.log("\tSEGMENT at index ", i, ":\t\tt runs ", this.segs[i].t_first, " to ", this.segs[i].t_last);
	    if ((i > 0) && (this.segs[i - 1].t_last + 1 != this.segs[i].t_first)) {  // check that t_last of prev seg vs. t_first of current
		throw new Error("ERROR 229687: time mismatch between segments!  Exiting...");
	    }
	    this.segs[i].output_info();  // output ST-specific info
	}
    }

    output_curr_vals() {
	console.log("CURRENT VALUES:")
	console.log("t =", this.t);
	console.log("p:", this.segs.at(-1).p);
	console.log("c:", this.get_cx());
	if (this.get_cx().output) this.get_cx().output();  // if applicable Coords class has .output() method, use it!
    }

    step_backward() {  // move on the trajectory one step backward, i.e., toward t_0
	if ( ! this.at_t_0()) {  // if there's "room" to decrement...
	    this.t--;
	    let move_prev_seg = (this.t < this.segs[this.csi].t_first);  // whether, after t--,  we've moved to the previous TrajSeg in the Trajectory
	    if (move_prev_seg) this.csi--;
	    if (move_prev_seg || this.sim.rs.params_changed) {
		this.push_curr_Params_to_UI();  // MUST BE AFTER csi-- !
	    }
	    return true;  // successfully stepped backward
	} else {
	    return false;  // no step made
	}
    }
    
    step_forward() {  // move on the trajectory one step forward (without recording), i.e., toward t_edge and t_max
	if (! this.at_t_edge()) {  // if there's "room" to increment...
	    this.t++;
	    let move_next_seg = (this.t > this.segs[this.csi].t_last);  // whether, after t++,  we've moved to the next TrajSeg in the Trajectory
	    if (move_next_seg) this.csi++;
	    if (move_next_seg || this.sim.rs.params_changed) {
		this.push_curr_Params_to_UI();  // MUST BE AFTER csi++ !
	    }
	    return true;  // successfully stepped forward
	} else {
	    return false;  // no step made
	}
    }
    
    step_forward_and_record() {  // add a step to the trajectory (only possible in the forward direction)
	if (this.at_t_max()) {  // if there's no "room", i.e., t_edge == t_max, alert the user
	    sim.ui.hv.show_view('HV_E_T_MAX_REACHED');
	    return false;  // no step/recording/changes made
	} else {  // otherwise, we can add a step, and figure out the details...
	    if (this.t < this.t_edge) {  // if we're back from the edge of the recorded trajectory...
		this.truncate_traj();  // ... truncate it in preparation for new segment

		if ((this.mc.model_is_stoch()) && CU.gcb("UI_CTRL_rng_recreate_traj")) {
		    this.set_rng_states_from_edge_vals();  // (add explanation of this functionality somewhere)
		}
	    }
	    if (this.sim.rs.params_changed) {  // if the parameters have been changed, we start a new segment...
		this.add_seg_at_edge();
		this.csi++;  // since we will "move onto" the new segment we've just created
		this.push_curr_Params_to_UI();
	    } else {  // if the parameters are the same, extend current segment...
		this.segs.at(-1).append_step();
	    }
	    this.t_edge++;
	    this.t = this.t_edge;
	    return true;  // successfully recorded one step forward
	}
    }

    truncate_traj() {
	this.segs.splice(this.csi + 1);  // trim unneeded segments from tail of segs array
	this.segs[this.csi].truncate_seg(this.t);  // trim unneeded x values from current segment
	this.t_edge = this.t;
	// this.csi does NOT need to be updated
    }

    jump_to(t_val) {
	this.t = CU.gmv(t_val, this.t_0, this.t_edge);  // gmv = get middle value; if input t_val is out of range, we go to nearest boundary
	let orig_csi = this.csi;
	this.csi = this.get_si(this.t);
	if ((this.csi != orig_csi) || this.sim.rs.params_changed) {
	    this.push_curr_Params_to_UI();
	}
    }

    get_si(t_val) {  // si = segment index
	let i = 0;
	while (t_val > this.segs[i].t_last) i++;
	return i;
    }

    get_xi(t_val) {  // xi = x index
	let si = this.get_si(t_val);
	return (t_val - this.segs[si].t_first);
    }

    get_cxi() {  // cxi = current x index
	return this.get_xi(this.t);
    }

    get_x(t_val) {  // x represents a generalized Coord
	return this.segs[this.get_si(t_val)].cv[this.get_xi(t_val)];
    }

    get_cx() {  // cx = current x; x represents a generalized Coord
	return this.get_x(this.t);
    }

    at_t_0() {
	return (this.t == this.t_0);
    }

    at_t_edge() {
	return (this.t == this.t_edge);
    }

    at_t_max() {
	return (this.t == this.t_max);
    }
    at_enabled_t_stop() {
	return (CU.gcb("UI_CTRL_use_t_stop") && (this.t == this.sim.ui.t_stop.v));  // t_stop is enabled and t == t_stop
    }
}
