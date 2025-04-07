
///////////////////////////////////////////////////////////////////////////////////////////////
////////  RW = Random Walk (from SP = Stochastic Processes)  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

class ModelCalc_RW extends ModelCalc_SP {

    constructor(N) {
	super();
	this.N = N;
    }

    get_P_step_R(p, x) {  // probability (conditioned on value of x) of x increasing
	return (x == this.N) ? 0.0 : p.r;  // (of course, in RW, there is no x-dependence)
    }

    get_P_step_L(p, x) {  // probability (conditioned on value of x) of x decreasing
	return (x == 0) ? 0.0 : p.l;  // (of course, in RW, there is no x-dependence)
    }
}

class Params_RW extends Params {

    static l = undefined;  // = new UINI_float(this, "UI_P_SP_RW_l", true);  assignment occurs in UserInterface(); see discussion there
    static r = undefined;  // = new UINI_float(this, "UI_P_SP_RW_r", true);  assignment occurs in UserInterface(); see discussion there

    constructor(l_val, r_val) {
	super();
	this.l = l_val;
	this.r = r_val;
    }

    push_vals_to_UI() {
	Params_RW.l.push_to_UI(this.l);
	Params_RW.r.push_to_UI(this.r);
    }

    get_info_str() {
	return "l, r = " + this.l + ", " + this.r;
    }
}

class Coords_RW extends Coords_SP_finite {

    static N = undefined;  // = new UINI_int(this, "UI_P_SP_RW_N", false);  assignment occurs in UserInterface(); see discussion there
    static x_0 = undefined;  // = new UINI_int(this, "UI_P_SP_RW_x_0", false);  assignment occurs in UserInterface(); see discussion there

    constructor(...args) {  // "..." is Javascript spread operator
	super(...args);
    }
}

class Trajectory_RW extends Trajectory_SP {

    constructor(sim) {
	super(sim);
    }

    gmc() {  // gmc = get ModelCalc object
	return new ModelCalc_RW(Coords_RW.N.v);
    }

    gp() {  // gp = get Params object
	return new Params_RW(Params_RW.l.v, Params_RW.r.v);
    }

    gc_ic(mc) {  // gc_ic = get Coords, initial condition
	return new Coords_RW(mc, [ Coords_RW.N.v, Coords_RW.x_0.v ]);
    }

    gc_nv(mc, p, c_prev) {  // gc_nv = get Coords, new value
	return new Coords_RW(mc, p, c_prev, [ Coords_RW.N.v ]);
    }

    get_max_num_t_steps() {
	return Trajectory.DEFAULT_MAX_NUM_T_STEPS
    }
}
