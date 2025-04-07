
///////////////////////////////////////////////////////////////////////////////////////////////
////////  CH = Chemical system (from SP = Stochastic Processes)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

class ModelCalc_CH extends ModelCalc_SP {

    constructor() {
	super();
    }

    get_P_step_R(p, x) {  // probability (conditioned on value of x) of x increasing
	return p.alpha;
    }

    get_P_step_L(p, x) {  // probability (conditioned on value of x) of x decreasing
	return p.beta * x;
    }
}

class Params_CH extends Params {

    static alpha = undefined;  // = new UINI_float(this, "UI_P_SP_CH_alpha", true);  assignment occurs in UserInterface(); see discussion there
    static beta = undefined;  // = new UINI_float(this, "UI_P_SP_CH_beta", true);  assignment occurs in UserInterface(); see discussion there

    constructor(alpha_val, beta_val) {
	super();
	this.alpha = alpha_val;
	this.beta = beta_val;
    }

    push_vals_to_UI() {
	Params_CH.alpha.push_to_UI(this.alpha);
	Params_CH.beta.push_to_UI(this.beta);
    }

    get_info_str() {
	return "alpha, beta = " + this.alpha + ", " + this.beta;

    }
}
class Coords_CH extends Coords_SP_semiinf {

    static x_0 = undefined;  // = new UINI_int(this, "UI_P_SP_CH_x_0", false);  assignment occurs in UserInterface(); see discussion there

    constructor(...args) {  // "..." is Javascript spread operator
	super(...args);
    }
}

class Trajectory_CH extends Trajectory_SP {

    constructor(sim) {
	super(sim);

	// determine how to balance alpha, beta values, max x value, x_0 value, size of "time step"
    }

    gmc() {  // gmc = get ModelCalc object
	return new ModelCalc_CH();
    }

    gp() {  // gp = get Params object
	return new Params_CH(Params_CH.alpha.v, Params_CH.beta.v);
    }

    gc_ic(mc) {  // gc_ic = get Coords, initial condition
	return new Coords_CH(mc, [ Coords_CH.x_0.v ]);
    }

    gc_nv(mc, p, c_prev) {  // gc_nv = get Coords, new value
	return new Coords_CH(mc, p, c_prev, []);
    }

    get_max_num_t_steps() {
	return Trajectory.DEFAULT_MAX_NUM_T_STEPS
    }
}
