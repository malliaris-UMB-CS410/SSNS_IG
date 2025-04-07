
///////////////////////////////////////////////////////////////////////////////////////////////
////////  GM = Gingerbread-man Map (from ND = Nonlinear Dynamics)  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

class ModelCalc_GM extends ModelCalc {

    constructor() {
	super();
    }

    model_is_stoch() {return false; }

    get_x_new(x, y) {
	return (1.0 - y + Math.abs(x));
    }

    get_y_new(x, y) {
	return x;
    }
}

class Params_GM extends Params {  // NOTE: this is a 0-parameter model!!!  x_0 and y_0 are "parameters" in interface, but kept in Coords_GM

    constructor() {
	super();
    }
    push_vals_to_UI() {}  // needed only because Params requires it
    get_info_str() {  // needed only because Params requires it
	return "** 0-PARAMETER MODEL **";
    }
}

class Coords_GM extends Coords {

    static x_0 = undefined;  // = new UINI_float(this, "UI_P_ND_GM_x_0", false);  assignment occurs in UserInterface(); see discussion there
    static y_0 = undefined;  // = new UINI_float(this, "UI_P_ND_GM_y_0", false);  assignment occurs in UserInterface(); see discussion there

    constructor(...args) {  // see discussion of # args at definition of abstract Coords()

	super(...args);

	if (this.constructing_init_cond) {

	    this.x = this.extra_args[0];
	    this.y = this.extra_args[1];

	} else {

	    this.x = this.mc.get_x_new(this.c_prev.x, this.c_prev.y);
	    this.y = this.mc.get_y_new(this.c_prev.x, this.c_prev.y);
	}
    }

    output() {
	console.log("x =", this.x);
	console.log("y =", this.y);
    }
}

class Trajectory_GM extends Trajectory {

    constructor(sim) {
	super(sim);
    }

    gmc() {  // gmc = get ModelCalc object
	return new ModelCalc_GM();
    }

    gp() {  // gp = get Params object
	return new Params_GM();
    }

    gc_ic(mc) {  // gc_ic = get Coords, initial condition
	return new Coords_GM(mc, [Coords_GM.x_0.v, Coords_GM.y_0.v]);  // packaging initial coords helps distinguish 2-argument constructor...
    }

    gc_nv(mc, p, c_prev) {  // gc_nv = get Coords, new value
	return new Coords_GM(mc, p, c_prev, []);  // ... from 3-argument constructor
    }

    get_max_num_t_steps() {
	return Trajectory.DEFAULT_MAX_NUM_T_STEPS
    }
}
