
///////////////////////////////////////////////////////////////////////////////////////////////
////////  MN = Moran Model (from SP = Stochastic Processes, population genetics, etc.)  ///////
///////////////////////////////////////////////////////////////////////////////////////////////

class ModelCalc_MN extends ModelCalc_SP {

    constructor(N) {
	super();
	this.N = N;
    }

    get_P_bear(p, x) {  // probability of choosing and individual with allele X to potentially bear offspring
	return ( x - p.s*x ) / ( this.N - p.s*x );  // this sampling is weighted by fitness
    }

    get_P_die(p, x) {  // probability of choosing and individual with allele X to depart the population, i.e., die
	return x / this.N;  // this sampling is unweighted
    }

    get_P_step_R(p, x) {  // probability (conditioned on value of x) of x increasing
	let P_bear = this.get_P_bear(p, x);
	let P_die = this.get_P_die(p, x);
	return (1.0 - P_die) * ( P_bear*(1.0 - p.mu) + (1.0 - P_bear)*p.mu );  // natural boundary condition built-in, i.e., P_step_R(x = N) = 0
    }

    get_P_step_L(p, x) {  // probability (conditioned on value of x) of x decreasing
	let P_bear = this.get_P_bear(p, x);
	let P_die = this.get_P_die(p, x);
	return P_die * ( P_bear*p.mu + (1.0 - P_bear)*(1.0 - p.mu) );  // natural boundary condition built-in, i.e., P_step_L(x = 0) = 0
    }
}

class Params_MN extends Params {

    // NOTE: N in 2-allele Moran model serves a dual role: it's not only the size of the domain, but also a parameter (as 1/N, the strength of
    //       genetic drift); either way, it is fixed for the Trajectory, and thus not included here (see Trajectory_SP constructor)

    static mu = undefined;  // = new UINI_float(this, "UI_P_SP_MN_mu", true);  assignment occurs in UserInterface(); see discussion there
    static s = undefined;  // = new UINI_float(this, "UI_P_SP_MN_s", true);  assignment occurs in UserInterface(); see discussion there

    constructor(mu_val, s_val) {
	super();
	this.mu = mu_val;
	this.s = s_val;
    }

    push_vals_to_UI() {
	Params_MN.mu.push_to_UI(this.mu);
	Params_MN.s.push_to_UI(this.s);
    }

    get_info_str() {
	return "mu, s = " + this.mu + ", " + this.s;
    }
}

class Coords_MN extends Coords_SP_finite {

    static N = undefined;  // = new UINI_int(this, "UI_P_SP_MN_N", false);  assignment occurs in UserInterface(); see discussion there
    static x_0 = undefined;  // = new UINI_int(this, "UI_P_SP_MN_x_0", false);  assignment occurs in UserInterface(); see discussion there

    constructor(...args) {  // "..." is Javascript spread operator
	super(...args);
    }
}

class Trajectory_MN extends Trajectory_SP {

    // it's not helpful to pass this constructor Coords_MN.N.v as an arg, store as this.N, and use as arg for other methods; see notes at Trajectory()
    constructor(sim) {
	super(sim);
    }

    gmc() {  // gmc = get ModelCalc object
	return new ModelCalc_MN(Coords_MN.N.v);
    }

    gp() {  // gp = get Params object
	return new Params_MN(Params_MN.mu.v, Params_MN.s.v);
    }

    gc_ic(mc) {  // gc_ic = get Coords, initial condition
	return new Coords_MN(mc, [ Coords_MN.N.v, Coords_MN.x_0.v ]);
    }

    gc_nv(mc, p, c_prev) {  // gc_nv = get Coords, new value
	return new Coords_MN(mc, p, c_prev, [ Coords_MN.N.v ]);
    }

    get_max_num_t_steps() {
	return Trajectory.DEFAULT_MAX_NUM_T_STEPS
    }
}
