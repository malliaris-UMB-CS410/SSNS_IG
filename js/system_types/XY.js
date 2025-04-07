
///////////////////////////////////////////////////////////////////////////////////////////////
////////  XY = XY Model (from SM = Statistical Mechanics)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

// NOTE: under the hood, we represent each XY spin value as a float on [0, 1), i.e., fraction of 1 revolution
class ModelCalc_XY extends ModelCalc_Spin {

    constructor(N) {
	super(N);
	this.J = 1.0;  // not currently letting the interaction strength vary, but that may change...
    }

    get_E_spin_pair(sA, sB) {
	return -1.0 * this.J * Math.cos( 2.0 * Math.PI * (sB - sA));  // for the moment, J = 1 and h = 0
    }
}

class Params_XY extends Params_Spin {

    static T = undefined;  // = new UINI_float(this, "UI_P_SM_XY_T", true);  assignment occurs in UserInterface(); see discussion there

    push_vals_to_UI() {
	Params_XY.T.push_to_UI(this.T);
    }
}

class Coords_XY extends Coords_Spin {

    static N = undefined;  // = new UINI_int(this, "UI_P_SM_XY_N", false);  assignment occurs in UserInterface(); see discussion there    

    constructor(...args) {  // see discussion of # args at definition of abstract Coords()

	super(...args);

	// TRY TO MOVE CHUNKS OF THIS UP TO Coords_Spin()?
	
	// this.spins, the data structure the holds the square grid of spins
	// * the spin is of unit length and rotates in the x-y plane; we represent it as a number on [0, 2*pi)
	this.prev_trans;  // will "point" to previous CoordTransition_Spin object for backward navigation
	this.next_trans;  // will "point" to next CoordTransition_Spin object for forward navigation

	let N = this.extra_args[0];

	if (this.constructing_init_cond) {

	    this.spins = zeros([ N, N ], {'dtype': 'float64'});
	    Coords_Spin.randomize_spins_arr(this.spins, this.get_rand_spin_val.bind(this));
	    this.prev_trans = null;  // since IC was not arrived at via a transition

	} else {

	    let x = this.get_rand_index(N);  // pick a random site, both x...
	    let y = this.get_rand_index(N);  // ... and y
	    let old_val = this.c_prev.spins.get(x, y);
	    let proposed_val = this.get_perturbed_spin_val(old_val);  // give the current spin a random kick to create proposed spin
	    let Delta_E = this.mc.get_Delta_E_proposed_move(proposed_val, this.c_prev.spins, x, y);
	    let move_accepted = this.mc.accept_move(Delta_E, this.p.T);
	    this.spins = copy(this.c_prev.spins);
	    if (move_accepted) {
		this.spins.set(x, y, proposed_val);
	    }
	    let new_val = this.spins.get(x, y);
	    let pts = new CoordTransition_Spin(x, y, old_val, new_val, move_accepted, Delta_E);
	    this.c_prev.next_trans = pts;
	    this.prev_trans = pts;
	}

	// regardless of how the data was calculated, we record the state of the rng(s) for potential future use
	this.unif01_rng_state = this.mc.unif01_rng.state;
	this.discunif_rng_state = this.mc.discunif_rng.state;
    }

    get_rand_spin_val() {
	return this.mc.unif01_rng();  // [0, 1), i.e., uniform probability over all angles (measured in fractions of 1 revolution)
    }

    get_perturbed_spin_val(s) {

        let perturbed_s = ((this.mc.unif01_rng() - 0.5) / 1.0);  // 1/4 is just heuristic that speeds up approach; consider refining...
        if (perturbed_s < 0) perturbed_s += 1.0;  // in case random kick sent theta outside the range [0.0, 1.0)
        else if (perturbed_s >= 1.0) perturbed_s -= 1.0;  // in case random kick sent theta outside the range [0.0, 1.0)
	return perturbed_s;
    }
}

// * note that, in this model, N is # rows/columns in square array of spins, whereas, in the 1D SP models, N + 1 sites ran [0, N]
class Trajectory_XY extends Trajectory_Stoch {

    constructor(sim) {

	super(sim);
    }

    gmc() {  // gmc = get ModelCalc object
	return new ModelCalc_XY(Coords_XY.N.v);
    }

    gp() {  // gp = get Params object
	return new Params_XY(Params_XY.T.v);
    }

    gc_ic(mc) {  // gc_ic = get Coords, initial condition
	return new Coords_XY(mc, [ Coords_XY.N.v ]);
    }

    gc_nv(mc, p, c_prev) {  // gc_nv = get Coords, new value
	return new Coords_XY(mc, p, c_prev, [ Coords_XY.N.v ]);
    }

    set_rng_states_from_edge_vals() {
	this.mc.unif01_rng.state = this.get_x(this.t_edge).unif01_rng_state;
	this.mc.discunif_rng.state = this.get_x(this.t_edge).discunif_rng_state;
    }	

    get_max_num_t_steps() {
	return Trajectory.DEFAULT_MAX_NUM_T_STEPS
    }
}
