
///////////////////////////////////////////////////////////////////////////////////////////////
////////  IS = Ising Model (from SM = Statistical Mechanics)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

// NOTE: we represent the two Ising spin values as 0,1 "under the hood" -- it's more computationally convenient to think of as binary
// and, in creating debugging output, etc., 0 and 1 have the same width; only in ModelCalc_IS.get_E_spin_pair() where we switch to energy
// quantities do we have to translate 0, 1 the to the more physically appropriate -1, 1 via as1n1(); all other methods dealing with energies
// take output from get_E_spin_pair(), and all methods that take spin values as inputs always expect 0, 1
class ModelCalc_IS extends ModelCalc_Spin {

    constructor(N) {
	super(N);
	this.J = 1.0;  // not currently letting the interaction strength vary, but that may change...
    }

    as1n1(s01) {  // as1n1 = as 1, -1: translate from the more computationally convenient 0,1 to the more physically appropriate -1, 1
	return (2*s01 - 1.0);
    }

    get_E_spin_pair(sA, sB) {
	return -1.0 * this.J * this.as1n1(sA) * this.as1n1(sB);  // for the moment, J = 1 and h = 0
    }
}

class Params_IS extends Params_Spin {

    static T = undefined;  // = new UINI_float(this, "UI_P_SM_IS_T", true);  assignment occurs in UserInterface(); see discussion there

    push_vals_to_UI() {
	Params_IS.T.push_to_UI(this.T);
    }
}

class Coords_IS extends Coords_Spin {

    static N = undefined;  // = new UINI_int(this, "UI_P_SM_IS_N", false);  assignment occurs in UserInterface(); see discussion there    

    constructor(...args) {  // see discussion of # args at definition of abstract Coords()

	super(...args);

	// TRY TO MOVE CHUNKS OF THIS UP TO Coords_Spin()?

	// this.spins, the data structure the holds the square grid of spins
	// * we represent spin down as 0 and spin up as 1
	// * dtype binary would have been perfect, but it's not implemented (yet?)
	// * also, the underlying array doesn't appear to be available at this.spins.data as it was with int32... stride issue?
	this.prev_trans;  // will "point" to previous CoordTransition_Spin object for backward navigation
	this.next_trans;  // will "point" to next CoordTransition_Spin object for forward navigation

	let N = this.extra_args[0];

	if (this.constructing_init_cond) {

	    this.spins = zeros([ N, N ], {'dtype': 'int8'});
	    Coords_Spin.randomize_spins_arr(this.spins, this.get_rand_spin_val.bind(this));
	    this.prev_trans = null;  // since IC was not arrived at via a transition

	} else {

	    let x = this.get_rand_index(N);  // pick a random site, both x...
	    let y = this.get_rand_index(N);  // ... and y
	    let old_val = this.c_prev.spins.get(x, y);
	    let proposed_val = this.flip_spin_val(old_val);  // only possible basic move is to flip the spin
	    let Delta_E = this.mc.get_Delta_E_proposed_move(proposed_val, this.c_prev.spins, x, y);
	    let move_accepted = this.mc.accept_move(Delta_E, this.p.T);
	    this.spins = copy(this.c_prev.spins);
	    if (move_accepted) {
		this.flip_spin(this.spins, x, y);
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
	return this.mc.discunif_rng(0, 1);  // [0, 1]
    }

    flip_spin_val(curr_val) {  // alternatively, we could return ((curr_val + 1) % 2)... safer?
	return 1 - curr_val;  // 1 --> 1 - 1 = 0  and  0 --> 1 - 0 = 1
    }

    flip_spin(sa, x, y) {  // sa = spins array
	sa.set(x, y, this.flip_spin_val(sa.get(x, y)));
    }

    flip_random_spin() {
	let x = get_rand_index(this.mc.N);
	let y = get_rand_index(this.mc.N);
	this.flip_spin(this.spins, x, y);
    }
}

// * note that, in this model, N is # rows/columns in square array of spins, whereas, in the 1D SP models, N + 1 sites ran [0, N]
class Trajectory_IS extends Trajectory_Stoch {

    constructor(sim) {

	super(sim);
    }

    gmc() {  // gmc = get ModelCalc object
	return new ModelCalc_IS(Coords_IS.N.v);
    }

    gp() {  // gp = get Params object
	return new Params_IS(Params_IS.T.v);
    }

    gc_ic(mc) {  // gc_ic = get Coords, initial condition
	return new Coords_IS(mc, [ Coords_IS.N.v ]);
    }

    gc_nv(mc, p, c_prev) {  // gc_nv = get Coords, new value
	return new Coords_IS(mc, p, c_prev, [ Coords_IS.N.v ]);
    }

    set_rng_states_from_edge_vals() {
	this.mc.unif01_rng.state = this.get_x(this.t_edge).unif01_rng_state;
	this.mc.discunif_rng.state = this.get_x(this.t_edge).discunif_rng_state;
    }	

    get_max_num_t_steps() {
	return Trajectory.DEFAULT_MAX_NUM_T_STEPS
    }
}
