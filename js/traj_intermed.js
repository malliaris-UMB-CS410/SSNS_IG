
///////////////////////////////////////////////////////////////////////////////////////////////
////////  SSNS          ///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////////  Javascript code for Simple Stochastic and Nonlinear Simulator (SSNS)  ///////////////
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
////////  INTERMEDIATE CLASSES  ///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////////  each of these facilitates code reuse among a subset of the models  //////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

class ModelCalc_Stoch extends ModelCalc {

    static rng_seed = undefined;  // = new UINI_int(this, "UI_CTRL_rng_seed", false);  assignment occurs in UserInterface(); see discussion there

    constructor() {

	super();

	this.unif01_rng = randu.factory({'seed': ModelCalc_Stoch.rng_seed.v });
	this.discunif_rng = discreteUniform.factory({'seed': ModelCalc_Stoch.rng_seed.v });
	this.binom_rng = binomial.factory({'seed': ModelCalc_Stoch.rng_seed.v });
	console.log("INFO:\tusing PRNG algorithm Mersenne Twister 19937 (the default) on all:", this.unif01_rng.PRNG.NAME, this.discunif_rng.PRNG.NAME, this.binom_rng.PRNG.NAME);
	console.log("INFO:\tusing seed value = " + ModelCalc_Stoch.rng_seed.v);
    }

    model_is_stoch() {return true; }

    // wrapper for stdlib.js's binomial sampling that correctly returns k = 0 for n = 0 draws rather than stdlib.js's NaN
    binom_rng_0(n, p) {
	if (n === 0) {
	    return 0;
	} else {
	    return this.binom_rng(n, p);
	}
    }

    test_rngs() {

	console.log("TESTING RNGs...");
	let TESTunif01_rng = randu.factory({'seed': ModelCalc_Stoch.rng_seed.v });
	let num_below = 0;
	for (let i = 0; i < 10000; i++) {
	    let rn01 = TESTunif01_rng();
	    //console.log(i, rn01);
	    //tic()[0]
	    num_below += (rn01 < 0.5) ? 1 : 0;
	}
	console.log("num_below = ", num_below);
	console.log("TESTS COMPLETED.");
    }
}

class Trajectory_Stoch extends Trajectory {

    constructor(sim) {

	super(sim);
	if (!this.set_rng_states_from_edge_vals) throw new Error("Derived Trajectory_Stoch must define set_rng_states_from_edge_vals()");
    }
}

// * CoordTransition_Spin is a specialty class -- it doesn't extend anything
// * it is used by spin "_Spin" systems _IS and _XY to link consecutive Coords objects via "next_trans", "prev_trans" pointers
// * may be useful to avoid storing entire NxN array of spins in each Coords if N is large and memory usage becomes an issue
// * each "Transition" is iteration of Metropolis-Hastings algorithm in which "move" may or may not occur, i.e., system state may be unchanged!
class CoordTransition_Spin {

    constructor(x, y, old_val, new_val, move_occurred, Delta_E) {

	this.x = x;  // x index of spin being considered
	this.y = y;  // y index of spin being considered
	this.old_val = old_val;  // value of spin being considered **before** transition, i.e., at smaller t val
	this.new_val = new_val;  // value of spin being considered **after** transition, i.e., at larger t val
	this.move_occurred = move_occurred;  // whether proposed move was accepted
	this.Delta_E = Delta_E;  // change in energy between proposed and current states -- useful for compiling statistics
    }
}

class ModelCalc_SP extends ModelCalc_Stoch {

    constructor() {

	super();

	if (!this.get_P_step_R) throw new Error("Derived ModelCalc_SP must define get_P_step_R()");
	if (!this.get_P_step_L) throw new Error("Derived ModelCalc_SP must define get_P_step_L()");
    }

    get_x_new(p, x) {
	let P_step_R = this.get_P_step_R(p, x);
	let P_step_L = this.get_P_step_L(p, x);
	let rn01 = this.unif01_rng();
	if (rn01 < P_step_R) return (x + 1);
	else if ( rn01 < (P_step_R + P_step_L) ) return (x - 1);
	else return x;
    }

    //    //arr.set(0, 0, 0.5*this.num_x_vals * this.num_ens_members_per_x);  // delta function
    //	TrajSeg_SP  get_x_arr() { return this.cv.map((item, index, array) => item.x);
}

class Coords_SP extends Coords {

    static num_IEM = undefined;  // = new UINI_int(this, "UI_CTRL_SP_NI", false);  assignment occurs in UserInterface(); see discussion there
    static num_GEM = undefined;  // = new UINI_int(this, "UI_CTRL_SP_NG", false);  assignment occurs in UserInterface(); see discussion there

    constructor(...args) {  // see discussion of # args at definition of abstract Coords()

	super(...args);
    }
}

class Coords_SP_finite extends Coords_SP {

    constructor(...args) {  // see discussion of # args at definition of abstract Coords()

	super(...args);

	let N = this.extra_args[0];
	let num_x = N + 1;  // # x sites = N + 1

	// data structure allocation
	if (Coords_SP.num_IEM.v > 0) {
	    this.x_indiv = zeros([ Coords_SP.num_IEM.v ], {'dtype': 'int32'});
	}
	if (Coords_SP.num_GEM.v > 0) {
	    this.H_x_group = zeros([ num_x ], {'dtype': 'int32'});
	}

	if (this.constructing_init_cond) {

	    let x_0 = this.extra_args[1];
	    if (x_0 > N) {
		throw new Error("ERROR 970132: invalid x_0 value in Coords_SP().  Exiting...");
	    }
	    if (Coords_SP.num_IEM.v > 0) {
		for (let i = 0; i < Coords_SP.num_IEM.v; i++) {
		    this.x_indiv.set(i, x_0);
		}
	    }
	    if (Coords_SP.num_GEM.v > 0) {
		this.H_x_group.set(x_0, Coords_SP.num_GEM.v);
	    }

	} else {

	    if (Coords_SP.num_IEM.v > 0) {
		for (let i = 0; i < Coords_SP.num_IEM.v; i++) {
		    let x_prev = this.c_prev.x_indiv.get(i);
		    this.x_indiv.set(i, this.mc.get_x_new(this.p, x_prev));
		}
	    }
	    if (Coords_SP.num_GEM.v > 0) {

		for (let x = 0; x < num_x; x++) {
		    
		    let num_at_x = this.c_prev.H_x_group.get(x);
		    let P_step_R = this.mc.get_P_step_R(this.p, x);
		    let P_step_L = this.mc.get_P_step_L(this.p, x);
		    let P_step_L_given_non_step_R = P_step_L / (1.0 - P_step_R);  // need this conditional probability since we're doing 2 binomial draws
		    let num_step_R = this.mc.binom_rng_0(num_at_x, P_step_R);  // one binomial draw to determine # moving R...
		    let num_stay_put_or_step_L = num_at_x - num_step_R;
		    let num_step_L = this.mc.binom_rng_0(num_stay_put_or_step_L, P_step_L_given_non_step_R);  // another to determine # moving L...
		    let num_stay_put = num_at_x - num_step_R - num_step_L;
		    this.H_x_group.data[x] += num_stay_put;
		    this.H_x_group.data[Math.max(0, x - 1)] += num_step_L;
		    this.H_x_group.data[Math.min(N, x + 1)] += num_step_R;
		}
	    }
	}

	// regardless of how the data was calculated, we record the state of the rng(s) for potential future use
	this.unif01_rng_state = this.mc.unif01_rng.state;
	this.binom_rng_state = this.mc.binom_rng.state;
    }
}

class Coords_SP_semiinf extends Coords_SP {

    constructor(...args) {  // see discussion of # args at definition of abstract Coords()

	super(...args);

	// data structure allocation
	if (Coords_SP.num_GEM.v > 0) {
	    this.H_x_group = new OrderedMap();
	}

	if (this.constructing_init_cond) {

	    let x_0 = this.extra_args[0];
	    if (Coords_SP.num_GEM.v > 0) {
		this.H_x_group.setElement(x_0, Coords_SP.num_GEM.v);
	    }

	} else {

	    if (Coords_SP.num_GEM.v > 0) {

		this.c_prev.H_x_group.forEach((element, index) => {

		    let x = element[0];
		    let num_at_x = element[1];
		    let P_step_R = this.mc.get_P_step_R(this.p, x);
		    let P_step_L = this.mc.get_P_step_L(this.p, x);
		    let P_step_L_given_non_step_R = P_step_L / (1.0 - P_step_R);  // need this conditional probability since we're doing 2 binomial draws
		    let num_step_R = this.mc.binom_rng_0(num_at_x, P_step_R);  // one binomial draw to determine # moving R...
		    let num_stay_put_or_step_L = num_at_x - num_step_R;
		    let num_step_L = this.mc.binom_rng_0(num_stay_put_or_step_L, P_step_L_given_non_step_R);  // another to determine # moving L...
		    let num_stay_put = num_at_x - num_step_R - num_step_L;
		    CU.add_to_entry_OM(this.H_x_group, x, num_stay_put);
		    CU.add_to_entry_OM(this.H_x_group, Math.max(0, x - 1), num_step_L);  // reflecting boundary at x = 0
		    CU.add_to_entry_OM(this.H_x_group, x + 1, num_step_R);  // no boundary (more or less) on other side
		});
		//CU.ppOM(this.H_x_group);///////
		//CU.print_stats_OM(this.H_x_group);/////////
	    }
	}

	// regardless of how the data was calculated, we record the state of the rng(s) for potential future use
	this.unif01_rng_state = this.mc.unif01_rng.state;
	this.binom_rng_state = this.mc.binom_rng.state;
    }
}

class Trajectory_SP extends Trajectory_Stoch {

    constructor(sim) {

	super(sim);

	console.log("INFO:\tconstructing Trajectory_SP with x on [0, N] or [0, +inf]");
	console.log(Coords_SP.num_IEM.v, ",", Coords_SP.num_GEM.v, " individual/group ensemble members");
    }

    get_max_num_t_steps() {
	return Trajectory.DEFAULT_MAX_NUM_T_STEPS
    }

    set_rng_states_from_edge_vals() {
	this.mc.unif01_rng.state = this.get_x(this.t_edge).unif01_rng_state;
	this.mc.binom_rng.state = this.get_x(this.t_edge).binom_rng_state;
    }	
}

// * 2D only
// * takes into account nearest neighbor interactions only
class ModelCalc_Spin extends ModelCalc_Stoch {

    constructor(N) {
	super();
	this.N = N;
    }

    get_Boltzmann_factor(Delta_E, T) {  // used via Metropolis-Hastings algorithm to build the canonical ensemble of SM
	return Math.exp(-1.0 * Delta_E / T);
    }

    gilw(z, N) {  // gilw = get index lower and wrap (periodic boundary conditions); works for both x and y
	return ((z + N - 1) % N);
    }

    gihw(z, N) {  // gihw = get index higher and wrap (periodic boundary conditions); works for both x and y
	return ((z + 1) % N);
    }

    // energy of some spin s (not necessarily equal to spin at spins(x, y)) if placed at x,y in spins array
    // (this separation is handy for evaluating the energy of a proposed new spin value at x, y, calculating Delta_E, etc.)
    get_E_spin_in_neighborhood(s, spins, x, y) {

	let N = this.N;  // for convenience
	let E = 0.0;
	E += this.get_E_spin_pair(s, spins.get(this.gilw(x, N), y));  // interaction with spin to the L
	E += this.get_E_spin_pair(s, spins.get(this.gihw(x, N), y));  // interaction with spin to the R
	E += this.get_E_spin_pair(s, spins.get(x, this.gilw(y, N)));  // interaction with spin below
	E += this.get_E_spin_pair(s, spins.get(x, this.gihw(y, N)));  // interaction with spin above
	return E;
    }

    get_E_at_site(spins, x, y) {  // energy of spin at x, y due to rest of spins in array

	return this.get_E_spin_in_neighborhood(spins.get(x, y), spins, x, y);
    }

    get_Delta_E_proposed_move(proposed_s, spins, x, y) {

	let E_f = this.get_E_spin_in_neighborhood(proposed_s, spins, x, y);  // final energy were proposed move to be made
	let E_i = this.get_E_at_site(spins, x, y);  // initial (i.e., current) energy
	return (E_f - E_i);
    }

    accept_move(Delta_E, T) {  // (rather than reject); used to carry out Metropolis-Hastings algorithm

	if (Delta_E <= 0) {  // "improvements" (and "lateral" steps) are always accepted

	    return true;

	} else {  // "setbacks" (i.e., steps up in energy) are accepted with T-dependent probability

	    let P_accept = this.get_Boltzmann_factor(Delta_E, T)
	    return (this.unif01_rng() < P_accept);
	}
    }
}

class Params_Spin extends Params {

    constructor(T_val) {
	super();
	this.T = T_val;  // absolute temperature (measured in units of energy)
    }

    get_info_str() {
	return "T = " + this.T;
    }
}

class Coords_Spin extends Coords {

    static randomize_spins_arr(sa, get_rand_spin_val_fxn) {

	let N = sa.shape[0];  // either 0 or 1 works, since array is N x N square
	for (let x = 0; x < N; x++) {
	    for (let y = 0; y < N; y++) {
		sa.set(x, y, get_rand_spin_val_fxn());
	    }
	}
    }

    static output_spins_arr(sa) {

	let N = sa.shape[0];  // either 0 or 1 works, since array is N x N square
	let output_str = ""
	for (let x = 0; x < N; x++) {
	    for (let y = 0; y < N; y++) {
		output_str += sa.get(x, y) + " ";
	    }
	    output_str += "\n";
	}
	console.log(output_str);
    }

    constructor(...args) {  // see discussion of # args at definition of abstract Coords()
	super(...args);
    }

    get_rand_index(N) {
	return this.mc.discunif_rng(0, N - 1);  // [0, N - 1]
    }

    output() {
	if (this.prev_trans != null) {
	    console.log("prev_trans x,y,mo,DE =", this.prev_trans.x, this.prev_trans.y, this.prev_trans.move_occurred, this.prev_trans.Delta_E);
	}
	if (this.next_trans != null) {
	    console.log("next_trans x,y,mo,DE =", this.next_trans.x, this.next_trans.y, this.next_trans.move_occurred, this.next_trans.Delta_E);
	}
    }
}
