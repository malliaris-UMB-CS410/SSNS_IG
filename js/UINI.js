
// UINI: User Interface Numeric Input
// helper class for connecting UI controls to internal parameters (in Params classes, etc.)
class UINI {

    static yields_valid_int(str_val) {

	if (str_val == "") {
	    console.log("INFO:\tyields_valid_int():\tignoring non-numeric input...");
	    return false;
	}
	let valid_input_pttn = /^0$|^-?[1-9][0-9]*$/;  // reject leading zeros, scientific notation, etc.
	if ( ! valid_input_pttn.test(str_val) ) {
	    console.log("INFO:\tyields_valid_int():\tignoring invalid numeric input (leading zeros, scientific notation, etc.)...");
	    return false;
	}
	if ( ! Number.isInteger(Number(str_val)) ) {
	    console.log("INFO:\tyields_valid_int():\tignoring valid numerical input which is not an integer...");
	    return false;
	}
	return true;
    }
    
    static yields_valid_float(str_val) {

	if (str_val == "") {
	    console.log("INFO:\tyields_valid_float():\tignoring non-numeric input...");
	    return false;
	}
	let valid_input_pttn = /^-?[0-9.+-e]+$/;  // consider making this more stringent/specific...
	if ( ! valid_input_pttn.test(str_val) ) {
	    console.log("INFO:\tyields_valid_float():\tignoring invalid numeric input...");
	    return false;
	}
	return true;
    }

    constructor(ui, id_str, indicate_params_changed) {

	if (!this.yields_valid_value) throw new Error("Derived UINI must define yields_valid_value()");
	if (!this.parse_from_str) throw new Error("Derived UINI must define parse_from_str()");

	this.ui = ui;
	this.id_str = id_str;
	this.indicate_params_changed = indicate_params_changed;
	this.v = undefined;  // holds the "official" (validated) value of the quantity

	// initial loading and validation of value
	let raw_default_str_val = CU.gs(this.id_str);
	if (this.yields_valid_value(raw_default_str_val)) {
	    this.set_min_max_or_val(this.parse_from_str(raw_default_str_val));
	} else {
	    throw new Error("UINI() failed to parse a valid default value for input " + this.id_str + ".  Exiting...");
	}
    }

    sv(val) {  // sv = set value method sets both internal variable and interface variable, thus keeping the two in sync
	this.v = val;
	CU.sv(this.id_str, val);
    }

    set_min_max_or_val(proposed_new_val) {

	let min = CU.gmn(this.id_str);  // for convenience
	let max = CU.gmx(this.id_str);  // for convenience
	if ((min == undefined) || (max == undefined)) {
	    throw new Error("Either min or max for " + this.id_str + " is undefined, so UINI can't check range.  Exiting...");
	}
	if (proposed_new_val < min) {
	    console.log("INFO:\t" + this.id_str + " UINI parsed value of " + proposed_new_val
			+ " which is out-of-range LOW... setting to min val of " + min + " instead...");
	    this.sv(min);  // we parsed a valid number, but it is out-of-range low, so we set to min, since it's closest valid value
	} else if (proposed_new_val > max) {
	    console.log("INFO:\t" + this.id_str + " UINI parsed value of " + proposed_new_val
			+ " which is out-of-range HIGH... setting to max val of " + max + " instead...");
	    this.sv(max);  // we parsed a valid number, but it is out-of-range high, so we set to max, since it's closest valid value
	} else {
	    this.sv(proposed_new_val);  // we parsed a valid number, and it's in range, so set it without correction
	}
    }
    
    try_update(str_val) {  // called, e.g., by "onChange()" in HTML input element; will update value if all conditions are met

	if (this.yields_valid_value(str_val)) {
	    this.set_min_max_or_val(this.parse_from_str(str_val));
	    if (this.indicate_params_changed) {
		this.ui.indicate_new_param_vals_ready_to_pull_UI_to_traj();
	    }
	} else {
	    console.log("INFO:\t" + this.id_str + " try_update() failed to parse a valid number... restoring previous valid value to UI...");
	    this.sv(this.v);  // if proposed new value is invalid, push current valid value back to UI
	}
    }

    push_to_UI(val) {  // take an already validated value (perhaps from a Params_*) and push to UI input
	this.sv(val);
    }
}

class UINI_int extends UINI {

    constructor(ui, id_str, indicate_params_changed) {

	super(ui, id_str, indicate_params_changed);
    }

    yields_valid_value(str_val) {
	return UINI.yields_valid_int(str_val);
    }

    parse_from_str(str_val) {
	return parseInt(str_val);
    }
}	

class UINI_float extends UINI {

    constructor(ui, id_str, indicate_params_changed) {

	super(ui, id_str, indicate_params_changed);
    }

    yields_valid_value(str_val) {
	return UINI.yields_valid_float(str_val);
    }

    parse_from_str(str_val) {
	return parseFloat(str_val);
    }
}	
