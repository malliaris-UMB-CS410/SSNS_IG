
//
// PlottingMachinery -- class that encapsulates various things related to plot output
//
// * single instance exists within UserInterface
// * contains, e.g., functionality for switching between different ways of displaying trajectory data
//
class PlottingMachinery {

    constructor(sim, bsvs) {

	this.sim = sim;

	// js object of objects containing all PlotType instances
	// * will be indexed as plots[sim.ST][sim.PT], e.g., to call plot() update
	// * used in ui.update_plot_type_buttons() to only show exiting plot types for given system type as non-disabled
	// * once a plot type for given system type is used by user, that "setting" is remembered in active_PT_by_ST js object
	// * listing order of plot types for given system type does matter in that the first listed will be active by default
	this.set_plot_container_dimensions(bsvs);
	this.plots = {};
	this.instantiate_all_PlotTypes();
	this.active_PT_by_ST = {};  // stores the active PT for each ST so it can be "remembered"
	for (let ST_val of Simulator.registered_STs.values()) {
	    this.active_PT_by_ST[ST_val] = Object.keys(this.plots[ST_val])[0];  // default value is the one listed first in pm.plots
	}
    }

    // listing order of plot types for given system type does matter in that the first listed will be active by default
    instantiate_PlotTypes_ST(ST_str) {

	switch(ST_str) {

	case "RW":
	    this.plots["RW"] = {};
	    this.plots["RW"]["HX"] = new PlotTypeHX_RW(this.sim.trjs["RW"]);
	    break;
	case "MN":
	    this.plots["MN"] = {};
	    this.plots["MN"]["XT"] = new PlotTypeXT_MN(this.sim.trjs["MN"]);
	    this.plots["MN"]["HX"] = new PlotTypeHX_MN(this.sim.trjs["MN"]);
	    break;
	case "CH":
	    this.plots["CH"] = {};
	    this.plots["CH"]["HX"] = new PlotTypeHX_CH(this.sim.trjs["CH"]);
	    break;
	case "IS":
	    this.plots["IS"] = {};
	    this.plots["IS"]["HM"] = new PlotTypeHM_IS(this.sim.trjs["IS"]);
	    break;
	case "XY":
	    this.plots["XY"] = {};
	    this.plots["XY"]["HM"] = new PlotTypeHM_XY(this.sim.trjs["XY"]);
	    break;
	case "LM":
	    this.plots["LM"] = {};
	    this.plots["LM"]["XT"] = new PlotTypeXT_LM(this.sim.trjs["LM"]);
	    this.plots["LM"]["PP"] = new PlotTypePP_LM(this.sim.trjs["LM"]);
	    break;
	case "GM":
	    this.plots["GM"] = {};
	    this.plots["GM"]["PP"] = new PlotTypePP_GM(this.sim.trjs["GM"]);
	    break;
	case "PF":
	    this.plots["PF"] = {};
	    this.plots["PF"]["XT"] = new PlotTypeXT_PF(this.sim.trjs["LM"]);  // TEMPORARY, TO ALLOW FOR TESTING!!!
	    break;  // DUMMY FOR TESTING.... FIX WHEN READY
	default:
	    throw new Error("ERROR 298290: Invalid ST_str!  Exiting...");
	}
    }

    instantiate_all_PlotTypes() {

	for (let ST_val of Simulator.registered_STs.values()) {
	    this.instantiate_PlotTypes_ST(ST_val);
	}
    }

    set_plot_container_dimensions(bsvs) {  // ALL LENGTHS ARE IN px UNLESS STATED OTHERWISE

	// gather info and calculate basic quantities
	let marg_L = parseInt($(".plot_target").css("margin-left"));
	let marg_R = parseInt($("#PLOT_CONTAINER_margin_right_px").val());
	let viewport_width = $(window).width();
	let max_usable_width = viewport_width - marg_L - marg_R;
	let viewport_height = $(window).height();
	let vh_1_2 = roundn(0.42 * viewport_height, 0);  // 1/2 viewport height, rounded to integer # px
	let vh_5_8 = roundn(0.625 * viewport_height, 0);  // 5/8 viewport height, rounded to integer # px
	let screen_is_XL = (bsvs == 5);  //  bsvs = Bootstrap viewport size            (viewport_width > 1200);  // 1200 px is Bootstrap4's breakpoint between "large" and "extra large" screen widths...

	// set 4 values based on viewport width; allowing future possibility of 8 distinct values, so leave as-is for straightforward editing
	let spw = (screen_is_XL ? vh_5_8 : max_usable_width);  // spw = square plot width
	let sph = (screen_is_XL ? vh_5_8 : max_usable_width);  // sph = square plot height
	let nspw = (screen_is_XL ? max_usable_width : max_usable_width);  // nspw = non-square plot width
	let nsph = (screen_is_XL ? vh_1_2 : vh_1_2);  // nsph = non-square plot height

	// store those 4 values
	PlotType.set_static_member_vals(spw, sph, nspw, nsph);  // used for on-the-fly adjustment of axis label element widths, etc.
    }
}
