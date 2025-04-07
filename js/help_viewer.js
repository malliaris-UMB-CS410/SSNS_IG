
// helper class that represents a single node in the hierarchy of help pages
class HelpViewerNode {

    static breadcrumbs_pre = '<span style="white-space: nowrap; "><i class="fas fa-xs fa-angles-right px-2"></i>';
    static breadcrumbs_post = '</span><wbr>';

    constructor(id_str, child_arr) {

	this.id_str = id_str;
	this.child_arr = child_arr;
	this.breadcrumbs_html = "";  // will be filled during recursion
	this.is_param_node = (this.id_str.substr(0, 5) == "HV_P_");  // id_str for a node/view corresponding to a system parameter begins "HV_P_..."
	this.is_ST_node = (this.id_str.substr(0, 6) == "HV_ST_");  // id_str for a node/view corresponding to a system type begins "HV_ST_..."
	let header_txt = $("#" + this.id_str).attr('data-hvmh');  // grab stored text in custom data-hvmh attribute; hvmh = helper view modal header
	let page_is_ready = parseInt($("#" + this.id_str).attr('data-ready'));  // TEMPORARY; REMOVE WHEN HELP PAGE WRITING IS COMPLETE
	/////console.log("YyYyYy", id_str, (page_is_ready ? "" : "not"), "ready");///////////////////////
	if (this.is_param_node) {  // nodes that represent system parameters are represented as LaTeX'd variables
	    this.header_html = katex.renderToString(header_txt, {throwOnError: false});
	} else {  // plain text for all others
	    this.header_html = header_txt;
	}
	this.header_html_as_link = this.get_link_str(this.header_html);
	//this.header_html_as_curr_page = '<span class="hv_curr_page">' + this.header_html + '</span>';
	if (page_is_ready) {  // TEMPORARY; REMOVE WHEN HELP PAGE WRITING IS COMPLETE
	    this.header_html_as_curr_page = '<span class="hv_curr_page">' + this.header_html + '</span>';  // TEMPORARY; REMOVE WHEN HELP PAGE WRITING IS COMPLETE
	} else {  // TEMPORARY; REMOVE WHEN HELP PAGE WRITING IS COMPLETE
	    this.header_html_as_curr_page = '<span class="hv_curr_page_not_ready">' + this.header_html + '</span>';  // TEMPORARY; REMOVE WHEN HELP PAGE WRITING IS COMPLETE
	}  // TEMPORARY; REMOVE WHEN HELP PAGE WRITING IS COMPLETE
	if (this.is_ST_node) {  // if it's a node corresponding to a system type (ST), add ability to click icon that closes viewer and loads that ST
	    let ST_abbrev = this.id_str.substr(6, 7);  // e.g., "LM" for Logistic Map
	    this.header_html_as_curr_page += String.raw`<a class="hv_link" onclick="window.sim.ui.close_HV_load_ST('` + ST_abbrev + String.raw`'); ">`;
	    this.header_html_as_curr_page += '<i class="fas fa-fw fa-sm fa-up-right-from-square" data-fa-transform="right-3 down-2"></i></a>';
	}
    }

    get_link_str(str_to_wrap) {
	return String.raw`<a class="hv_link" onclick="window.sim.ui.hv.show_view('` + this.id_str + String.raw`'); ">` + str_to_wrap + String.raw`</a>`;
    }

    // all are links except the last one which is the current node, i.e., the one being viewed
    get_breadcrumbs_html_frm_arr(breadcrumbs_arr) {
	let arr_to_iter = breadcrumbs_arr.concat(this.header_html_as_curr_page);  // append html for current page
	arr_to_iter = arr_to_iter.slice(1);  // trim first element, which is root node Home (handled separately)
	let result_html = "";
	for (let node_html of arr_to_iter.values()) {
	    result_html += HelpViewerNode.breadcrumbs_pre + node_html + HelpViewerNode.breadcrumbs_post;
	}
	return result_html;
    }

    get_all_child_links() {
	let child_links_html = "";
	for (let cn of this.child_arr.values()) {  // cn = child node
	    child_links_html += cn.header_html_as_link + "<br/>\n";
	}
	return child_links_html;
    }

    recurse_and_fill_data(hvn_map, breadcrumbs_arr) {

	//console.log("visiting", this.id_str);////////////
	hvn_map[this.id_str] = this;  // add this node to the map of nodes that will be used by HelpViewer
	this.breadcrumbs_html = this.get_breadcrumbs_html_frm_arr(breadcrumbs_arr);
	breadcrumbs_arr.push(this.header_html_as_link);  // append link to this node to breadcrumbs array
	if (this.child_arr.length > 0) {
	    for (let cn of this.child_arr.values()) {  // cn = child node
		cn.recurse_and_fill_data(hvn_map, breadcrumbs_arr);
	    }
	}
	breadcrumbs_arr.pop();  // remove this node's link now that recursion involving this node is complete
    }
}




// helper class to coordinate storage and display of help pages
//
// * the (very rudimentary) back button utilizes a page history of length 1, i.e., just the previously visited page;
//   note that this means two consecutive presses, and you're already "bouncing" between the same two pages!
class HelpViewer {

    constructor(sim, bsvs) {

	// basic/miscellaneous settings
	this.sim = sim;
	this.initial_view = "HV_P_ND_LM_x_0";  // default setting
	this.curr_view = "";
	this.prev_view = "";
	this.show_on_load = false;  // whether to show HelpViewer on app loading
	$("#md_container").on("hidden.bs.modal", function ()   { this.deployed = false;  });
	$("#md_container").on("shown.bs.modal", function ()   { this.deployed = true;  });

	// set image widths, first for all full width "fw", then the single screenshot; ALL LENGTHS ARE IN px UNLESS STATED OTHERWISE
	let fw_max_width = parseInt($("#hv_fw_img_max_width_px").val());  // read in max value, which is about width of mobile viewport
	let suggested_fw_width_based_on_viewport = roundn(0.90 * $(window).width(), 0);  // % of viewport width, rounded to integer # px
	let fw_width = Math.min(fw_max_width, suggested_fw_width_based_on_viewport);  // fill screen on smaller devices only
	$(".hv_fw_img").css("width", fw_width);
	let hv_ssns_screenshot_min_width = parseInt($("#hv_ssns_screenshot_min_width_px").val());  // read in screenshot min width val and...
	let hv_ssns_screenshot_width = (bsvs <= 3) ? hv_ssns_screenshot_min_width : fw_width;  // use it on xs, sm, and md viewports only
	$("#hv_ssns_screenshot").css("width", hv_ssns_screenshot_width);

	// use literals to specify structure of HV page network
	let hv_root =

	    new HelpViewerNode("HV_HOME", [

		new HelpViewerNode("HV_USING", [

		    new HelpViewerNode("HV_PLOT_AREA", []),
		    new HelpViewerNode("HV_MAIN_PANEL", []),
		    new HelpViewerNode("HV_C_CONTINUOUS_CMDS", []),
		    new HelpViewerNode("HV_PRMS_DROPDOWN", []),
		    new HelpViewerNode("HV_CTRL_DROPDOWN", [

			new HelpViewerNode("HV_C_ONE_STEP_CMDS", []),
			new HelpViewerNode("HV_C_PLOT_TYPES", []),
			new HelpViewerNode("HV_C_T_DT_JUMP", []),
			new HelpViewerNode("HV_C_DELAY_WINDOW", []),
			new HelpViewerNode("HV_C_T_0_MAX_STOP", []),
			new HelpViewerNode("HV_C_RNG_RECR_TRAJ", []),
			new HelpViewerNode("HV_C_SP_ENSEMBLE", []),
			new HelpViewerNode("HV_C_RELOAD_CMDS", []),

		    ]),

		    new HelpViewerNode("HV_SYS_DROPDOWN", []),
		    new HelpViewerNode("HV_HELP_VIEWER", []),

		]),

		new HelpViewerNode("HV_COMPUTATION", [

		    new HelpViewerNode("HV_CODE_DESIGN", []),
		    new HelpViewerNode("HV_CODE_ORGANIZATION", []),
		    new HelpViewerNode("HV_CODE_DEPENDENCIES", []),
		    new HelpViewerNode("HV_CODE_IO", []),
		    new HelpViewerNode("HV_CODE_FUTURE", []),
		    new HelpViewerNode("HV_E_T_MAX_REACHED", []),
		    new HelpViewerNode("HV_E_NOT_IMPLEMENTED", []),

		]),

		new HelpViewerNode("HV_CONCEPTS", [

		    new HelpViewerNode("HV_ENSEMBLES", []),
		    new HelpViewerNode("HV_MARKOV_ONE_STEP", []),
		    new HelpViewerNode("HV_SPIN_SYSTEMS", []),
		    new HelpViewerNode("HV_PHASE_TRANSITIONS", []),

		    new HelpViewerNode("HV_SA_SP", [

			new HelpViewerNode("HV_ST_RW", [

			    new HelpViewerNode("HV_P_SP_RW_l", []),
			    new HelpViewerNode("HV_P_SP_RW_r", []),
			    new HelpViewerNode("HV_P_SP_RW_N", []),
			    new HelpViewerNode("HV_P_SP_RW_x_0", []),

			]),

			new HelpViewerNode("HV_ST_MN", [

			    new HelpViewerNode("HV_P_SP_MN_mu", []),
			    new HelpViewerNode("HV_P_SP_MN_s", []),
			    new HelpViewerNode("HV_P_SP_MN_N", []),
			    new HelpViewerNode("HV_P_SP_MN_x_0", []),

			]),

		    ]),

		    new HelpViewerNode("HV_SA_SM", [

			new HelpViewerNode("HV_ST_IS", [

			    new HelpViewerNode("HV_P_SM_IS_T", []),
			    new HelpViewerNode("HV_P_SM_IS_h", []),
			    new HelpViewerNode("HV_P_SM_IS_N", []),

			]),

			new HelpViewerNode("HV_ST_XY", [

			    new HelpViewerNode("HV_P_SM_XY_T", []),
			    new HelpViewerNode("HV_P_SM_XY_h", []),
			    new HelpViewerNode("HV_P_SM_XY_N", []),

			]),

		    ]),

		    new HelpViewerNode("HV_SA_ND", [

			new HelpViewerNode("HV_ST_LM", [

			    new HelpViewerNode("HV_P_ND_LM_r", []),
			    new HelpViewerNode("HV_P_ND_LM_x_0", []),

			]),

			new HelpViewerNode("HV_ST_GM", [

			    new HelpViewerNode("HV_P_ND_GM_x_0", []),
			    new HelpViewerNode("HV_P_ND_GM_y_0", []),

			]),

		    ]),

		    new HelpViewerNode("HV_SA_FD", [

			new HelpViewerNode("HV_ST_EU", []),
			new HelpViewerNode("HV_ST_PF", []),

		    ]),
		]),

	    ]);

	// use hv_root structure to recursively fill data fields, create map, etc.
	this.hvn_map = {};
	hv_root.recurse_and_fill_data(this.hvn_map, []);
	//console.log(this.hvn_map);///////////
	if (this.show_on_load) {  // if indicated, show HelpViewer on app loading
	    this.show_view(this.initial_view);
	}
    }

    help_btn_pressed() {  // remember last shown view, or show Home

	if (this.curr_view == "") {
	    this.show_view('HV_HOME');
	} else {
	    $("#md_container").modal("show");
	}
    }

    show_view(v) {

	this.sim.rs.set_NR();  // exit from any running mode while using help viewer

	// hide outgoing view if it is showing
	if (this.curr_view != "") {
	    $("#" + this.curr_view).hide();
	}
	$('#md_container').modal('handleUpdate');  // reload modal to "scroll to top" of page about to be shown
	
	// update appearance/behavior of Home icon
	if (v == "HV_HOME") {
	    $("#hv_home_icon").removeClass("hv_link");
	    $("#hv_home_icon").addClass("hv_home_curr");
	} else {
	    $("#hv_home_icon").removeClass("hv_home_curr");
	    $("#hv_home_icon").addClass("hv_link");
	}

	// show incoming view (and modal, if it's hidden)
	let hvn = this.hvn_map[v];
	CU.sh("md_breadcrumbs", hvn.breadcrumbs_html);
	$("#" + v).show();
	if ( ! this.deployed) $("#md_container").modal("show");

	// if we're about to gain our one step of "history," enable back button (a one-time transition)
	if ((this.prev_view == "") && (this.curr_view != "")) {
	    $("#hv_back_btn_link").removeClass("hv_disabled_link");
	    $("#hv_back_btn_link").addClass("hv_link");
	}

	// update view variables and back button link
	this.prev_view = this.curr_view;  // save outgoing view to update back button
	this.curr_view = v;  // save incoming view
	if (this.prev_view != "") {  // update back button
	    $("#hv_back_btn_link").attr("onclick", "window.sim.ui.hv.show_view('" + this.prev_view + "'); ");
	}

	//console.log(hvn.get_all_child_links());  ////////////////
    }
}
