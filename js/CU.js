
///////////////////////////////////////////////////////////////////////////////////////////////
////////  UTILITY CLASSES  ////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////////  classes that are handy but have no strong relationship to the overall project  //////
///////////////////////////////////////////////////////////////////////////////////////////////


// convenience utility (CU) methods -- all static, for easy access without instantiation
// NOTE: for those that are CSS/jQuery-based and take string identifier argument, leave leading "#" off -- it will be added!

class CU {

    static gs(str_id) {  // gs = get string
	return $("#" + str_id).val();
    }
    static gi(str_id) {  // gi = get int
	return parseInt($("#" + str_id).val());
    }
    static gf(str_id) {  // gf = get float
	return parseFloat($("#" + str_id).val());
    }
    static sv(str_id, new_val) {  // s = set value (int or float input, but stored as a string I think)
	$("#" + str_id).val(new_val);
    }
    static gmn(str_id) {  // gmn = get min
	return $("#" + str_id).attr("min");
    }
    static gmx(str_id) {  // gmx = get max
	return $("#" + str_id).attr("max");
    }
    static smn(str_id, new_min_val) {  // smn = set min
	$("#" + str_id).attr("min", new_min_val);
    }
    static smx(str_id, new_max_val) {  // smx = set max
	$("#" + str_id).attr("max", new_max_val);
    }
    static sh(str_id, new_html) {  // sh = set html, e.g., between open/close pair of html tags
	$("#" + str_id).html(new_html);
    }
    static sk(str_id, tex_str) {  // sk = set katex, i.e., render a dynamically generated LaTeX string
	katex.render(tex_str, document.getElementById(str_id), {throwOnError: false});  // heads up: only works with js, not jQuery
    }
    static gcb(str_id) {  // gcb = get check box value
	return $("#" + str_id).is(":checked");
    }
    // NOTE: grb() arg is **name**, while srb() arg is **id** in, e.g., <input type="radio" name="UI_ST" id="UI_ST_ND_LM">
    static grb(str_id) {  // grb = get radio button value
	return $("input[name='" + str_id + "']:checked").val();
    }
    // NOTE: grb() arg is **name**, while srb() arg is **id** in, e.g., <input type="radio" name="UI_ST" id="UI_ST_ND_LM">
    static srb(str_id) {  // srb = set radio button value
	$("#" + str_id).prop("checked", true);
    }
    static eb(str_id) {  // NEEDED?  eb = enable button
    	$("#" + str_id).removeClass("disabled");  // NOTE: enabling an already enabled button is not an error -- it has no net effect
    }
    static db(str_id) {  // NEEDED?  db = disable button
    	$("#" + str_id).addClass("disabled");  // NOTE: disabling an already disabled button is not an error -- it has no net effect
    }
    static gmv(va, vb, vc) {  // gmv = get middle value, i.e., the one that's between the other two numerically
    	let arr_to_sort = [ va, vb, vc ];  // assemble array of all three
	return arr_to_sort.sort( (U,V) => U - V )[1];  // after sorting numerically, the middle value (index [1]) is the one we want...
    }
    static add_to_entry_OM(om, k, v) {  // for OrderedMap om (from js-sdsl library), if om[k] exists, do om[k] += v, else insert om[k] = v
	if (v > 0) {  // don't do anything if v == 0
	    let orig_v = om.getElementByKey(k);
	    if (orig_v == undefined) {
		om.setElement(k, v);
	    } else {
		om.setElement(k, orig_v + v);
	    }
	}
    }
    static print_stats_OM(om) {
	let sum = 0;
	om.forEach((element, index) => {
	    sum += element[1];
	});
	console.log("OM x_min, x_max, x_span_p1, #_x_entries, sum:\t", om.front()[0], om.back()[0], (om.back()[0] - om.front()[0] + 1), om.size(), sum);
    }
    static ppOM(om) {  // ppOM = pretty print OrderedMap (from js-sdsl library)
	om.forEach((element, index) => {
	    console.log("ppOM:\t", element[0], element[1], index);
	});
    }
}
