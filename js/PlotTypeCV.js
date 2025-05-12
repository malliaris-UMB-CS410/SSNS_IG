
///////////////////////////////////////////////////////////////////////////////////////////////
////////  PlotTypeCV_*  ///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////////  a PlotType that is built on the HTML Canvas "CV" element using basic shapes, etc.  //
///////////////////////////////////////////////////////////////////////////////////////////////

class PlotTypeCV extends PlotType {

    constructor() {

        super();

        this.canv_dim = undefined; // "declaration" (set within)
    }

    get_plot_width() {
        return this.canv_dim;
    }

    get_plot_height() {
        return this.canv_dim;
    }

    setup_canvas() {
        this.cc = document.getElementById(this.get_html_targ_id_str()).getContext("2d"); // cc = canvas context, for plotting calls
        this.plotCanvas = document.getElementById(this.get_html_targ_id_str());
        this.plotCanvas.style.border = "1px solid black";
        $("#" + this.get_html_targ_id_str()).attr("width", this.canv_dim); // NOTE: CV canvas plot_targets use **attr* not *css* for h/w
        $("#" + this.get_html_targ_id_str()).attr("height", this.canv_dim); // NOTE: CV canvas plot_targets use **attr* not *css* for h/w
    }

    clear_canvas() {
        this.cc.clearRect(0, 0, this.canv_dim, this.canv_dim);
    }

    // rtoa = relative to absolute; utility to transform coordinate from value on [0, 1] to [0, canvas_size]
    rtoa(rc) {
        return this.canv_dim * rc;
    }

    // ator = absolute to relative; utility to transform coordinate from value on [0, canvas_size] to [0, 1]
    ator(ac) {
        return ac / this.canv_dim;
    }

    // fyc = flip y coordinate; utility to transform (absolute) y coordinate so that the origin for x,y is lower left corner
    fyc(yc) {
        return this.canv_dim - yc;
    }
}

class PlotTypeCV_IG extends PlotTypeCV {

    constructor(trj) {
        super();

        console.warn("PLOT TYPE CONSTRUCTOR");
        this.trj = trj;
        this.display = "ATOM";
        this.canv_dim = parseInt(document.getElementById('UI_P_SM_IG_V').value) * 100;
        this.setup_canvas();

        this.vInput = document.getElementById("UI_P_SM_IG_V");
        if (this.vInput) {
            this.updateCanvasSize();
            this.vInput.addEventListener("input", () => {
                this.updateCanvasSize();
                this.setup_canvas();
                this.plot(this.trj.t);
            });
        }

        this.collisionCheckbox = document.getElementById("UI_P_SM_IG_C");
        if (this.collisionCheckbox) {
            this.collisionCheckbox.addEventListener("change", this.handleCollisionToggle.bind(this));
        }

        this.velocityDistributionCheckbox = document.getElementById("UI_P_SM_IG_VD");
        if (this.velocityDistributionCheckbox) {
            this.velocityDistributionCheckbox.addEventListener("change", this.handleVelocityDistributionToggle.bind(this));
        }
    }

    handleCollisionToggle(event) {
        const isChecked = event.target.checked;
        window.interactionType = isChecked;

        if (isChecked) {
            const particlesObj = this.trj.get_x(this.trj.t)?.particles;
            if (particlesObj && typeof particlesObj.handleCollisions === "function") {
                particlesObj.handleCollisions();  // ✅ No argument
            }
        }
    }

    handleVelocityDistributionToggle(event) {
        const isChecked = event.target.checked;
        window.interactionType = isChecked;
        if (isChecked) {
            console.warn("CHECKED FOR DISPLAY VELOCITY DISTRIBUTION");
            this.display = "Distribution";
        } else {
            this.display = "ATOM";
            console.warn("VELOCITY DISTRIBUTION DISABLED");
        }
    }

    updateCanvasSize() {
        const vValue = parseFloat(this.vInput.value);
        if (!isNaN(vValue) && vValue > 0) {
            this.canv_dim = Math.floor(vValue * 100);
        } else {
            this.canv_dim = PlotType.square_plot_width;
        }
    }

    draw_circle(x, y, radius) {
        let scaleX = this.rtoa(x);
        let scaleY = this.rtoa(y);
        let scaleR = this.rtoa(radius);
        this.cc.beginPath();
        this.cc.arc(scaleX, scaleY, scaleR, 0, Math.PI * 2);
        this.cc.fillStyle = 'red';
        this.cc.fill();
    }

    draw_bar(x, height, bar_width, bar_height){
    this.cc.fillStyle = "blue";
    this.cc.fillRect(x, height, bar_width, bar_height);
    }

    get_ext_x_axis_lbl_str() {
        return "x";
    }

    get_ext_y_axis_lbl_str() {
        return "y";
    }

    get_html_targ_id_str() {
        return "plot_CV_IG";
    }

    update_canvas(t) {
        this.clear_canvas();
        
        if (this.display == "ATOM") {
            console.warn("SHOULD DISPLAY ATOMS IN CANVAS");
            for (let i = 0; i < Params_IG.N.v; i++) {
                let cp = this.trj.get_x(t).particles.particles[i];

                let halfBox = Params_IG.boxSize / 2;
                let normX = (cp.x + halfBox) / Params_IG.boxSize;
                let normY = (cp.y + halfBox) / Params_IG.boxSize;
                this.draw_circle(normX, normY, Params_IG.R);
            }
        } else {
            let { velocities, freq } = this.trj.get_x(t).particles.calculateVelocityDistribution();
            console.warn("SHOULD DISPLAY FREQUENCY IN CANVAS | Total: ", velocities.length, " - ", freq.length);

            if (velocities.length && freq.length) {
                const canvasW = this.cc.canvas.width;
                const canvasH = this.cc.canvas.height;

                const maxFrequency = Math.max(...freq);
                const barWidth = Math.max(1, (canvasW * 0.8) / freq.length);
                const spacing = Math.max(1, (canvasW * 0.2) / (freq.length + 1));
                const minHeight = 2;  // Minimum visible bar height

                let x = spacing;

                for (let i = 0; i < freq.length; i++) {
                    let rawHeight = (freq[i] / maxFrequency) * (canvasH * 0.8);
                    const barHeight = Math.max(minHeight, rawHeight);
                    const y = canvasH - barHeight;

                    console.log(`DRAWING BAR: x=${x}, y=${y}, w=${barWidth}, h=${barHeight}`);
                    this.draw_bar(x, y, barWidth, barHeight);
                    x += barWidth + spacing;
                }
            }
        }
    }

    plot(t) {     
        this.update_canvas(t);   
    }
}

// NOTE: heatmaps are rendered on an html canvas and maintain state, so their get_html_targ_id_str() must be ST-dependent... explain this a bit more
class PlotTypeCV_Spin extends PlotTypeCV {

    constructor(trj) {

        super();

        this.trj = trj;
        if (!this.get_color_from_spin_val)
            throw new Error("Derived PlotTypeCV must define get_color_from_spin_val()");
        this.tile_dim = undefined; // "declaration" (set within)
        this.determine_tile_canv_dims(this.trj.mc.N);
    }

    determine_tile_canv_dims(N) {
        let rough_tile_dim = Math.floor(PlotType.square_plot_width / N); // "target" a plot dim that is just under PlotType.square_plot_width
        this.tile_dim = parseInt(Math.max(rough_tile_dim, 1.0)); // require at least 1 pixel per tile (which may put plot dim above target!)
        this.canv_dim = this.tile_dim * N;
    }

    set_pixel(x, y, v) {

        let xc = y * this.tile_dim; // xc = x on canvas; NOTE THE TRANSPOSE OPERATION y --> xc and x --> yc
        let yc = x * this.tile_dim; // yc = y on canvas; NOTE THE TRANSPOSE OPERATION y --> xc and x --> yc
        this.cc.fillStyle = this.get_color_from_spin_val(v);
        this.cc.fillRect(xc, yc, this.tile_dim, this.tile_dim);
    }

    draw_entire_canv_from_spin_arr(sa) {

        for (let x = 0; x < sa.shape[0]; x++) {
            for (let y = 0; y < sa.shape[1]; y++) {
                this.set_pixel(x, y, sa.get(x, y));
            }
        }
    }

    setup_spin_array() {
        let init_sa = this.trj.get_x(this.trj.t_0).spins; // spin array associated with initial time step t_0 used to initialize canvas
        this.draw_entire_canv_from_spin_arr(init_sa);
        this.last_t_displayed = this.trj.t_0;
    }

    update_canvas(t) {

        let moving_fwd = (this.last_t_displayed < t);
        if (moving_fwd) {
            for (let s = this.last_t_displayed; s < t; s++) {
                let nt = this.trj.get_x(s).next_trans; // nt = next transition object (for coords at time s)
                if (nt.move_occurred) {
                    this.set_pixel(nt.x, nt.y, nt.new_val);
                }
            }
        } else { // moving backward (or at t == this.last_t_displayed, in which case for loop will execute 0 times)
            for (let s = this.last_t_displayed; s > t; s--) {
                let pt = this.trj.get_x(s).prev_trans; // pt = previous transition object (for coords at time s)
                if (pt.move_occurred) {
                    this.set_pixel(pt.x, pt.y, pt.old_val);
                }
            }
        }
        this.last_t_displayed = t; // update for next call
    }

    get_ext_x_axis_lbl_str() {
        return "\\mathrm{spin} \\hspace{0.4em} x \\hspace{0.4em} \\mathrm{index}";
    }

    get_ext_y_axis_lbl_str() {
        return "\\mathrm{spin} \\hspace{0.4em} y \\hspace{0.4em} \\mathrm{index}";
    }

    plot(t) {
        this.update_canvas(t);
    }
}

class PlotTypeCV_IS extends PlotTypeCV_Spin {

    constructor(trj) {

        super(trj);

        this.trj = trj;
        this.setup_canvas();
        this.setup_spin_array();
    }

    get_html_targ_id_str() {
        return "plot_CV_IS"; // heatmaps are rendered on an html canvas and maintain state, so must be ST-dependent
    }

    get_color_from_spin_val(v) { // translates a numeric spin value into an html color string

        return ((v == 0) ? "hsl(29, 85%, 44%)" : "hsl(52, 100%, 51%)"); // dark orange hsl(29, 85%, 44%), bright yellow hsl(52, 100%, 51%)
    }
}

class PlotTypeCV_XY extends PlotTypeCV_Spin {

    constructor(trj) {

        super(trj);

        this.trj = trj;
        this.setup_canvas();
        this.setup_spin_array();
    }

    get_html_targ_id_str() {
        return "plot_CV_XY"; // heatmaps are rendered on an html canvas and maintain state, so must be ST-dependent
    }

    get_color_from_spin_val(v) { // translates a numeric spin value into an html color string

        let hue_val = v * 360.0; // convert fraction of revolution --> degrees
        let hsl_str = "hsl(" + hue_val + ", 55%, 50%)";
        return hsl_str;
    }
}
