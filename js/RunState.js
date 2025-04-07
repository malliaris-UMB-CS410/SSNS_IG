
// helper class for Simulator
class RunState {

    constructor(sim) {
	this.sim = sim;
	this.params_changed = false;  // i.e., a **user-initiated** change of parameters
	this.is_running_forward = false;  // dummy value
	this.is_running_backward = false;  // dummy value
	this.is_recording = false;  // dummy value
	this.set_NR();  // state is Not Running when instantiated
    }
    is_RN() {  // "is_running" (backward, forward, or recording)
	return (this.is_running_backward || this.is_running_forward || this.is_recording);
    }
    is_RB() {
	return this.is_running_backward;
    }
    is_RF() {
	return this.is_running_forward;
    }
    is_RC() {
	return this.is_recording;
    }
    set_RB() {
	this.is_running_backward = true;
    }
    set_RF() {
	this.is_running_forward = true;
    }
    set_RC() {
	this.is_recording = true;
    }
    set_NR() {  // NR = Not Running
	this.is_running_backward = false;
	this.is_running_forward = false;
	this.is_recording = false;
    }
}
