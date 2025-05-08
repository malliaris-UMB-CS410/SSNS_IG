
///////////////////////////////////////////////////////////////////////////////////////////////
////////  IG = Insane gas (from SM = Statistical Mechanics)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
const MassType = {
    air: 4.8*Math.pow(10,-26),
    helium: 6.6e-27,
    oxygen: 5.32e-23,
    nitrogen: 4.7e-26
};
// NOTE: we represent the two Ising spin values as 0,1 "under the hood" -- it's more computationally convenient to think of as binary
// and, in creating debugging output, etc., 0 and 1 have the same width; only in ModelCalc_IG.get_E_spin_pair() where we switch to energy
// quantities do we have to translate 0, 1 the to the more physically appropriate -1, 1 via as1n1(); all other methods dealing with energies
// take output from get_E_spin_pair(), and all methods that take spin values as inputs always expect 0, 1
class ModelCalc_IG extends ModelCalc {

    constructor() {
	super();

	this.unif01_rng = randu.factory({'seed': ModelCalc_Stoch.rng_seed.v });
	this.normal_rng = normal.factory({'seed': ModelCalc_Stoch.rng_seed.v });
	console.log("INFO:\tusing PRNG algorithm Mersenne Twister 19937 (the default) on all:", this.unif01_rng.PRNG.NAME, this.normal_rng.PRNG.NAME);
	console.log("INFO:\tusing seed value = " + ModelCalc_Stoch.rng_seed.v);
	console.log("INFO:\tNOTE: ModelCalc_IG **does not** extend ModelCalc_Stoch!  While PRNGs are used for initial condition, all time evolution is deterministic!");
    console.log(4*1.380649*Math.pow(10,-23)*1/Params_IG.boxSize);
	}
    model_is_stoch(){
		return false;
	}
}

class Atom {
		constructor(x, y, vx, vy, radius=0, mass=MassType.air) {
			this.x = x;
			this.y = y;
			this.vx = vx;
			this.vy = vy;
			this.radius = radius;
			this.mass = mass;
            this.xenergy = 0;
            this.yenergy = 0;
	}

	updatePosition(timeStep, boxSize, boundaryType) {
        //console.log(this.x, this.y, this.vx, this.vy, Math.random());
        let newX = this.x + this.vx * timeStep;
        let newY = this.y + this.vy * timeStep;

        const halfBox = boxSize / 2;

        if (boundaryType) {
            // Wrap mode
            if (newX < -halfBox) {
                newX = newX + boxSize;
            } else if (newX > halfBox) {
                newX = newX - boxSize;
            }

            if (newY < -halfBox) {
                newY = newY + boxSize;
            } else if (newY > halfBox) {
                newY = newY - boxSize;
            }
        } else {
            // Bounce mode
            const radiusOffset = this.radius;
            // Check X bounce
            if (newX - radiusOffset < -halfBox || newX + radiusOffset > halfBox) {
                this.vx *= -1;
                Params_IG.total_pressure += 2 * this.mass * Math.abs(this.vx);
                newX = Math.max(-halfBox + radiusOffset, Math.min(halfBox - radiusOffset, newX));
            }

            if (newY - radiusOffset < -halfBox || newY + radiusOffset > halfBox) {
                this.vy *= -1;
                Params_IG.total_pressure += 2 * this.mass * Math.abs(this.vy);

                // 防止粒子穿过边界后卡住
                newY = Math.max(-halfBox + radiusOffset, Math.min(halfBox - radiusOffset, newY));
            }
            //console.log(Params_IG.total_pressure);
            this.x = newX;
            this.y = newY;
        }
		//console.log("update position: x:", this.x)
    }

    

    // Calculate distance between two atoms
    static distance(atom1, atom2) {
        const dx = atom1.x - atom2.x;
        const dy = atom1.y - atom2.y;
        return Math.sqrt(dx * dx + dy * dy);
    }

    // Check for collision with another atom
    checkCollision(other) {
        // Calculate actual distance between particles in the simulation space
        const distance = Atom.distance(this, other);

        // Collision occurs when distance is less than sum of radii
        // We divide by 100 to convert from pixel space to simulation space
        return distance < (this.radius + other.radius) / 100;
    }

    // Handle collision physics with another atom
    resolveCollision(other) {
        // Calculate velocity difference and position difference
        const xVelocityDiff = this.vx - other.vx;
        const yVelocityDiff = this.vy - other.vy;

        const xDist = other.x - this.x;
        const yDist = other.y - this.y;

        // Prevent division by zero
        if (xVelocityDiff * xDist + yVelocityDiff * yDist >= 0) {
            // Calculate collision angle
            const angle = -Math.atan2(yDist, xDist);

            // Store masses in variables for readability
            const m1 = this.mass;
            const m2 = other.mass;

            // Calculate velocity magnitudes and angles before collision
            const u1 = this.rotateVector(this.vx, this.vy, angle);
            const u2 = this.rotateVector(other.vx, other.vy, angle);

            // Calculate velocity after collision using conservation of momentum
            // and elastic collision formula
            const v1 = {
                x: u1.x * (m1 - m2) / (m1 + m2) + u2.x * 2 * m2 / (m1 + m2),
                y: u1.y
            };

            const v2 = {
                x: u2.x * (m2 - m1) / (m1 + m2) + u1.x * 2 * m1 / (m1 + m2),
                y: u2.y
            };

            // Rotate velocities back to original coordinate system
            const vFinal1 = this.rotateVector(v1.x, v1.y, -angle);
            const vFinal2 = this.rotateVector(v2.x, v2.y, -angle);

            // Update velocities
            this.vx = vFinal1.x;
            this.vy = vFinal1.y;

            other.vx = vFinal2.x;
            other.vy = vFinal2.y;

            // Prevent particles from getting stuck together
            this.separateParticles(this, other);
        }
    }

    // Helper method to rotate vectors for collision calculation
    rotateVector(x, y, angle) {
        return {
            x: x * Math.cos(angle) - y * Math.sin(angle),
            y: x * Math.sin(angle) + y * Math.cos(angle)
        };
    }

    // Separate particles that are overlapping to prevent sticking
    separateParticles(particle1, particle2) {
        const distance = Atom.distance(particle1, particle2);
        const minDistance = (particle1.radius + particle2.radius) / 100;

        // Only separate if particles are overlapping
        if (distance < minDistance) {
            const overlap = (minDistance - distance) / 2;

            // Calculate separation vector
            const dx = particle2.x - particle1.x;
            const dy = particle2.y - particle1.y;

            // Normalize direction vector
            const magnitude = Math.sqrt(dx * dx + dy * dy) || 1; // Avoid division by zero
            const unitX = dx / magnitude;
            const unitY = dy / magnitude;

            // Move particles apart along separation vector
            particle1.x -= unitX * overlap;
            particle1.y -= unitY * overlap;
            particle2.x += unitX * overlap;
            particle2.y += unitY * overlap;
        }
    }
	
}
//const particleRadius = 6;       // set size of particle 
//const minDistance = 0;// (this.radius * 2) / 100; // Converted to simulation units
//let lastTimestamp = 0;
//const maxTimeStep = 1.0 / 30.0;; // 30 FPS
//const interactionType = 0; //parseInt(document.getElementById('interactionType').value); // interaction with particles
/*

// values
let seed = SeededRNG(1);
//const userBoxWidth = parseFloat(document.getElementById('boxWidth').value);
//const userBoxHeight = parseFloat(document.getElementById('boxHeight').value);
const userVelocity = 4; //parseFloat(document.getElementById('velocityMagnitude').value) / 10;
//let numParticles = Coords_IG.N.v;//  5; //parseInt(document.getElementById('numParticles').value);
const boundaryType = 0; //parseInt(document.getElementById('boundaryType').value);		// interaction with walls
const interactionType = 0; //parseInt(document.getElementById('interactionType').value); // interaction with particles
//const initialSeed = parseInt(document.getElementById('initialSeed').value);
const particles = [];

const minDistance = (particleRadius * 2) / 100; // Converted to simulation units
const userBoxHeight = 4;
const userBoxWidth = 4;
let animationId = null;
// replace with RAF optimize the animate playing
let lastTimestamp = 0;
const maxTimeStep = 1.0 / 30.0;; // 30 FPS



*/
// Check and resolve collisions between all particles
function handleCollisions() {
	// if interaction type is set to collision, collide. else ideal gas
	if (interactionType){
		// Check each pair of particles exactly once
		for (let i = 0; i < particles.length; i++) {
			for (let j = i + 1; j < particles.length; j++) {
				const particle1 = particles[i];
				const particle2 = particles[j];

				if (particle1.checkCollision(particle2)) {
					particle1.resolveCollision(particle2);
				}
			}
		}
	}
	//else do nothing (ideal gas)
}

function createParticles(n, particles, func, particleSize = .025, particleMass = 1) {
	// create numParticles number of particles in a random starting position moving in a random direction
    const minDistance = particleSize / 2;

	for (let i = 0; i < n; i++) {
		const temp = func();
        console.log("temp", temp);
        const angle = temp * Math.PI * 2; //getRandomAngle();  // Random direction (angle)
		
        //const angle = this.mc.unif01_rng() * Math.PI * 2; 

        console.log("angle ", angle);
        

			// Calculate the x and y components of the velocity based on the random angle
			const vx = Math.sqrt(2*1.380649*Math.pow(10,-23)*Params_IG.T.v/ MassType.air) * Math.cos(angle); // X velocity component
			const vy = Math.sqrt(2*1.380649*Math.pow(10,-23)*Params_IG.T.v/ MassType.air) * Math.sin(angle); // Y velocity component

		// Try to find a valid non-overlapping position
		let x = 0, y = 0;
		let attempts = 0;
		const maxAttempts = 1000;
		
		do {
			//x = getRandomPosition(Params_IG.boxWidth / 2);
			//y = getRandomPosition(Params_IG.boxHeight /2);
		        
            //x = func() * Params_IG.boxWidth - (Params_IG.boxWidth / 2);
            //y = func() * Params_IG.boxHeight - (Params_IG.boxHeight / 2);

                x = func() * Params_IG.boxSize - (Params_IG.boxSize / 2);
				y = func() * Params_IG.boxSize - (Params_IG.boxSize / 2);
				attempts++;
		
				// Break out after too many attempts to prevent infinite loop
				if (attempts > maxAttempts) {
					console.warn('Could not find non-overlapping position after', maxAttempts, 'attempts');
					x = 0;
					y = 0;
					break;
				}
			} while (!isPositionValid(x, y, particles, minDistance, Params_IG.boxSize));
		
		console.log(`Creating particle: x=${x}, y=${y}, vx=${vx}, vy=${vy}`);

			const atom = new Atom(
				x, y, vx, vy, particleSize, MassType.air // x, y, vx, vy, radius, mass
			);
			//particles.push(atom);
			//console.log("create particles end:", atom);
            particles.push(atom);
            //console.log("After push: particles[n].x =", particles[particles.length - 1].x, "particles[n].y =", particles[particles.length - 1].y);
	    }
}


// Initialize with random starting positions
function getRandomPosition(range) {
	//console.log('getRandomPosition range:', userBoxWidth / 2, userBoxHeight / 2);
    const randomVar = seed.nextFloat() * range * 2 - range;
    console.log('randomVar:', randomVar);
    return randomVar;
}

// Generate random angle (in radians) between 0 and 2 * PI (360)
function getRandomAngle() {
    return seed.nextFloat() * Math.PI * 2; // Random angle between 0 and 2π (360 degrees)
}

// Helper function to ensure particles don't overlap initially
function isPositionValid(x, y, particles, minDistance, boxSize) {
	if(Math.abs(x) > boxSize / 2 - minDistance || Math.abs(y) > boxSize / 2 - minDistance) {
		return false;
	}
    console.log(particles);
	for (const particle of particles) {
		const dx = x - particle.x;
		const dy = y - particle.y;
		const distance = Math.sqrt(dx * dx + dy * dy);
		if (distance < minDistance) {
			return false;
		}
	}
	return true;
}

class Params_IG extends Params {

    static T = undefined;  // = new UINI_float(this, "UI_P_SM_IG_T", false);  assignment occurs in UserInterface(); see discussion there
	static V = undefined;  // = new UINI_float(this, "UI_P_SM_IG_V", false);  assignment occurs in UserInterface(); see discussion there			!!!! not sure this is right -jg !!!!
    static N = undefined;  // = new UINI_int(this, "UI_P_SM_IG_N", false);  assignment occurs in UserInterface(); see discussion there
    static timeStep = 1.0 / (1000*1000);
    static boxSize = 1;
    //static total_pressure = 0;
    static total_time = 0;
    static sim_pressure = 0;
    push_vals_to_UI() {
		Params_IG.T.push_to_UI(this.T);
    }

	get_info_str() {
		return "T = " + this.T;
	}
}

/*
function animate(timestamp) {
	if (!lastTimestamp) lastTimestamp = timestamp;

	const deltaTime = timestamp - lastTimestamp;
	

	// limited 30 FPS
	const timeStep = Math.min(deltaTime / 1000, maxTimeStep);

	handleCollisions();

	//console.log("animate: deltatime:", deltaTime);
	//console.log("animate: timestep:", timeStep);
	//console.log("animate: timestamp:", timestamp);
	//console.log("animate: maxtimestep:", maxTimeStep);

	// update particles position

	
	Coords_IG.particles.forEach(atom => {
        atom.updatePosition(timeStep, userBoxWidth, userBoxHeight, boundaryType);

            if (!boundaryType) {
	            const radiusOffset = atom.radius / 100;
	            atom.x = Math.max(-userBoxWidth / 2 + radiusOffset, Math.min(userBoxWidth / 2 - radiusOffset, atom.x));
	            atom.y = Math.max(-userBoxHeight / 2 + radiusOffset, Math.min(userBoxHeight / 2 - radiusOffset, atom.y));
            }
        }
    );    

	//draw();
	lastTimestamp = timestamp;
	animationId = requestAnimationFrame(animate);
}
*/
//let timeStep = 0; 
class Coords_IG extends Coords {

    //static N = undefined;  // = new UINI_Int(this, "UI_P_SM_IG_N", false);  assignment occurs in UserInterface(); see discussion there
	

    constructor(...args) {  // see discussion of # args at definition of abstract Coords()

	    super(...args);
	
    	/// TEMP CODE ///
	    //const computedVelocity = 5; // !!!!!!!!!!!!! This is TEMP code until max boltz is implemented !!!!!!!!!!!!!!!!
	    //const boxSize = 4;			// !!!!!!!!!!!!! TEMP code until we have a uniform size

	    let numParticles = Params_IG.N.v;
    	let tempK = Params_IG.T.v;
	    //console.log("IG.js coords_IG: numParticles", numParticles);
	    //console.log("IG.js coords_IG: temp", tempK); 
	    //console.log(timestamp);
	    //console.log("Coords_IG");
	    //console.log(particles);


	

    	if (this.constructing_init_cond) {
	    	//console.log("Coords_IG if");
		    //console.log("numParticles:", numParticles);
            this.particles = [];
            /*
            for (let i = 0; i < numParticles; i++) {
		        const temp = this.mc.unif01_rng();
                this.particles.push(new Atom(0.1 + i * 0.1 + 0.1*temp, 0.1 + i * 0.1, 0.3, 0.4, 0.005, 1));
                console.log("temp: ", temp);
            }
            console.log(this.particles);
            */
            //console.log("CIG:", this.particles);

		    createParticles(numParticles, this.particles, this.mc.unif01_rng);

    		//console.log("particles created:", JSON.stringify(particles, null, 2));
            //console.log("After created: particles[n].x =", particles[particles.length - 1].x, "particles[n].y =", particles[particles.length - 1].y);
            Params_IG.total_pressure = 0;
            Params_IG.total_time = 0;
	    } else {
		    //console.log("Coords_IG else", this.c_prev.particles[0].x);
            this.particles = copy(this.c_prev.particles);
	    	//console.log(particles);
		    for (let i = 0; i < numParticles; i++) {
			    if (this.particles[i] == undefined) {		// If N is increased by the user create a new particle
				    createParticles(1, this.particles, this.mc.unif01_rng);
    			} 
			
    			//console.log(i, this.particles[i]);
                //console.log(i, this.particles[i].x, this.particles[i].y, this.particles[i].vx, this.particles[i].vy);
		    	Params_IG.total_time += Params_IG.timeStep;
                this.particles[i].updatePosition(Params_IG.timeStep, Params_IG.boxSize, false);
                Params_IG.sim_pressure = Params_IG.total_pressure / (Params_IG.total_time * Params_IG.boxSize);
                console.log(Params_IG.sim_pressure);		
		    }
    		//timeStep++;
	        // this.x = this.mc.get_x_new(this.p, this.c_prev.x);
		
		}
		//animate();
    }

    //output() {
	//console.log("x =", this.x);
    //}


}

class Trajectory_IG extends Trajectory {

    constructor(sim) {
	super(sim);
    }

    gmc() {  // gmc = get ModelCalc object
	return new ModelCalc_IG();
    }

    gp() {  // gp = get Params object
	return new Params_IG(Params_IG.T.v); // T is uini object .v is the value
    }

    gc_ic(mc) {  // gc_ic = get Coords, initial condition
	return new Coords_IG(mc, [ Params_IG.N.v ]);
    }

    gc_nv(mc, p, c_prev) {  // gc_nv = get Coords, new value
	return new Coords_IG(mc, p, c_prev, []);
    }

    get_max_num_t_steps() {
	return Trajectory.DEFAULT_MAX_NUM_T_STEPS
    }
}


