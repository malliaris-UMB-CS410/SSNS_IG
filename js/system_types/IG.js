
///////////////////////////////////////////////////////////////////////////////////////////////
////////  IG = Insane gas (from SM = Statistical Mechanics)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
const MassType = {
    air: 4.8 * Math.pow(10, -26),
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

        this.unif01_rng = randu.factory({
            'seed': ModelCalc_Stoch.rng_seed.v
        });
        this.normal_rng = normal.factory({
            'seed': ModelCalc_Stoch.rng_seed.v
        });
        console.log("INFO:\tusing PRNG algorithm Mersenne Twister 19937 (the default) on all:", this.unif01_rng.PRNG.NAME, this.normal_rng.PRNG.NAME);
        console.log("INFO:\tusing seed value = " + ModelCalc_Stoch.rng_seed.v);
        console.log("INFO:\tNOTE: ModelCalc_IG **does not** extend ModelCalc_Stoch!  While PRNGs are used for initial condition, all time evolution is deterministic!");
        console.log(4 * 1.380649 * Math.pow(10, -23) * 1 / Params_IG.boxSize);
    }
    model_is_stoch() {
        return false;
    }
}

class Atom {
    constructor(x, y, vx, vy, radius = 0, mass = MassType.air) {
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
            this.x = newX;
            this.y = newY;
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
                newY = Math.max(-halfBox + radiusOffset, Math.min(halfBox - radiusOffset, newY));
            }
            this.x = newX;
            this.y = newY;
        }
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
        return distance < (this.radius + other.radius);
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

// Check and resolve collisions between all particles
function handleCollisions(particles) {
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

// function to to make the 2D boltz dist to sample the prng velocity 
function gaussian(func, mean = 0, stdDev = 1) {
    // this uses the box-muller transform
    let u1 = func(); 
    let u2 = func();

    let randomStdNormal = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);

    return mean + stdDev *  randomStdNormal;
}



function createParticles(n, particles, func, particleSize = .025, particleMass = 1) {
    // create numParticles number of particles in a random starting position moving in a random direction
    const minDistance = particleSize / 2;

    for (let i = 0; i < n; i++) {              
        const angle = func() * Math.PI * 2;   // Random direction (angle) using the prng

        //assign a random velocity based on the temp using boltzman dist. 
        const boltz = Math.sqrt(1.38064852e-23 * Params_IG.T.v / particleMass);

        const velo = gaussian(func, 0, boltz);
        const vx = velo * Math.cos(angle);
        const vy = velo * Math.sin(angle);

        /*
        // Calculate the x and y components of the velocity based on the random angle
        const vx = Math.sqrt(2 * 1.380649 * Math.pow(10, -23) * Params_IG.T.v / particleMass) * Math.cos(angle); // X velocity component
        const vy = Math.sqrt(2 * 1.380649 * Math.pow(10, -23) * Params_IG.T.v / particleMass) * Math.sin(angle); // Y velocity component
        */

        // Try to find a valid non-overlapping position
        let x = 0,
        y = 0;
        let attempts = 0;
        const maxAttempts = 1000;

        do {
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

        const atom = new Atom(
                x, y, vx, vy, particleSize, particleMass // x, y, vx, vy, radius, mass
            );
        particles.push(atom);
    }
}

// Helper function to ensure particles don't overlap initially
function isPositionValid(x, y, particles, minDistance, boxSize) {
    if (Math.abs(x) > boxSize / 2 - minDistance || Math.abs(y) > boxSize / 2 - minDistance) {
        return false;
    }
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

    static T = undefined;
    static V = undefined;
    static N = undefined;
    static timeStep = 1.0 / (1000 * 1000);
    static boxSize = 1;
    static total_time = 0;
    static sim_pressure = 0;
    push_vals_to_UI() {
        Params_IG.T.push_to_UI(this.T);
    }

    get_info_str() {
        return "T = " + this.T;
    }
}

class Coords_IG extends Coords {

    constructor(...args) { // see discussion of # args at definition of abstract Coords()

        super(...args);

        let numParticles = Params_IG.N.v;
        let tempK = Params_IG.T.v;

        const boundaryType = false;      // bounce  = false 
 


        let particleDisplaySize = 0.025;
        let particleMass = MassType.air;

        if (this.constructing_init_cond) {
            this.particles = [];
            createParticles(numParticles, this.particles, this.mc.unif01_rng, particleDisplaySize, particleMass);
            Params_IG.total_pressure = 0;
            Params_IG.total_time = 0;
        } else {
            this.particles = copy(this.c_prev.particles);
            for (let i = 0; i < numParticles; i++) {
                if (this.particles[i] == undefined) { // If N is increased by the user create a new particle
                    createParticles(1, this.particles, this.mc.unif01_rng, particleDisplaySize, particleMass);
                }

                Params_IG.total_time += Params_IG.timeStep;
                this.particles[i].updatePosition(Params_IG.timeStep, Params_IG.boxSize, boundaryType);

                if (Params_IG.C.v) {
                    handleCollisions(this.particles);
                }
            }
            Params_IG.sim_pressure = Params_IG.total_pressure / (Params_IG.total_time * Params_IG.boxSize);
            console.log(Params_IG.sim_pressure);
        }
    }
}

class Trajectory_IG extends Trajectory {

    constructor(sim) {
        super(sim);
    }

    gmc() { // gmc = get ModelCalc object
        return new ModelCalc_IG();
    }

    gp() { // gp = get Params object
        return new Params_IG(Params_IG.T.v); // T is uini object .v is the value
    }

    gc_ic(mc) { // gc_ic = get Coords, initial condition
        return new Coords_IG(mc, [Params_IG.N.v]);
    }

    gc_nv(mc, p, c_prev) { // gc_nv = get Coords, new value
        return new Coords_IG(mc, p, c_prev, []);
    }

    get_max_num_t_steps() {
        return Trajectory.DEFAULT_MAX_NUM_T_STEPS
    }
}
