
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

        this.collisionEnabled = false;
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

        console.log("TIME -> ", timeStep);

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

    clone() {
    return new Atom(this.x, this.y, this.vx, this.vy, this.radius, this.mass);
}

}
class Particles {
    constructor(coords) {
        this.N = Params_IG.N.v;
        this.T = Params_IG.T.v;
        this.V = Params_IG.V.v;
        this.S = Params_IG.particleDisplaySize;
        this.M = Params_IG.particleMass;
        this.R = Params_IG.R;
        this.particles = [];
        //this.createParticles(this.N, coords.mc.unif01_rng, this.particleDisplaySize, this.particleMass);
    }

    handleCollisions() {
        if (typeof interactionType !== 'undefined' && interactionType) {
            console.warn("HANDLE COLLISIONS CALLED");
            for (let i = 0; i < this.particles.length; i++) {
                for (let j = i + 1; j < this.particles.length; j++) {
                    const p1 = this.particles[i];
                    const p2 = this.particles[j];

                    if (p1.checkCollision(p2)) {
                        p1.resolveCollision(p2);
                    }
                }
            }
        }
    }

    createParticles(n, func, particleSize = .025, particleMass = 1) {
       const minDistance = particleSize / 2;

        for (let i = 0; i < n; i++) {              
            const angle = func() * Math.PI * 2;   // Random direction (angle) using the prng

            //assign a random velocity based on the temp using boltzman dist. 
            const boltz = Math.sqrt(1.38064852e-23 * Params_IG.T.v / particleMass);

            const velo = Params_IG.gaussian(func, 0, boltz);
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
            } while (!this.isPositionValid(x, y, minDistance));

            const atom = new Atom(
                    x, y, vx, vy, particleSize, particleMass // x, y, vx, vy, radius, mass
                );
                this.particles.push(atom);
        }
    }

    isPositionValid(x, y, minDistance) {
        const box = Params_IG.boxSize;

        if (Math.abs(x) > box / 2 - minDistance || Math.abs(y) > box / 2 - minDistance) {
            return false;
        }

        for (const p of this.particles) {
            const dx = x - p.x;
            const dy = y - p.y;
            const dist = Math.sqrt(dx * dx + dy * dy);
            if (dist < minDistance) {
                return false;
            }
        }
        return true;
    }

    calculateVelocityDistribution() {
        const velocities = this.particles.map(p => Math.sqrt(p.vx ** 2 + p.vy ** 2));
        const freq = [];
        if (velocities.length === 0) return { velocities: [], freq: [] };

        const maxV = Math.max(...velocities);
        const minV = Math.min(...velocities);
        const binCount = 50;
        const binSize = (maxV - minV) / binCount;

        for (let i = 0; i < binCount; i++) freq[i] = 0;

        for (let v of velocities) {
            const index = Math.min(Math.floor((v - minV) / binSize), binCount - 1);
            freq[index]++;
        }

        const binCenters = [];
        for (let i = 0; i < binCount; i++) {
            binCenters[i] = minV + binSize * (i + 0.5);
        }

        return { velocities: binCenters, freq };
    }

    redistributeThermalEnergy(rngFunc = Math.random) {
        const k = 1.380649e-23;
        const m = Params_IG.particleMass;
        const targetTemp = Params_IG.T.v;

        const targetE = this.N * k * targetTemp;

        // Compute current total kinetic energy
        let currentE = 0;
        for (let p of this.particles) {
            const v2 = p.vx ** 2 + p.vy ** 2;
            currentE += 0.5 * m * v2;
        }

        const deltaE_total = targetE - currentE;

        // If T = 0 or target energy is 0, zero all velocities
        if (targetE <= 0) {
            for (let p of this.particles) {
                p.vx = 0;
                p.vy = 0;
            }
            return;
        }

        // Generate Gaussian weights (positive-only)
        let weights = [];
        let totalWeight = 0;
        for (let i = 0; i < this.N; i++) {
            let w = Math.abs(Params_IG.gaussian(rngFunc));
            weights.push(w);
            totalWeight += w;
        }

        // Apply energy change to each particle
        for (let i = 0; i < this.N; i++) {
            const p = this.particles[i];
            const angle = Math.atan2(p.vy, p.vx);
            const v2 = p.vx ** 2 + p.vy ** 2;
            const deltaE_i = (weights[i] / totalWeight) * deltaE_total;
            const v2_new = v2 + (2 * deltaE_i / m);
            const v_new = Math.sqrt(Math.max(v2_new, 0));
            p.vx = v_new * Math.cos(angle);
            p.vy = v_new * Math.sin(angle);
        }
    }



}


class Params_IG extends Params {

    static T = undefined;
    static V = undefined;
    static N = undefined;
    static R = 0.025;
    static timeStep = 1.0 / (1000 * 1000);
    static boxSize = 1;
    static total_time = 0;
    static sim_pressure = 0;
    static particleDisplaySize = 0.025;
    static particleMass = MassType.air;
    static expectedVelocity = 0;

    push_vals_to_UI() {
        Params_IG.T.push_to_UI(this.T);
    }

    get_info_str() {
        return "T = " + this.T;
    }

    static gaussian(func, mean = 0, stdDev = 1) {
        let u1 = func(); 
        let u2 = func();
        let randomStdNormal = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
        return mean + stdDev * randomStdNormal;
    }
}

class Coords_IG extends Coords {
    constructor(...args) {
        super(...args);

        const numParticles = Params_IG.N.v;
        const particleDisplaySize = Params_IG.particleDisplaySize;
        const particleMass = Params_IG.particleMass;
        this.temperature = Params_IG.T.v;
        const boundaryType = false; // bounce mode
        let sum = 0;

        if (this.constructing_init_cond) {
            // Initial condition
            this.particles = new Particles();
            this.particles.createParticles(numParticles, this.mc.unif01_rng, particleDisplaySize, particleMass);
            Params_IG.total_pressure = 0;
            Params_IG.total_time = 0;
            console.warn("CREATING NEW PARTICLES (init)");
        } else {
            // Evolution step: clone particles from c_prev
            this.particles = new Particles();
            this.particles.particles = this.c_prev.particles.particles.map(atom => atom.clone());  // deep clone

            for (let i = 0; i < numParticles; i++) {
                // If N increased
                 sum += Math.sqrt(this.particles.particles[i].vx ** 2 + this.particles.particles[i].vy ** 2);
                if (this.particles.particles[i] === undefined) {
                    this.particles.createParticles(1, this.mc.unif01_rng, particleDisplaySize, particleMass);
                }

                Params_IG.total_time += Params_IG.timeStep;
                this.particles.particles[i].updatePosition(Params_IG.timeStep, Params_IG.boxSize, boundaryType);
            }
            Params_IG.expectedVelocity = sum/Params_IG.N.v;

            console.log("Temp: ", Params_IG.T.v, " | N: ", Params_IG.N.v, " | Volume: ", Params_IG.V.v, " | Sim_Pressure: ", Params_IG.sim_pressure, 
                        " | Expected Velocity: ", Params_IG.expectedVelocity);
            if (this.mc.collisionEnabled) {
                    this.particles.handleCollisions();
                }

            if(this.c_prev.temperature != this.temperature) this.particles.redistributeThermalEnergy(this.mc.unif01_rng);

            Params_IG.sim_pressure = Params_IG.total_pressure / (Params_IG.total_time * Params_IG.boxSize);
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
