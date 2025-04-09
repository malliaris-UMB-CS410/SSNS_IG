class Atom {
    constructor(x, y, vx, vy, radius=5, mass=1) {
        this.x = x;
        this.y = y;
        this.vx = vx;
        this.vy = vy;
        this.radius = radius;
        this.mass = mass;
        this.m = mass;
    }

    updatePosition(timeStep, BoxSize, boundaryType) {
        let newX = this.x + this.vx * timeStep;
        let newY = this.y + this.vy * timeStep;

        if (boundaryType) {
            // Wrap mode
            if (newX < -BoxSize / 2) {
                this.x = newX + BoxSize;
            } else if (newX > BoxSize / 2) {
                this.x = newX - BoxSize;
            } else {
                this.x = newX;
            }

            if (newY < -BoxSize / 2) {
                this.y = newY + BoxSize;
            } else if (newY > BoxSize / 2) {
                this.y = newY - BoxSize;
            } else {
                this.y = newY;
            }
        } else {
            // Bounce mode
            const radiusOffset = this.radius / 100;

            if (newX - radiusOffset < -BoxSize / 2 || newX + radiusOffset > BoxSize / 2) {
                this.vx *= -1; // Reverse velocity in X

                // prevents particles from getting stuck in boundaries
                this.x = Math.max(-BoxSize / 2 + radiusOffset, Math.min(BoxSize / 2 - radiusOffset, newX));
            }
            if (newY - radiusOffset < -BoxSize / 2 || newY + radiusOffset > BoxSize / 2) {
                this.vy *= -1; // Reverse velocity in Y

                // prevents particles from getting stuck in boundaries
                this.y = Math.max(-BoxSize / 2 + radiusOffset, Math.min(BoxSize / 2 - radiusOffset, newY));
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

// User input for boundary, number of particles, and velocity magnitude
let animationId = null;
document.getElementById('startSimulation').addEventListener('click', function () {
    // clean previous animation frame resource
    if (animationId) {
        cancelAnimationFrame(animationId);
        const canvasContainer = document.getElementById('canvasContainer');
        canvasContainer.innerHTML = '';
        animationId = null;
    }

    // validate the user input
    function validateInput(input) {
        if (isNaN(input) || input < 0.0) {
            return false;
        }
        return true;
    }

    const userBoxSize = parseFloat(document.getElementById('boxSize').value);
    const temperature = parseFloat(document.getElementById('temperature').value);
    const numParticles = parseInt(document.getElementById('numParticles').value);
    const boundaryType = parseInt(document.getElementById('boundaryType').value);				// interaction with walls
    const interactionType = parseInt(document.getElementById('interactionType').value); // interaction with particles

    const air_mass = 5.32 * Math.pow(10, -26);
    console.log(air_mass);
    const m = air_mass;
    const k = 1.38 * Math.pow(10, -23);;
    var userVelocity = Math.sqrt((3*k*temperature)/m)/1000;
    //var userVelocity = temperature / 10;
    var bounce = 0;
    

    // validate the input
    if (!validateInput(userBoxSize) || !validateInput(userBoxSize) || !validateInput(userVelocity)
        || !validateInput(numParticles)) {
        alert('Please enter the valid number for all fields.');
        return;
    }

    // Initialize with random starting positions
    function getRandomPosition(range) {
        return Math.random() * range * 2 - range;
    }

    // Generate random angle (in radians) between 0 and 2 * PI (360)
    function getRandomAngle() {
        return Math.random() * Math.PI * 2; // Random angle between 0 and 2π (360 degrees)
    }

    // Helper function to ensure particles don't overlap initially
    function isPositionValid(x, y, particles, minDistance, BoxSize) {
        if(Math.abs(x) > BoxSize / 2 - minDistance || Math.abs(y) > BoxSize / 2 - minDistance) {
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

    const particles = [];
    const particleRadius = 5;
    const minDistance = (particleRadius * 2) / 100; // Converted to simulation units

    // create numParticles number of particles in a random starting position moving in a random direction
    for (let i = 0; i < numParticles; i++) {
        const angle = getRandomAngle();  // Random direction (angle)

        // Calculate the x and y components of the velocity based on the random angle
        const vx = userVelocity * Math.cos(angle); // X velocity component
        const vy = userVelocity * Math.sin(angle); // Y velocity component

        // Try to find a valid non-overlapping position
        let x, y;
        let attempts = 0;
        const maxAttempts = 1000;

        do {
            x = getRandomPosition(userBoxSize / 2);
            y = getRandomPosition(userBoxSize / 2);
            attempts++;

            // Break out after too many attempts to prevent infinite loop
            if (attempts > maxAttempts) {
                console.warn('Could not find non-overlapping position after', maxAttempts, 'attempts');
                break;
            }
        } while (!isPositionValid(x, y, particles, minDistance, userBoxSize));

        const atom = new Atom(
            x, y, vx, vy, particleRadius, 1 // x, y, vx, vy, radius, mass
        );
        particles.push(atom);
    }

    // Canvas setup
    const canvas = document.createElement('canvas');
    canvas.id = "simulationCanvas";
    canvas.width = userBoxSize * 100;
    canvas.height = userBoxSize * 100;
    canvas.style.border = "1px solid black";

    const button = document.getElementById('startSimulation');
    const buttonTop = button.getBoundingClientRect().top;
    const buttonHeight = button.offsetHeight

    const canvasContainer = document.getElementById('canvasContainer');
    canvasContainer.style.position = "absolute";
    canvasContainer.style.top = `${buttonTop + buttonHeight}px`;
    canvasContainer.innerHTML = '';
    canvasContainer.style.width = `${canvas.width * 1.5}px`;
    canvasContainer.style.height = `${canvas.height * 1.5}px`;
    canvasContainer.appendChild(canvas);

    const ctx = canvas.getContext("2d");

    function draw() {
        ctx.fillStyle = 'rgba(255, 255, 255, 0.3)';
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        // Draw each atom
        particles.forEach(atom => {
            const drawX = (canvas.width / 2) + atom.x * 100;
            const drawY = (canvas.height / 2) - atom.y * 100;

            ctx.beginPath();
            ctx.arc(drawX, drawY, atom.radius, 0, Math.PI * 2);
            ctx.fillStyle = 'red';
            ctx.fill();
        });
    }

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

    // replace with RAF optimize the animate playing
    let lastTimestamp = 0;
    const maxTimeStep = 1.0 / 30.0;; // 30 FPS

    function animate(timestamp) {
        if (!lastTimestamp) lastTimestamp = timestamp;

        const deltaTime = timestamp - lastTimestamp;

        // limited 30 FPS
        const timeStep = Math.min(deltaTime / 1000, maxTimeStep);

        handleCollisions();

        // update particles position

        
        particles.forEach(atom => {
    atom.updatePosition(timeStep, userBoxSize, boundaryType);
    
    if (!boundaryType) {
        const radiusOffset = atom.radius / 100;
        atom.x = Math.max(-userBoxSize / 2 + radiusOffset, Math.min(userBoxSize / 2 - radiusOffset, atom.x));
        atom.y = Math.max(-userBoxSize / 2 + radiusOffset, Math.min(userBoxSize / 2 - radiusOffset, atom.y));
    }
});

        

        draw();
        lastTimestamp = timestamp;
        animationId = requestAnimationFrame(animate);
    }
    animationId = requestAnimationFrame(animate);
});
function collision_pressure_calc(mass, magnitude, bounce) {
    Bounce += 1;


}



