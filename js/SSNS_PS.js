class ParticleSimulate {
    constructor(params) {
        this.params = params;
        this.particles = [];
        this.radius = 5; // Define particle radius
        this.init();
    }

    init() {
        // Initialize particles based on parameters
        const { N } = this.params;
        this.updateParticles(N); // Use the updateParticles method to initialize
    }

    updateParticles(N) {
        if (!N || N <= 0) {
            console.error("Invalid number of particles:", N);
            return;
        }

        // Clear existing particles
        this.particles = [];

        // Calculate root mean square speed based on temperature
        const rmsSpeed = Math.sqrt(3 * this.params.T); // v_rms = sqrt(3 * T)

        // Initialize new particles
        for (let i = 0; i < N; i++) {
            // Generate random direction and scale by rmsSpeed
            const theta = Math.random() * 2 * Math.PI; // Random angle in [0, 2π)
            let speed = rmsSpeed * Math.sqrt(-2 * Math.log(Math.random())); // Maxwell-Boltzmann distribution

            this.particles.push({
                x: Math.random(), // Random x position in [0, 1)
                y: Math.random(), // Random y position in [0, 1)
                vx: speed * Math.cos(theta), // Velocity x-component
                vy: speed * Math.sin(theta), // Velocity y-component
                mass: 1 // Add a default mass property
            });
        }
    }

    update(dt) {
        const canvasWidth = document.getElementById('plot_HM_IG').width;
        const canvasHeight = document.getElementById('plot_HM_IG').height;
        const normalizedRadius = this.radius / Math.min(canvasWidth, canvasHeight);

        // Update particle positions based on velocity
        this.particles.forEach((particle, i) => {
            particle.x += particle.vx * dt;
            particle.y += particle.vy * dt;

            // Handle boundary conditions (bounce back with radius consideration)
            if (particle.x - normalizedRadius < 0 || particle.x + normalizedRadius > 1) {
                particle.vx = -particle.vx; // Reverse x velocity
                particle.x = Math.max(normalizedRadius, Math.min(1 - normalizedRadius, particle.x)); // Clamp position
            }
            if (particle.y - normalizedRadius < 0 || particle.y + normalizedRadius > 1) {
                particle.vy = -particle.vy; // Reverse y velocity
                particle.y = Math.max(normalizedRadius, Math.min(1 - normalizedRadius, particle.y)); // Clamp position
            }

            // Handle particle-particle collisions if enabled
            if (this.params.isCollision) {
                for (let j = i + 1; j < this.particles.length; j++) {
                    const other = this.particles[j];
                    const dx = particle.x - other.x;
                    const dy = particle.y - other.y;
                    const distance = Math.sqrt(dx * dx + dy * dy);

                    // Check if particles overlap
                    if (distance < 2 * normalizedRadius) {
                        const epsilon = 0.0001; // Minimum distance threshold
                        if (distance < epsilon) {
                            distance = epsilon; // Prevent division by zero or very small values
                        }

                        // Calculate collision response
                        const angle = Math.atan2(dy, dx);
                        const speed1 = Math.sqrt(particle.vx * particle.vx + particle.vy * particle.vy);
                        const speed2 = Math.sqrt(other.vx * other.vx + other.vy * other.vy);

                        const direction1 = Math.atan2(particle.vy, particle.vx);
                        const direction2 = Math.atan2(other.vy, other.vx);

                        const velocityX1 = speed1 * Math.cos(direction1 - angle);
                        const velocityY1 = speed1 * Math.sin(direction1 - angle);
                        const velocityX2 = speed2 * Math.cos(direction2 - angle);
                        const velocityY2 = speed2 * Math.sin(direction2 - angle);

                        const finalVelocityX1 = ((particle.mass - other.mass) * velocityX1 + 2 * other.mass * velocityX2) / (particle.mass + other.mass);
                        const finalVelocityX2 = ((other.mass - particle.mass) * velocityX2 + 2 * particle.mass * velocityX1) / (particle.mass + other.mass);

                        particle.vx = Math.cos(angle) * finalVelocityX1 + Math.cos(angle + Math.PI / 2) * velocityY1;
                        particle.vy = Math.sin(angle) * finalVelocityX1 + Math.sin(angle + Math.PI / 2) * velocityY1;
                        other.vx = Math.cos(angle) * finalVelocityX2 + Math.cos(angle + Math.PI / 2) * velocityY2;
                        other.vy = Math.sin(angle) * finalVelocityX2 + Math.sin(angle + Math.PI / 2) * velocityY2;

                        // Validate velocities
                        if (isNaN(particle.vx) || isNaN(particle.vy)) {
                            particle.vx = 0;
                            particle.vy = 0;
                        }
                        if (isNaN(other.vx) || isNaN(other.vy)) {
                            other.vx = 0;
                            other.vy = 0;
                        }

                        // Move particles apart to prevent overlap
                        const overlap = 2 * normalizedRadius - distance;
                        particle.x += overlap * dx / distance / 2;
                        particle.y += overlap * dy / distance / 2;
                        other.x -= overlap * dx / distance / 2;
                        other.y -= overlap * dy / distance / 2;

                        // Re-check boundary conditions after collision adjustments
                        if (particle.x - normalizedRadius < 0 || particle.x + normalizedRadius > 1) {
                            particle.vx = -particle.vx; // Reverse x velocity
                            particle.x = Math.max(normalizedRadius, Math.min(1 - normalizedRadius, particle.x)); // Clamp position
                        }
                        if (particle.y - normalizedRadius < 0 || particle.y + normalizedRadius > 1) {
                            particle.vy = -particle.vy; // Reverse y velocity
                            particle.y = Math.max(normalizedRadius, Math.min(1 - normalizedRadius, particle.y)); // Clamp position
                        }
                        if (other.x - normalizedRadius < 0 || other.x + normalizedRadius > 1) {
                            other.vx = -other.vx; // Reverse x velocity
                            other.x = Math.max(normalizedRadius, Math.min(1 - normalizedRadius, other.x)); // Clamp position
                        }
                        if (other.y - normalizedRadius < 0 || other.y + normalizedRadius > 1) {
                            other.vy = -other.vy; // Reverse y velocity
                            other.y = Math.max(normalizedRadius, Math.min(1 - normalizedRadius, other.y)); // Clamp position
                        }
                    }
                }
            }
        });
    }

    draw(ctx) {
        // Set canvas background to white
        ctx.fillStyle = "white";
        ctx.fillRect(0, 0, ctx.canvas.width, ctx.canvas.height);

        // Draw container border
        ctx.strokeStyle = "black"; // Border color
        ctx.lineWidth = 2; // Border thickness
        ctx.strokeRect(0, 0, ctx.canvas.width, ctx.canvas.height);

        // Draw particles in red
        ctx.fillStyle = "red";
        this.particles.forEach(particle => {
            const x = particle.x * ctx.canvas.width;
            const y = particle.y * ctx.canvas.height;
            ctx.beginPath();
            ctx.arc(x, y, this.radius, 0, 2 * Math.PI);
            ctx.fill();
        });
    }
}

// Ensure DOM elements are loaded before accessing them
document.addEventListener("DOMContentLoaded", () => {
    const inputN = document.getElementById('UI_P_SM_IG_N');
    const inputT = document.getElementById('UI_P_SM_IG_T');
    const inputCollision = document.getElementById('UI_P_SM_IG_COLLISION');

    // Default parameters
    let params = {
        T: parseFloat(inputT?.value || 0),
        N: parseInt(inputN?.value || 4),
        isCollision: document.getElementById('UI_P_SM_IG_COLLISION')?.checked || false
    };

    const particleSimulate = new ParticleSimulate(params);
    const canvas = document.getElementById('plot_HM_IG');
    canvas.width = 400;
    canvas.height = 400;
    const ctx = canvas.getContext('2d');

    // Listen for changes in T
    inputT.addEventListener("input", () => {
        const newT = parseFloat(inputT.value);
        if (!isNaN(newT)) {
            params.T = newT;
            particleSimulate.params.T = newT; // Update temperature in the simulation
            particleSimulate.updateParticles(params.N); // Reinitialize particles with new speed
        }
    });

    // Listen for changes in N
    inputN.addEventListener("input", () => {
        const newN = parseInt(inputN.value);
        if (!isNaN(newN)) {
            params.N = newN;
            particleSimulate.updateParticles(newN);
        }
    });

    // Listen for changes in Collision switch
    inputCollision.addEventListener("change", () => {
        params.isCollision = inputCollision.checked;
    });

    function animate() {
        particleSimulate.update(0.01);
        particleSimulate.draw(ctx);
        requestAnimationFrame(animate);
    }

    animate();
});
