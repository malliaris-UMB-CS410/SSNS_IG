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
            const speed = rmsSpeed * Math.sqrt(-2 * Math.log(Math.random())); // Maxwell-Boltzmann distribution
            this.particles.push({
                x: Math.random(), // Random x position in [0, 1)
                y: Math.random(), // Random y position in [0, 1)
                vx: speed * Math.cos(theta), // Velocity x-component
                vy: speed * Math.sin(theta) // Velocity y-component
            });
        }
    }


    update(dt) {
        // Update particle positions based on velocity
        const canvasWidth = document.getElementById('plot_HM_IG').width;
        const canvasHeight = document.getElementById('plot_HM_IG').height;

        // Normalize radius to [0, 1] range
        const normalizedRadius = this.radius / Math.min(canvasWidth, canvasHeight);

        this.particles.forEach(particle => {
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

    // Default parameters
    let params = {
        T: parseFloat(inputT?.value || 0),
        N: parseInt(inputN?.value || 4)
    };

    const particleSimulate = new ParticleSimulate(params);
    const canvas = document.getElementById('plot_HM_IG');
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

    function animate() {
        particleSimulate.update(0.01);
        particleSimulate.draw(ctx);
        requestAnimationFrame(animate);
    }

    animate();
});
