class Atom {
    constructor(x, y, vx, vy, radius=5, m=1) {
        this.x = x;
        this.y = y;
        this.vx = vx;
        this.vy = vy;
        this.radius = radius;
        this.m = m;
    }

    updatePosition(timeStep, boxWidth, boxHeight, mode)  {
        let newX = this.x + this.vx * timeStep;
        let newY = this.y + this.vy * timeStep;

        if (mode) {
            // Wrap mode
            if (newX < -boxWidth / 2) {
                this.x = newX + boxWidth;
            } else if (newX > boxWidth / 2) {
                this.x = newX - boxWidth;
            } else {
                this.x = newX;
            }

            if (newY < -boxHeight / 2) {
                this.y = newY + boxHeight;
            } else if (newY > boxHeight / 2) {
                this.y = newY - boxHeight;
            } else {
                this.y = newY;
            }
        } else {
            // Bounce mode
            const radiusOffset = this.radius / 100;

            if (newX - radiusOffset < -boxWidth / 2 || newX + radiusOffset > boxWidth / 2) {
                this.vx *= -1; // Reverse velocity in X

                // prevents particles from getting stuck in boundaries
                this.x = Math.max(-boxWidth / 2 + radiusOffset, Math.min(boxWidth / 2 - radiusOffset, newX));
            }
            if (newY - radiusOffset < -boxHeight / 2 || newY + radiusOffset > boxHeight / 2) {
                this.vy *= -1; // Reverse velocity in Y

                // prevents particles from getting stuck in boundaries
                this.y = Math.max(-boxHeight / 2 + radiusOffset, Math.min(boxHeight / 2 - radiusOffset, newY));
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
        const m1 = this.m;
        const m2 = other.m;

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

        // Ensure particles are separated
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

// User input for boundary and initial conditions
document.getElementById('startSimulation').addEventListener('click', function() {
    const width = 5; 
    const height = 5;
    const maxSteps = parseInt(document.getElementById('maxSteps').value);
    const mode = parseInt(document.getElementById('mode').value);
    const N = parseInt(document.getElementById('particles').value);

    //let atom = new Atom(startX, startY, startVx, startVy, 1);
    const timeStep = 0.50;
    let stepCount = 0;

    // Create a canvas for visualization
    const canvas = document.createElement('canvas');
    canvas.id = "simulationCanvas";
    canvas.width = width * 100;
    canvas.height = height * 100;
    canvas.style.border = "1px solid black";

//         Math.random()*userBoxWidth/2, 
// Math.random()*userBoxHeight/2,


    const particles = [];

    for(let i=0; i< N; i++){
        const angle = getRandomAngle();  // Random direction (angle)

        // Calculate the x and y components of the velocity based on the random angle
        const vx = initialVelocity() * Math.cos(angle)*0.01; // X velocity component
        const vy = initialVelocity() * Math.sin(angle)*0.01; // Y velocity component

        const atom = new Atom(
            getRandomPositionX(),
            getRandomPositionY(),
            vx,  // Random X velocity component
            vy   // Random Y velocity component
        );
        particles.push(atom);
    }

    // get simulation button object for position top and height
    const canvasContainer = document.getElementById('canvasContainer');
    const distributionContainer = document.getElementById('atomContainer');
    const anotherCanvasDiv = document.createElement('div');
    anotherCanvasDiv.id = 'distributionContainer'

    const atomCanvas = document.createElement('canvas');
    atomCanvas.id = 'simulationCanvas';
    atomCanvas.width = width*100; // Set your desired width
    atomCanvas.height = height*100; // Set your desired height
    atomCanvas.style.border = "1px solid black";
    distributionContainer.appendChild(atomCanvas);

    const distributionCanvas = document.createElement('canvas');
    distributionCanvas.id = 'distributionCanvas';
    distributionCanvas.width = width*200; // Set your desired width
    distributionCanvas.height = height*100; // Set your desired height
    distributionCanvas.style.border = "1px solid blue";
    canvasContainer.appendChild(distributionCanvas); // Append to the main container

    const button = document.getElementById('startSimulation');
    const buttonTop = button.getBoundingClientRect().top;
    const buttonHeight = button.offsetHeight

    // Atoms container

    // Distribution Container


    const ctx = atomCanvas.getContext("2d");
    const dtr = distributionCanvas.getContext("2d");

    function initialVelocity(){
        const m = 6.44*Math.pow(10, -24); // Considering Helium
        const k = 1.381*Math.pow(10, -23);
        const T = parseInt(document.getElementById('temperature').value);
        const pi = Math.PI;
        const vel = (4/pi)*Math.pow((k*T/m), 0.5)*Math.random()*10;
        return vel;
    }

    function getRandomPositionX() {
        return Math.random() * width - (width / 2);// random between -width/2 and +width/2
    }
    function getRandomPositionY() {
        return Math.random() * height - (height / 2);
    }
    
    // Generate random angle (in radians) between 0 and 2 * PI (360)
    function getRandomAngle() {
        return Math.random() * Math.PI * 2; // Random angle between 0 and 2π (360 degrees)
    }

   // let velocities = []; // Declare vel outside the drawDistribution function
    //let freq = []; // Declare freq outside the drawDistribution function
    
    function calculateVelocityDistribution(particles) {
        velocities = particles.map(p => Math.sqrt(p.vx * p.vx + p.vy * p.vy));
        velocities.sort((a, b) => a - b);
    
        let frequencyMap = new Map();
        for (let velocity of velocities) {
            let roundedVelocity = Math.round(velocity * 100) / 100; // Round velocity to two decimal places
            frequencyMap.set(roundedVelocity, (frequencyMap.get(roundedVelocity) || 0) + 1);
        }
    
        // Convert map to arrays for drawing
        velocities = [...frequencyMap.keys()];
        freq = [...frequencyMap.values()];
    }
    
    function drawDistribution() {
        dtr.clearRect(0, 0, distributionCanvas.width, distributionCanvas.height);
        if (!velocities.length || !freq.length) return;
    
        const maxFrequency = Math.max(...freq);
        const barWidth = (distributionCanvas.width / velocities.length) * 0.5;
        let x = (distributionCanvas.width * 0.1) / 24 // Starting position with margin
    
        for (let i = 0; i < velocities.length; i++) {
            const barHeight = (freq[i] / maxFrequency) * (distributionCanvas.height * 0.8);
            dtr.fillStyle = "blue";
            dtr.fillRect(x, distributionCanvas.height - barHeight, barWidth, barHeight);
            x += barWidth + 5; // Move to next bar position
        }
    }
    

    function draw() {
        ctx.fillStyle = 'rgba(255, 255, 255, 0.9)';
        ctx.fillRect(0, 0, atomCanvas.width, atomCanvas.height);
        
         particles.forEach(atom => { // Iterate through each atom in the array
         // Transform coordinates to fit canvas where (0,0) is at center and scale by 100
            let drawX = (atomCanvas.width / 2) + (atom.x * 100);
            let drawY = (atomCanvas.height / 2) - (atom.y * 100);
        
      // Ensure atoms are within the visible area (optional debug step)
            //console.log(`Drawing atom at: (${drawX.toFixed(2)}, ${drawY.toFixed(2)}) within canvas: (${canvas.width}, ${canvas.height})`);
        
         // Draw atom
            ctx.beginPath();
            ctx.arc(drawX, drawY, atom.radius, 0, Math.PI * 2);
            ctx.fillStyle = "red";
            ctx.fill();
        });
    }

    function handleCollisions() {
        // Check each pair of particles exactly once
        for (let i = 0; i < particles.length-1; i++) {
            for (let j = i + 1; j < particles.length; j++) {
                const particle1 = particles[i];
                const particle2 = particles[j];

                if (particle1.checkCollision(particle2)) {
                    particle1.resolveCollision(particle2);
                }
            }
        }
    }

    const interval = setInterval(() => {
        if (stepCount >= maxSteps) {
            clearInterval(interval);
            console.log("Simulation stopped after " + maxSteps + " time steps.");
            return;
        }

        handleCollisions();
        for (let i = 0; i < N; i++) {
            particles[i].updatePosition(timeStep, width, height, mode);
            const radiusOffset = particles[i].radius / 100;
            particles[i].x = Math.max(-width / 2 + radiusOffset, Math.min(width / 2 - radiusOffset, particles[i].x));
            particles[i].y = Math.max(-height / 2 + radiusOffset, Math.min(height / 2 - radiusOffset, particles[i].y));
            //console.log(`Atom Position: (${particles[i].x.toFixed(2)}, ${particles[i].y.toFixed(2)})`);
        }
        draw();
        calculateVelocityDistribution(particles); // Calculate the velocity distribution
        drawDistribution(); // Draw the distribution graph
        stepCount++;
    }, 100);

})
