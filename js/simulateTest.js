class Atom {
    constructor(x, y, vx, vy, m) {
        this.x = x;
        this.y = y;
        this.vx = vx;
        this.vy = vy;
        this.m = m;
    }

    updatePosition(timeStep, boxWidth, boxHeight, mode) {
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
            if (newX < -boxWidth / 2 || newX > boxWidth / 2) {
                this.vx *= -1;
            }
            if (newY < -boxHeight / 2 || newY > boxHeight / 2) {
                this.vy *= -1;
            }
            this.x = newX;
            this.y = newY;
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
    const vel = initialVelocity();

    for(let i=0; i< N; i++){
        const angle = getRandomAngle();  // Random direction (angle)

        // Calculate the x and y components of the velocity based on the random angle
        const vx = vel * Math.cos(angle)*0.01; // X velocity component
        const vy = vel * Math.sin(angle)*0.01; // Y velocity component

        const atom = new Atom(
            getRandomPositionX(),
            getRandomPositionY(),
            vx,  // Random X velocity component
            vy   // Random Y velocity component
        );
        particles.push(atom);
    }

    // get simulation button object for position top and height
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


    function initialVelocity(){
        const m = 6.44*Math.pow(10, -24); // Considering Helium
        const k = 1.381*Math.pow(10, -23);
        const T = parseInt(document.getElementById('temperature').value);
        const pi = Math.PI;
        const vel = (4/pi)*Math.pow((k*T/m), 0.5);
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
    function draw() {
         ctx.clearRect(0, 0, canvas.width, canvas.height);
        
         particles.forEach(atom => { // Iterate through each atom in the array
         // Transform coordinates to fit canvas where (0,0) is at center and scale by 100
            let drawX = (canvas.width / 2) + (atom.x * 100);
            let drawY = (canvas.height / 2) - (atom.y * 100);
        
      // Ensure atoms are within the visible area (optional debug step)
            console.log(`Drawing atom at: (${drawX.toFixed(2)}, ${drawY.toFixed(2)}) within canvas: (${canvas.width}, ${canvas.height})`);
        
         // Draw atom
         ctx.beginPath();
            ctx.arc(drawX, drawY, 4, 0, Math.PI * 2);
            ctx.fillStyle = "red";
            ctx.fill();
            ctx.closePath();
        });
    }

    const interval = setInterval(() => {
        if (stepCount >= maxSteps) {
            clearInterval(interval);
            console.log("Simulation stopped after " + maxSteps + " time steps.");
            return;
        }
        //arrayAtoms.forEach(atom => {
        //atom.updatePosition(timeStep, userBoxWidth, userBoxHeight, mode);
        //});
        for (let i = 0; i < N; i++) {
            particles[i].updatePosition(timeStep, width, height);
            console.log(`Atom Position: (${particles[i].x.toFixed(2)}, ${particles[i].y.toFixed(2)})`);
        }
        draw();
        stepCount++;
    }, 100);
})
