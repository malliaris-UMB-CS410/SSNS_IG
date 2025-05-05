class SeededRandom {

    constructor(seed) {
        if (!seed) {
            seed = "ideal gas";
        }

        if (typeof seed === 'string') {
            this.state = this.stringToSeed(seed);
        } else if (typeof seed === 'number' && Number.isInteger(seed)) {
            this.state = seed;
        } 
    }

    stringToSeed(str) {
        let seed = 0;
        for (let i = 0; i < str.length; i++) {
            seed ^= str.charCodeAt(i) * (i + 1);
        }
        return seed;
    }

    next() {
        this.state ^= (this.state << 21);
        this.state ^= (this.state >>> 35);
        this.state ^= (this.state << 4);
        return this.state >>> 0;
    }

    getRandom(min, max) {
        if (min >= max) {
            throw new Error('Min should be less than max'); 
        }
        const rand = this.next();
        const range = max - min;
        return min + (rand / 0xFFFFFFFF) * range;  
    }

    getRandomInt(min, max) {
        if (min >= max) {
            throw new Error('Min should be less than max');
        }
        return Math.floor(this.getRandom(min, max + 1)); 
    }
}
