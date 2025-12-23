const Calcsize = 1000;
const Calcdim = 100;
const calc = window.createCalculator(Calcsize, Calcdim);

// Define potential
const M = 1;
const Free_Potential = new Float64Array(calc.size).fill(0.0);
const Square_Potential = new Float64Array(calc.size);
for (let i = 0; i < calc.size; i++) {
    Square_Potential[i] = 0.00001*(i-calc.size/2)**2;
}
const Well_potential = new Float64Array(calc.size).fill(50.0);
for (let i = 0; i < calc.size; i++) {
    if (i > 3*calc.size/8 && i < 5*calc.size/8) {
        Well_potential[i] = 0.0;
    }
}
const Mexican_potential = new Float64Array(calc.size);
for (let i = 0; i < calc.size; i++) {
    const x = (i - calc.size/2)/50;
    Mexican_potential[i] = 0.01*(x**2 - 1)**2;
}

// Select potential to use
const V = Square_Potential; 

// Initialize random Krylov vector (guess) with random numbers and normalize
let re = new Float64Array(calc.size);
let im = new Float64Array(calc.size);
for(let i = 0; i < calc.size; i++) {
    re[i] = Math.random() -0.5;
    im[i] = Math.random() -0.5;
}
let KrylovVector = [re, im];
calc.normalize(KrylovVector);


const iterations = 10;
const energyfocus = 0; // Zero is good
for (let iteration = 0; iteration < iterations; iteration++) {
    // k is KrylovVector
    // Q is n x m orthonormal. Transforms vectors from main space to Krylov subspace, and Qt does reverse
    // energies are the energies of corresponding vectors
    // ritzVectors are eigenvectors of T in krylov space. Need to be returned back into main space
    
    let hamiltonian = calc.createHamiltonian(M, V);
    var [k, Q, energies, ritzVectors] = calc.InitializeLanczos(KrylovVector, hamiltonian);
    
    let wavefunctionGroundRe = calc.MM(Q[0], [ritzVectors[energyfocus]])[0]; // Re
    let wavefunctionGroundIm = calc.MM(Q[1], [ritzVectors[energyfocus]])[0]; // Im

    KrylovVector = [wavefunctionGroundRe, wavefunctionGroundIm];
    calc.normalize(KrylovVector);
}

const check = false;// After iterations, check orthogonality of Q basis returned by last InitializeLanczos
if (Q && Q[0] && check) {
    const qcheck = calc.checkOrthogonality(Q[0], Q[1], 1e-8);
    console.log(`Q orthogonality check: maxOverlap=${qcheck.maxOverlap}, pair=${qcheck.pair}, isOrthogonal=${qcheck.isOrthogonal}`);
}

// Return Ritz vectors to main space
const wavefunctionsRe = calc.MM(Q[0], ritzVectors); // Re
const wavefunctionsIm = calc.MM(Q[1], ritzVectors); // Im
const magnitudes = [];
for (let i = 0; i < wavefunctionsRe.length; i++) {
    let mag = new Float64Array(calc.size);
    for (let j = 0; j < calc.size; j++) {
        mag[j] = Math.sqrt(wavefunctionsRe[i][j]**2 + wavefunctionsIm[i][j]**2);
    }
    magnitudes.push(mag);
}


const ctx = document.getElementById('myChart');
const Xaxis = [];
for(let i = 0; i < calc.size; i++) {
    Xaxis.push(i);
}

const plot = new Chart(ctx, {
    type: 'line',
    data: {
        labels: Xaxis,
        datasets: [{
            label: 'E0 = '+energies[0].toFixed(5),
            data: Array.from(magnitudes[0])
        },{
            label: 'E1 = '+energies[1].toFixed(5),
            data: Array.from(magnitudes[1])
        },{
            label: 'E2 = '+energies[2].toFixed(5),
            data: Array.from(magnitudes[2])
        },{
            label: 'E3 = '+energies[3].toFixed(5),
            data: Array.from(magnitudes[3])
        },{
            label: 'Potential',
            data: Array.from(V)
        }]
    },
    options: {
        scales: {
            y: {
            min: -0.6,
            max: 0.6
            }
        },
    }
});


