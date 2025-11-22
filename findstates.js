const Calcsize = 512;
const Calcdim = 50;
const calc = window.createCalculator(Calcsize, Calcdim);

const M = 1;
const Free_Potential = new Float64Array(calc.size).fill(0.0);

const V = Free_Potential;

// Initialize Krylov vector with random numbers and normalize
const re = new Float64Array(calc.size);
const im = new Float64Array(calc.size);
for(let i = 0; i < calc.size; i++) {
    re[i] = Math.random();
    im[i] = Math.random();
}
const KrylovVector = [re, im] ;
calc.normalize(KrylovVector);

// k is KrylovVector
// Q is n x m orthonormal. Transforms vectors from main space to Krylov subspace, and Qt does reverse
// energies are the energies of corresponding vectors
// ritzVectors are eigenvectors of T in krylov space. Need to be returned back into main space
const [k, Q, energies, ritzVectors] = calc.InitializeLanczos(KrylovVector, M, V);

// Return Ritz vectors to main space
const wavefunctionsRe = calc.MM(Q[0], ritzVectors); // Re
const wavefunctionsIm = calc.MM(Q[1], ritzVectors); // Im

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
            label: 'real. E0 = '+energies[0].toFixed(5),
            data: Array.from(wavefunctionsRe[0])
        },{
            label: 'real. E49 = '+energies[49].toFixed(5),
            data: Array.from(wavefunctionsRe[49])
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
        animation: true,
        plugins: {
            subtitle: {
                display: true,
                text: 'T = 0'
            }
        },
    }
});


