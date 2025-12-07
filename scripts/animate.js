const ctx = document.getElementById('myChart');

// create a calculator instance bound to a specific size
const CalcSize = 512; 
const CalcDim = 20;
const calc = window.createCalculator(CalcSize, CalcDim);

const Xaxis = [];
for(let i = 0; i < calc.size; i++) {
    Xaxis.push(i);
}

// Inputs later
const M = 1;
let controlX = 256;
let controlS = 5;
let controlP = 1;
let controlVtype = 'free';
let controlTStep = 0.2;

const plot = new Chart(ctx, {
    type: 'line',
    data: {
        labels: Xaxis,
        datasets: [{
            label: 'Real',
            data: new Array(calc.size).fill(0)
        },{
            label: 'Imaginary',
            data: new Array(calc.size).fill(0)
        },{
            label: 'Potential',
            data: new Array(calc.size).fill(0)
        }]
    },
    options: {
        scales: {
            y: {
            min: -0.6,
            max: 0.6
            }
        },
        animation: false,
        plugins: {
            subtitle: {
                display: true,
                text: 'T = 0'
            }
        },
    }
});

const batchSteps = 20;
let step = 0;
let batch = 1;
let totalTime = 0;

let running = false;

let state = null;
let V = null;
let mainInitializer = null; 

function buildPotential(name) {
    const n = calc.size;
    const V = new Float64Array(n);
    if (name === 'free') {
        V.fill(0);
    } else if (name === 'low') {
        V.fill(0);
        const a = Math.floor(n * 0.52734375); // ~270/512
        const b = Math.floor(n * 0.546875); // ~280/512
        V.fill(0.25, a, b);
    } else if (name === 'high') {
        V.fill(0);
        const a = Math.floor(n * 0.52734375);
        const b = Math.floor(n * 0.546875);
        V.fill(0.5, a, b);
    } else if (name === 'well') {
        V.fill(0.0);
        const a = Math.floor(n * 0.3046875); // ~156/512
        const b = Math.floor(n * 0.6953125); // ~356/512
        V.fill(1, 0, a);
        V.fill(1, b, n);
    }
    return V;
}

function readControls() {
    controlX = Math.floor(document.getElementById('inputX').value);
    controlS = parseFloat(document.getElementById('inputS').value);
    controlP = parseFloat(document.getElementById('inputP').value);
    controlVtype = document.getElementById('inputV').value;
    controlTStep = parseFloat(document.getElementById('inputSpeed').value);
    document.getElementById('xVal').textContent = controlX;
    document.getElementById('sVal').textContent = controlS;
    document.getElementById('pVal').textContent = controlP.toFixed(2);
}

function updatePreview() {
    if (running) return;
    readControls();
    V = buildPotential(controlVtype);
    const preview = calc.gaussian(controlX, controlP, controlS);
    // update chart datasets
    plot.data.datasets[0].data = Array.from(preview[0]);
    plot.data.datasets[1].data = Array.from(preview[1]);
    plot.data.datasets[2].data = Array.from(V);
    plot.options.plugins.subtitle.text = 'Preview (t = 0)';
    plot.update();
    // store preview state and potential for start
    state = preview;
}

function startSimulation() { 
    if (running) return;
    if (!state) updatePreview();
    // disable controls
    ['inputX','inputS','inputP','inputV','inputSpeed'].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.disabled = true;
    });
    document.getElementById('btnStart').disabled = true;
    document.getElementById('btnStop').disabled = false;

    // initialize Lanczos for first batch with current state, potential
    const potential = V || buildPotential(controlVtype);
    mainInitializer = calc.InitializeLanczos(state, M, potential); 

    step = 0;
    totalTime = 0;
    running = true;

    const dt = controlTStep;
    const updateInterval = 50; //ms

    interval = setInterval(() => {
        // Propagate by one tStep each frame
        const currentState = calc.evolveLanczos(mainInitializer, step * dt);

        // Render
        plot.data.datasets[0].data = Array.from(currentState[0]);
        plot.data.datasets[1].data = Array.from(currentState[1]);
        plot.options.plugins.subtitle.text = 'T = ' + (totalTime + dt).toFixed(3);
        plot.update();

        // Every batchSteps replace with the nextbatch
        if (step >= batchSteps) {
            mainInitializer = calc.InitializeLanczos(currentState, M, potential);
            step = 0;
        } else {
            step ++;
        }
        totalTime += dt;
    }, updateInterval);
}

function stopSimulation() {
    if (!running) return;
    clearInterval(interval);
    interval = null;
    running = false;
    ['inputX','inputS','inputP','inputV','inputSpeed'].forEach(id => {
        const el = document.getElementById(id);
        if (el) el.disabled = false;
    });
    document.getElementById('btnStart').disabled = false;
    document.getElementById('btnStop').disabled = true;
    document.getElementById('btnStart').innerHTML = "Restart";
}

function setupControls() {
    const inputIds = ['inputX','inputS','inputP','inputV','inputSpeed'];
    inputIds.forEach(id => {
        const el = document.getElementById(id);
        if (!el) return;
        // Preview update in real time for sliders/select (input event)
        el.addEventListener('input', updatePreview);
    });

    document.getElementById('btnStart').addEventListener('click', startSimulation);
    document.getElementById('btnStop').addEventListener('click', stopSimulation);

    // initialize preview
    updatePreview();
}

window.addEventListener('load', () => {
    // call your setup if present
    if (typeof setupControls === 'function') {
        setupControls();
        updatePreview?.();
        return;
    }

    // fallback wiring if setupControls is missing
    document.getElementById('btnStart')?.addEventListener('click', startSimulation);
    document.getElementById('btnStop')?.addEventListener('click', stopSimulation);

    // try to preview immediately
    try { updatePreview(); } catch (e) { console.warn('preview failed:', e); }
});

window.addEventListener('resize', function() {
    if (plot) {
        plot.resize();
    }
});

