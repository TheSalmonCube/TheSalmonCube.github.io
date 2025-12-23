document.addEventListener('DOMContentLoaded', function() {
    const canvas = document.getElementById('2d-plot');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');

    // Canvas dynamic size (internal buffer) will be set in resizeCanvas()
    let width = 0, height = 0;
    let imageData = null;
    let data = null;

    function resizeCanvas() {
        const dpr = window.devicePixelRatio || 1;
        const rect = canvas.getBoundingClientRect();
        const w = Math.max(1, Math.floor(rect.width * dpr));
        const h = Math.max(1, Math.floor(rect.height * dpr));
        if (w === width && h === height) return;
        width = w; height = h;
        canvas.width = width;
        canvas.height = height;
        canvas.style.width = rect.width + 'px';
        canvas.style.height = rect.height + 'px';
        ctx.setTransform(1,0,0,1,0,0);
        imageData = ctx.createImageData(width, height);
        data = imageData.data;
    }

    window.addEventListener('resize', resizeCanvas);
    window.addEventListener('orientationchange', resizeCanvas);
    resizeCanvas();

    // Calculation grid 
    const CalcX = 100;
    const CalcY = 50;
    const CalcSize = CalcX * CalcY;
    const CalcDim = 20;
    const calc = window.createCalculator ? window.createCalculator(CalcSize, CalcDim) : null;
    const BatchT = 8;

    const M = 1;
    let controlX = 25;
    let controlY = 25;
    let controlS = 3;
    let controlPX = 1;
    let controlPY = 0;
    let controlVtype = 'oneslit';
    let controlTStep = 0.3;

    function buildPotential(vtype) {
        const V = new Array(CalcSize).fill(0);

        if (vtype === 'highwall') {
            for (let yy = 0; yy < CalcY; yy++) {
                for (let xx = 0; xx < CalcX; xx++) {
                    const index = yy * CalcX + xx;
                    if (xx > CalcX * 0.5 - 3 && xx < CalcX * 0.5 + 3) {
                        V[index] = 1;
                    }
                }
            }
        } else if (vtype === 'lowwall') {
            for (let yy = 0; yy < CalcY; yy++) {
                for (let xx = 0; xx < CalcX; xx++) {
                    const index = yy * CalcX + xx;
                    if (xx > CalcX * 0.5 - 3 && xx < CalcX * 0.5 + 3) {
                        V[index] = 0.3;
                    }
                }
            }
        } else if (vtype === 'oneslit') {
            for (let yy = 0; yy < CalcY; yy++) {
                for (let xx = 0; xx < CalcX; xx++) {
                    const index = yy * CalcX + xx;
                    if ((xx > CalcX * 0.5 - 3 && xx < CalcX * 0.5 + 3)
                    && (yy < CalcY * 0.5 - 2 || yy > CalcY * 0.5 + 2)) {    
                        V[index] = 1;
                    }
                }
            }
        } else if (vtype === 'twoslit') {
            for (let yy = 0; yy < CalcY; yy++) {
                for (let xx = 0; xx < CalcX; xx++) {
                    const index = yy * CalcX + xx;
                    if ((xx > CalcX * 0.5 - 3 && xx < CalcX * 0.5 + 3)
                    && (yy < CalcY * 0.5 - 6 || (yy > CalcY * 0.5 - 2 && yy < CalcY * 0.5 + 2) || yy > CalcY * 0.5 + 6)) {    
                        V[index] = 1;
                    }
                }
            }
        }
        return V;
    }

    function buildInitialState(x0, y0, px, py, standard_dev) {
        const re = new Float64Array(CalcSize);
        const im = new Float64Array(CalcSize);
        for (let yy = 0; yy < CalcY; yy++) {
            for (let xx = 0; xx < CalcX; xx++) {
                const index = yy * CalcX + xx;
                const dx = xx - x0;
                const dy = yy - y0;
                const density = Math.exp(-0.5 * (dx * dx + dy * dy) / (standard_dev * standard_dev));
                re[index] = density * Math.cos(px * dx + py * dy);
                im[index] = density * Math.sin(px * dx + py * dy);
            }
        }
        const state = [re, im];
        if (calc && calc.normalize) calc.normalize(state);
        return state;
    }

    // Initialize Hamiltonian and Lanczos
    let state0 = buildInitialState(controlX, controlY, controlPX, controlPY, controlS);
    const V = buildPotential(controlVtype);
    const H = calc.create2DHamiltonian(M, V, CalcX, CalcY);
    let Initializer = calc.InitializeLanczos(state0, H);

    // Buffers for calculation
    const magnitude = new Float64Array(CalcSize);
    const phase = new Float64Array(CalcSize);

    let t = 0;
    const colorBoost = 5.0;

    function Simulate(tval) {
        if (!calc || !Initializer) return;

        // COMPUTE STATE 
        const currentState = calc.evolveLanczos(Initializer, tval);

        // RENDERING
        // Convert cartesian to polar for calc grid
        for (let yy = 0; yy < CalcY; yy++) {
            for (let xx = 0; xx < CalcX; xx++) {
                const cidx = yy * CalcX + xx;
                const re = currentState[0][cidx];
                const im = currentState[1][cidx];
                magnitude[cidx] = Math.hypot(re, im);
                phase[cidx] = Math.atan2(im, re);
            }
        }

        // Interpolate calc grid values into pixel buffer
        const maxXIndex = CalcX - 1;
        const maxYIndex = CalcY - 1;
        for (let py = 0; py < height; py++) {
            // Find the closest calc grid points and interpolation weights based on distance to that point
            const fy = (py / (height - 1)) * maxYIndex * (1 - 1e-8);
            const y00 = Math.floor(fy);
            const yfrac = fy - y00;
            const y00Offset = y00 * CalcX;
            const y10Offset = Math.min(y00 + 1, maxYIndex) * CalcX;
            for (let px = 0; px < width; px++) {
                const fx = (px / (width - 1)) * maxXIndex * (1 - 1e-8);
                const x00 = Math.floor(fx);
                const xfrac = fx - x00;

                const idx00 = y00Offset + x00;
                const idx10 = y00Offset + Math.min(x00 + 1, maxXIndex);
                const idx01 = y10Offset + x00;
                const idx11 = y10Offset + Math.min(x00 + 1, maxXIndex);

                const w00 = (1 - xfrac) * (1 - yfrac);
                const w10 = xfrac * (1 - yfrac);
                const w01 = (1 - xfrac) * yfrac;
                const w11 = xfrac * yfrac;

                // 00      | 01
                //         y.
                //         fr
                //         |
                // -xfrac--x
                // 10        11

                const lerpmag = (w00 * magnitude[idx00] + w10 * magnitude[idx10] + w01 * magnitude[idx01] + w11 * magnitude[idx11]) * colorBoost;

                // Circlular linear interpolation for phase. Based around phase at idx00, if things are more than pi away, wrap them round
                if (phase[idx10] > phase[idx00] + Math.PI) {phase[idx10] -= 2 * Math.PI};
                if (phase[idx01] > phase[idx00] + Math.PI) {phase[idx01] -= 2 * Math.PI};
                if (phase[idx11] > phase[idx00] + Math.PI) {phase[idx11] -= 2 * Math.PI};
                if (phase[idx10] < phase[idx00] - Math.PI) {phase[idx10] += 2 * Math.PI};
                if (phase[idx01] < phase[idx00] - Math.PI) {phase[idx01] += 2 * Math.PI};
                if (phase[idx11] < phase[idx00] - Math.PI) {phase[idx11] += 2 * Math.PI}
                const lerpphase = w00 * phase[idx00] + w10 * phase[idx10] + w01 * phase[idx01] + w11 * phase[idx11]; 
                

                // Convert to RGB
                const i = (py * width + px) * 4; // Index in data array, same as normal index
                const s1 = Math.sin(lerpphase); // Three phase sine waves for hue based on phase
                const s2 = Math.sin(lerpphase + 2/3 * Math.PI);
                const s3 = Math.sin(lerpphase + 4/3 * Math.PI);
                data[i] = Math.min(Math.floor((s1 + 1) * 0.5 * 255 * lerpmag),255); // Multiply by magnitude for brightness (and colorboost) and add to R
                data[i + 1] = Math.min(Math.floor((s2 + 1) * 0.5 * 255 * lerpmag),255); // G
                data[i + 2] = Math.min(Math.floor((s3 + 1) * 0.5 * 255 * lerpmag),255); // and B
                data[i + 3] = 255; // Opacity
            }
        }

        // REINITIALIZE
        if (t > BatchT) {
            t = 0;
            state0 = currentState; // Self evident
            Initializer = calc.InitializeLanczos(state0, H); // Same hamiltonian
        }
    }

    // Animation loop 
    let lastTs = performance.now();
    let acc = 0;
    const computeInterval = 100; // ms

    function loop(ts) {
        resizeCanvas(); // ensures imageData exists and matches displayed size
        const dt = ts - lastTs;
        lastTs = ts;
        acc += dt;

        if (acc >= computeInterval) {
            Simulate(t);
            t += controlTStep;
            acc = 0;
        }

        if (imageData) ctx.putImageData(imageData, 0, 0);
        requestAnimationFrame(loop);
    }

    requestAnimationFrame(loop);
});