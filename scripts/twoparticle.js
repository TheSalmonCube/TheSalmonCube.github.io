document.addEventListener('DOMContentLoaded', function() {
    const canvas = document.getElementById('2p-plot');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');

    // Canvas size
    let width = 0, height = 0;
    let imageData = null;
    let data = null;
    
    // Constants - Simulation grid
    // In this context, X and Y refer to two different particles' 1D positions in configuration space.
    const CalcX = 50;
    const CalcY = 50;
    const CalcSize = CalcX * CalcY;
    const CalcDim = 20;
    const calc = window.createCalculator ? window.createCalculator(CalcSize, CalcDim) : null;
    const BatchT = 8;

    // Constants - Controlled parameters
    const M = 1;
    let controlX = 25;
    let controlY = 25;
    let controlSX = 3;
    let controlSY = 3;
    let controlPX = 0.75;
    let controlPY = 0.75;
    let controlVtype = 'free';
    let controlTStep = 0.3;

    // Constants - Rendering and animation
    let t = 0;
    const colorBoost = 7.0;

    // ensure imageData and data are initialized and canvas sized to aspect ratio
    function resizeAndInit() {
        // Make the canvas half the width of the <main> section and keep it square
        const mainEl = document.querySelector('main') || document.body;
        const mainRect = mainEl.getBoundingClientRect();
        const canvasSize = Math.max(200, Math.floor(mainRect.width * 0.5));
        if (canvas.width === canvasSize && canvas.height === canvasSize && imageData) {
            width = canvas.width;
            height = canvas.height;
            return;
        }
        canvas.width = canvasSize;
        canvas.height = canvasSize;
        width = canvas.width;
        height = canvas.height;
        // disable smoothing
        if (ctx.imageSmoothingEnabled !== undefined) ctx.imageSmoothingEnabled = false;
        imageData = ctx.createImageData(width, height);
        data = imageData.data;
    }

    // Recompute buffer 
    window.addEventListener('resize', function() {
        resizeAndInit();
        if (!running) updatePreview();
    });

    // Now we copy paste the 2D version which is equivalent to two particles in 1D each. Now nice is that?

    // Buffers for calculation
    const magnitude = new Float64Array(CalcSize);
    const phase = new Float64Array(CalcSize);

    // Render from current magnitude/phase buffers, interpolate calculated grid, apply overlay, and add to pixel buffer
    function renderFromGrid() {
        if (!imageData) {
            resizeAndInit();
            if (!imageData) return;
        }
         const maxXIndex = CalcX - 1;
         const maxYIndex = CalcY - 1;

        // 1. Determine interpolation indices and weights (2D lerp)
        for (let py = 0; py < height; py++) {
            const fy = (py / (height - 1)) * maxYIndex;
            const y00 = Math.floor(fy);
            const yfrac = fy - y00;
            const y00Offset = y00 * CalcX;
            const y10Offset = Math.min(y00 + 1, maxYIndex) * CalcX; // Prevents y10 from overflowing index
            for (let px = 0; px < width; px++) {
                const fx = (px / (width - 1)) * maxXIndex;
                const x00 = Math.floor(fx);
                const xfrac = fx - x00;

                const idx00 = y00Offset + x00;
                const idx10 = y00Offset + Math.min(x00 + 1, maxXIndex); // Prevents x10 from overflowing index
                const idx01 = y10Offset + x00;
                const idx11 = y10Offset + Math.min(x00 + 1, maxXIndex);

                const w00 = (1 - xfrac) * (1 - yfrac);
                const w10 = xfrac * (1 - yfrac);
                const w01 = (1 - xfrac) * yfrac;
                const w11 = xfrac * yfrac;

                // 2. Interpolate magnitude and phase (phase circular)
                let p00 = phase[idx00];
                let p10 = phase[idx10];
                let p01 = phase[idx01];
                let p11 = phase[idx11];

                const lerpmag = (w00 * magnitude[idx00] + w10 * magnitude[idx10] + w01 * magnitude[idx01] + w11 * magnitude[idx11]) * colorBoost;

                // Circular linear interpolation for phase by moving everything near p00
                if (p10 > p00 + Math.PI) p10 -= 2 * Math.PI;
                if (p01 > p00 + Math.PI) p01 -= 2 * Math.PI;
                if (p11 > p00 + Math.PI) p11 -= 2 * Math.PI;
                if (p10 < p00 - Math.PI) p10 += 2 * Math.PI;
                if (p01 < p00 - Math.PI) p01 += 2 * Math.PI;
                if (p11 < p00 - Math.PI) p11 += 2 * Math.PI;
                const lerpphase = w00 * p00 + w10 * p10 + w01 * p01 + w11 * p11;
                
                // 3. Convert to RGB and write to pixel buffer
                const i = (py * width + px) * 4;
                const s1 = Math.sin(lerpphase); // 3 phase sin: Phase determines hue
                const s2 = Math.sin(lerpphase + 2/3 * Math.PI);
                const s3 = Math.sin(lerpphase + 4/3 * Math.PI);
                data[i]     = Math.min(Math.floor((s1 + 1) * 0.5 * 255 * lerpmag),255); // Magnitude (and colorboost) determines brightness
                data[i + 1] = Math.min(Math.floor((s2 + 1) * 0.5 * 255 * lerpmag),255);
                data[i + 2] = Math.min(Math.floor((s3 + 1) * 0.5 * 255 * lerpmag),255);
                data[i + 3] = 255;
                
                // 4. Overlay potential if enabled
                if (overlayEnabled) {
                    const pot = Math.max(0, w00 * V[idx00] + w10 * V[idx10] + w01 * V[idx01] + w11 * V[idx11]);
                    if (pot > 1e-6) {
                        const alpha = Math.min(pot, 1);
                        data[i]   = Math.floor((1 - alpha) * data[i]   + alpha * 128);
                        data[i+1] = Math.floor((1 - alpha) * data[i+1] + alpha * 128);
                        data[i+2] = Math.floor((1 - alpha) * data[i+2] + alpha * 128);
                    }
                }
            }
        }
        ctx.putImageData(imageData, 0, 0);
    }

    // Current state
    let running = false;
    let state0 = null;
    let V = null;
    let Initializer = null;
    let H = null;
    let overlayEnabled = false;
    let controlsDirty = true;
    let lastState2D = null;

    function buildPotential(vtype) {
        const V = new Array(CalcSize).fill(0);
        // Interactive potentials V(X,Y) or V(x1,x2). If they are coupled, they represent interaction potentials between the two particles. 
        // Otherwise they are independent (seperable) and V = Vx(x)Vy(y)
        if (vtype === 'free') {
            // Free particle, V = 0 everywhere
        } else if (vtype === 'well') {
            // Independent well for each particle
            for (let xx = 0; xx < CalcX; xx++) {
                for (let yy = 0; yy < CalcY; yy++) {
                    const index = yy * CalcX + xx;
                    if ((Math.abs(xx-25) > 15) || (Math.abs(yy-25) > 15)) {V[index] = 1}  
                }
            }
        } else if (vtype === 'attraction') {
            // Coupled inverse distance attraction
            for (let xx = 0; xx < CalcX; xx++) {
                for (let yy = 0; yy < CalcY; yy++) {
                    const index = yy * CalcX + xx;
                    const r = Math.abs(xx - yy)
                    V[index] = -2/(r+1)
                }
            }
        } else if (vtype === 'repulsion') {
            // Coupled inverse distance attraction
            for (let xx = 0; xx < CalcX; xx++) {
                for (let yy = 0; yy < CalcY; yy++) {
                    const index = yy * CalcX + xx;
                    const r = Math.abs(xx - yy)
                    V[index] = 2/(r+1)
                }
            }
        }
        return V;
    }

    V = buildPotential(controlVtype); // init

    function buildInitialState(x0, y0, px, py, sx, sy) {
        const rex = new Float64Array(CalcX);
        const imx = new Float64Array(CalcX);
        const rey = new Float64Array(CalcY);
        const imy = new Float64Array(CalcY);
        
        const re = new Float64Array(CalcSize);
        const im = new Float64Array(CalcSize);
        for (let xx = 0; xx < CalcX; xx++) {
            const dx = xx - x0;
            const density = Math.exp(-0.5 * (dx * dx) / (sx * sx));
            rex[xx] = density * Math.cos(px * dx);
            imx[xx] = density * Math.sin(px * dx);
        }
        for (let yy = 0; yy < CalcY; yy++) {
            const dy = yy - y0;
            const density = Math.exp(-0.5 * (dy * dy) / (sy * sy));
            rey[yy] = density * Math.cos(py * dy);
            imy[yy] = density * Math.sin(py * dy);
        }
        
        for (let yy = 0; yy < CalcY; yy++) {
            for (let xx = 0; xx < CalcX; xx++) {
                const index = yy * CalcX + xx;
                // (a+ib)*(c+id) = (ac - bd) + i(ad + bc)
                const a = rex[xx];
                const b = imx[xx];
                const c = rey[yy];
                const d = imy[yy];
                re[index] = a * c - b * d;
                im[index] = a * d + b * c;
                 // Outer product of two 1D wavefunctions to create 2D wavefunction. Particles are unentangled at t=0,
                 // Demonstrated by seperability of wavefunction: Psi(x1,x2) = Psi1(x1) * Psi2(x2)
                 // But of course we collapse it down to 1D arrays.
             }
         }

         const state = [re, im];
         if (calc && calc.normalize) calc.normalize(state);
         return state;
     }

    state0 = buildInitialState(controlX, controlY, controlPX, controlPY, controlSX, controlSY); // init

    function readControls() {
        controlX = parseFloat(document.getElementById('inputX').value);
        controlY = parseFloat(document.getElementById('inputY').value);
        controlSX = parseFloat(document.getElementById('inputSX').value);
        controlSY = parseFloat(document.getElementById('inputSY').value);
        controlPX = parseFloat(document.getElementById('inputPX').value);
        controlPY = parseFloat(document.getElementById('inputPY').value);
        controlVtype = document.getElementById('inputV').value;
        controlTStep = parseFloat(document.getElementById('inputSpeed').value);
        document.getElementById('xVal').textContent = controlX.toFixed(1);
        document.getElementById('yVal').textContent = controlY.toFixed(1);
        document.getElementById('sxVal').textContent = controlSX.toFixed(1);
        document.getElementById('syVal').textContent = controlSY.toFixed(1);
        document.getElementById('pxVal').textContent = controlPX.toFixed(2);
        document.getElementById('pyVal').textContent = controlPY.toFixed(2);
    }   

    function updatePreview() {
        if (running) return;
        readControls();
        V = buildPotential(controlVtype);
        const preview = buildInitialState(controlX, controlY, controlPX, controlPY, controlSX, controlSY);

        // update image (Same as in render loop)

        for (let yy = 0; yy < CalcY; yy++) {
            for (let xx = 0; xx < CalcX; xx++) {
                const cidx = yy * CalcX + xx;
                const re = preview[0][cidx];
                const im = preview[1][cidx];
                magnitude[cidx] = Math.hypot(re, im);
                phase[cidx] = Math.atan2(im, re);
            }
        }

        renderFromGrid();

        // update side-panel charts from preview state immediately
        try {
            if (window.chart1Obj && window._chartXSamples) {
                const margX = new Float64Array(CalcX);
                const margY = new Float64Array(CalcY);
                for (let i = 0; i < CalcX; i++) margX[i] = 0;
                for (let i = 0; i < CalcY; i++) margY[i] = 0;
                for (let yy = 0; yy < CalcY; yy++) {
                    for (let xx = 0; xx < CalcX; xx++) {
                        const index = yy * CalcX + xx;
                        const re = preview[0][index];
                        const im = preview[1][index];
                        const prob = re * re + im * im;
                        margX[xx] += prob;
                        margY[yy] += prob;
                    }
                }
                const xs = window._chartXSamples;
                if (xs && xs.length > 0) {
                    const len = xs.length;
                    const arr1 = new Array(len);
                    const arr2 = new Array(len);
                    for (let i = 0; i < len; i++) {
                        const xi = Math.round(i * (CalcX - 1) / (len - 1));
                        const yi = Math.round(i * (CalcY - 1) / (len - 1));
                        arr1[i] = { x: xs[i], y: margX[xi] };
                        arr2[i] = { x: xs[i], y: margY[yi] };
                    }
                    window.chart1Obj.data.datasets[0].data = arr1;
                    window.chart2Obj.data.datasets[0].data = arr2;
                    try { window.chart1Obj.update('none'); } catch (e) {}
                    try { window.chart2Obj.update('none'); } catch (e) {}
                }
            }
        } catch (e) {}

        // store preview state for start
        state0 = preview;
    }    

    function Simulate(tval) {
        if (!calc || !Initializer) return;

        // COMPUTE STATE 
        const currentState = calc.evolveLanczos(Initializer, tval);
        // remember latest state so Stop can preserve it
        lastState2D = currentState;

        // RENDERING
        // Convert complex numbers from a+bi to polar for calc grid
        for (let yy = 0; yy < CalcY; yy++) {
            for (let xx = 0; xx < CalcX; xx++) {
                const cidx = yy * CalcX + xx;
                const re = currentState[0][cidx];
                const im = currentState[1][cidx];
                magnitude[cidx] = Math.hypot(re, im);
                phase[cidx] = Math.atan2(im, re);
            }
        }

        // SIDE PANEL
        if (window.chart1Obj && window._chartXSamples) {
            // compute marginals (probability density sums)
            const margX = new Float64Array(CalcX);
            const margY = new Float64Array(CalcY);
            for (let xx = 0; xx < CalcX; xx++) margX[xx] = 0;
            for (let yy = 0; yy < CalcY; yy++) margY[yy] = 0;
            for (let yy = 0; yy < CalcY; yy++) {
                for (let xx = 0; xx < CalcX; xx++) {
                    const index = yy * CalcX + xx;
                    const prob = magnitude[index] * magnitude[index];
                    margX[xx] += prob;
                    margY[yy] += prob;
                }
            }

            // Map marginals into {x,y} arrays matching chart x-samples
            const xs = window._chartXSamples;
            if (xs && xs.length > 0) {
                const len = xs.length;
                const arr1 = new Array(len);
                const arr2 = new Array(len);
                for (let i = 0; i < len; i++) {
                    // map sample index to nearest grid index
                    const xi = Math.round(i * (CalcX - 1) / (len - 1));
                    const yi = Math.round(i * (CalcY - 1) / (len - 1));
                    arr1[i] = { x: xs[i], y: margX[xi] };
                    arr2[i] = { x: xs[i], y: margY[yi] };
                }
                window.chart1Obj.data.datasets[0].data = arr1;
                window.chart2Obj.data.datasets[0].data = arr2;
                // update immediately without animation
                try { window.chart1Obj.update('none'); } catch (e) {}
                try { window.chart2Obj.update('none'); } catch (e) {}
            }
        }

        // STATUS ROW
        if (entanglementEl) {
            const yProjRe = new Float64Array(CalcY).fill(0);
            const yProjIm = new Float64Array(CalcY).fill(0);
            const xProjRe = new Float64Array(CalcX).fill(0);
            const xProjIm = new Float64Array(CalcX).fill(0);
            
            for (let xx = 0; xx < CalcX; xx++) {
                for (let yy = 0; yy < CalcY; yy++) {
                    const index = yy * CalcX + xx;
                    yProjRe[yy] += currentState[0][index];
                    yProjIm[yy] += currentState[1][index];
                    xProjRe[xx] += currentState[0][index];
                    xProjIm[xx] += currentState[1][index];
                }
            }

            const R1approx = [new Float64Array(CalcSize), new Float64Array(CalcSize)];
            
            for (let yy = 0; yy < CalcY; yy++) {
                for (let xx = 0; xx < CalcX; xx++) {
                    const index = yy * CalcX + xx;
                    // (a+ib)*(c+id) = (ac - bd) + i(ad + bc)
                    const a = xProjRe[xx];
                    const b = xProjIm[xx];
                    const c = yProjRe[yy];
                    const d = yProjIm[yy];
                    R1approx[0][index] = a * c - b * d;
                    R1approx[1][index] = a * d + b * c;
                }
            }
            calc.normalize(R1approx);

            const [diffRe, diffIm] = calc.dotProduct(currentState, R1approx);
            const entanglement = 1 - Math.sqrt(diffRe * diffRe + diffIm * diffIm);
            entanglementEl.textContent = entanglement.toFixed(2);
        }

        renderFromGrid();

        // REINITIALIZE
        if (t > BatchT) {
            t = 0;
            state0 = currentState; // Self evident
            Initializer = calc.InitializeLanczos(state0, H); // Same hamiltonian
        }
    }

    // START SIMULATION
    function startSimulation() {
        // Only rebuild H/Initializer when controls changed
        readControls();
        if (controlsDirty || !H || !Initializer) {
            V = buildPotential(controlVtype);
            H = calc.create2DHamiltonian(M, V, CalcX, CalcY);
            Initializer = calc.InitializeLanczos(state0, H);
            controlsDirty = false;
        }
        running = true;
        const btnStart = document.getElementById('btnStart');
        const btnStop = document.getElementById('btnStop');
        if (btnStart) btnStart.disabled = true;
        if (btnStop) btnStop.disabled = false;
    }

    function stopSimulation() {
        // stop but preserve H/Initializer so start resumes without rebuild
        running = false;
        // if we have a preserved recent state, make it the new state0 and reset time
        if (lastState2D) {
            state0 = lastState2D;
            // rebuild initializer so resuming starts NnOwfg
            if (H) Initializer = calc.InitializeLanczos(state0, H);
            t = 0;
        }
        const btnStart = document.getElementById('btnStart');
        const btnStop = document.getElementById('btnStop');
         if (btnStart) btnStart.disabled = false;
        if (btnStop) btnStop.disabled = true;
    }

    // Wire up UI controls
    (function attachUI(){
    const inputs = ['inputX','inputY','inputSX','inputSY','inputPX','inputPY','inputV','inputSpeed'];
        // include overlay checkbox afterwards
        inputs.forEach(id => {
            const el = document.getElementById(id);
            if (!el) return;
            const ev = (el.tagName === 'SELECT') ? 'change' : 'input';
            // mark controls dirty when user changes a control and update preview
            el.addEventListener(ev, function(){
                controlsDirty = true;
                updatePreview();
            });
        });
        // add overlay checkbox separately below
        // overlay checkbox also toggles overlayEnabled
        const overlayEl = document.getElementById('inputOverlay');
        if (overlayEl) {
            overlayEl.addEventListener('change', function(){
                overlayEnabled = !!overlayEl.checked;
                renderFromGrid();
            });
        }
        const btnStart = document.getElementById('btnStart');
        const btnStop = document.getElementById('btnStop');
        if (btnStart) btnStart.addEventListener('click', function(){
            // only update preview if controls were changed since last run
            if (controlsDirty) updatePreview();
            startSimulation();
        });
        if (btnStop) btnStop.addEventListener('click', stopSimulation);
        // initialize canvas, preview and button states
        resizeAndInit();
        // ensure preview uses correct buffer sizes
        updatePreview();
        if (btnStart) btnStart.disabled = false;
        if (btnStop) btnStop.disabled = true;
    })();

    // Cache status elements
    const statusTimeEl = document.getElementById('statusTime');
    const entanglementEl = document.getElementById('statusEntanglement');

    // Animation loop 

    let lastTs = performance.now();
    let acc = 0;
    const computeInterval = 100; // ms

    function loop(ts) {
        const dt = ts - lastTs;
        lastTs = ts;
        acc += dt;

        if (acc >= computeInterval) {
            if (running) {
                Simulate(t);
                t += controlTStep;
            }
            // update status readouts on each compute tick
            try {
                if (statusTimeEl) statusTimeEl.textContent = t.toFixed(2);
            } catch (e) {console.log(e)}

            acc = 0;
        }

        if (imageData) ctx.putImageData(imageData, 0, 0);
        requestAnimationFrame(loop);
    }

    requestAnimationFrame(loop);
});