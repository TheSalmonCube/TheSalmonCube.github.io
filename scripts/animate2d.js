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

    // Constants - Simulation grid
    const CalcX = 100;
    const CalcY = 50;
    const CalcSize = CalcX * CalcY;
    const CalcDim = 20;
    const calc = window.createCalculator ? window.createCalculator(CalcSize, CalcDim) : null;
    const BatchT = 8;

    // Constants - Controlled parameters
    const M = 1;
    let controlX = 25;
    let controlY = 25;
    let controlS = 3;
    let controlPX = 0.75;
    let controlPY = 0;
    let controlVtype = 'oneslit';
    let controlTStep = 0.3;

    // Constants - Rendering and animation
    let t = 0;
    const colorBoost = 5.0;

    // Buffers for calculation
    const magnitude = new Float64Array(CalcSize);
    const phase = new Float64Array(CalcSize);

    // Render from current magnitude/phase buffers, interpolate calculated grid, apply overlay, and add to pixel buffer
    function renderFromGrid() {
        if (!imageData) return;
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

    V = buildPotential(controlVtype); // init

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

    state0 = buildInitialState(controlX, controlY, controlPX, controlPY, controlS); // init

    function readControls() {
        controlX = parseFloat(document.getElementById('inputX').value);
        controlY = parseFloat(document.getElementById('inputY').value);
        controlS = parseFloat(document.getElementById('inputS').value);
        controlPX = parseFloat(document.getElementById('inputPX').value);
        controlPY = parseFloat(document.getElementById('inputPY').value);
        controlVtype = document.getElementById('inputV').value;
        controlTStep = parseFloat(document.getElementById('inputSpeed').value);
        document.getElementById('xVal').textContent = controlX.toFixed(1);
        document.getElementById('yVal').textContent = controlY.toFixed(1);
        document.getElementById('sVal').textContent = controlS.toFixed(1);
        document.getElementById('pxVal').textContent = controlPX.toFixed(2);
        document.getElementById('pyVal').textContent = controlPY.toFixed(2);
    }   

    function updatePreview() {
        if (running) return;
        readControls();
        V = buildPotential(controlVtype);
        const preview = buildInitialState(controlX, controlY, controlPX, controlPY, controlS);

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
            // rebuild initializer so resuming starts from this instant
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
        const inputs = ['inputX','inputY','inputS','inputPX','inputPY','inputV','inputSpeed'];
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
                 // do not mark controls dirty; do not trigger full preview render
                 // only redraw current frame so overlay toggles visually
                 // update image using current magnitude/phase buffers
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
        // initialize preview and button states
        updatePreview();
        if (btnStart) btnStart.disabled = false;
        if (btnStop) btnStop.disabled = true;
    })();

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
            if (running) {
                Simulate(t);
                t += controlTStep;
            }
            acc = 0;
        }

        if (imageData) ctx.putImageData(imageData, 0, 0);
        requestAnimationFrame(loop);
    }

    requestAnimationFrame(loop);
});