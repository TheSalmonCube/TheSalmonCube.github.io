const size = 512;
const d = 20; // Krylov Dim

function norm(state) {
    let totalSquared = 0;
    for(let i = 0; i < size; i++) {
        totalSquared += state[0][i] * state[0][i] + state[1][i] * state[1][i]
    }
    return Math.sqrt(totalSquared);
}

function normalize(state) {
    // in place
    let totalSquared = 0;
    for(let i = 0; i < size; i++) {
        totalSquared += state[0][i] * state[0][i] + state[1][i] * state[1][i]
    }
    const norm = Math.sqrt(totalSquared);
    for(let i = 0; i < size; i++) {
        state[0][i] /= norm;
        state[1][i] /= norm;
    }
}

function dotProduct(stateA, stateB) {
    let realDot = 0.0;
    let imagDot = 0.0;
    for (let i = 0; i < size; i++) {
        realDot += stateA[0][i] * stateB[0][i] + stateA[1][i] * stateB[1][i];
        imagDot += stateA[0][i] * stateB[1][i] - stateA[1][i] * stateB[0][i];
    }
    return [realDot, imagDot];
}

function oneMGS(state, statesRe, statesIm, normal = 0) {
    // Modified Gram-Schmidt orthogonalization for state against states
    // In place
    // statesRe and statesIm are Float64Array. state is [Float64Array, Float64Array]
    const n = statesRe[0].length;
    const i = statesRe.length;

    for (let j = 0; j < i; j++) {
        let dot = dotProduct([statesRe[j], statesIm[j]], state);
        for (let idx = 0; idx < n; idx++) {
            const a = statesRe[j][idx];
            const b = statesIm[j][idx];

            state[0][idx] -= dot[0] * a - dot[1] * b;
            state[1][idx] -= dot[0] * b + dot[1] * a;
        }
    }

    const projectedNorm = norm(state);

    if (normal === 0) {
        for(let idx = 0; idx < size; idx++) {
            state[0][idx] /= projectedNorm;
            state[1][idx] /= projectedNorm;
        }
    }

    return projectedNorm;
}

class SparseMatrix {
    constructor(shape, rows, cols, vals) {
        const [m, n] = shape;

        if (!(rows.length === cols.length && cols.length === vals.length )) {
            throw new Error("Unmatching Input Lengths.");
        }

        this.m = m;
        this.n = n;
        const entries = [];

        // Store all entries into an array
        for (let i = 0; i < vals.length; i++) {
            const r = rows[i];
            const c = cols[i];
            const v = vals[i];
            
            if (r < 0 || r >= m || c < 0 || c >= n) {
                continue;
            }

            if (v !== 0) {
                entries.push({r, c, v});
            }
        }
        
        // Sort entries by row, then column
        entries.sort((a, b) => {
            if (a.r !== b.r) return a.r - b.r;
            return a.c - b.c;
        })

        // Fill up Values, Column Indices, and Row Pointers. Row pointers show where to access each row.
        this.values = [];
        this.columnIndices = [];
        this.rowPointers = [];

        // Initialize rowPointers with starting pointer for row 0
        this.rowPointers.push(0);
        let currentRow = 0;

        for (const {r, c, v} of entries) {
            // Advance row pointers for any empty rows between currentRow and r
            while (r > currentRow) {
                this.rowPointers.push(this.values.length);
                currentRow ++;
            }
            // Now r === currentRow
            this.values.push(v);
            this.columnIndices.push(c);
        }

        // Fill the rest of the empty row pointers after all the data has been entered
        while (currentRow < this.m) {
            this.rowPointers.push(this.values.length);
            currentRow ++;
        }

        this.values = new Float64Array(this.values);
        this.columnIndices = new Int32Array(this.columnIndices);
        this.rowPointers = new Int32Array(this.rowPointers);
   }

    vectorMultiplyComplex(state) {
        const [realX, imagX] = state;
        if (!(realX.length === this.n && imagX.length === this.n)) {
            throw new Error("Incompatible Vector.");
        }

        const realY = new Float64Array(this.m);
        const imagY = new Float64Array(this.m);

        for (let r = 0; r < this.m; r++) {
            let realSum = 0.0;
            let imagSum = 0.0;
            const start = this.rowPointers[r];
            const end = this.rowPointers[r+1];

            for (let i = start; i < end; i++) {
                const value = this.values[i];
                const col = this.columnIndices[i];
                realSum += value * realX[col];
                imagSum += value * imagX[col];
            }
            realY[r] = realSum;
            imagY[r] = imagSum;
        }
        return [realY, imagY];
    }

    expMultiply(state, imagCoefficient) { // Taylor, Legacy
        const [realX, imagX] = state;
        if (!(this.m === this.n)) {
            throw new Error ("Matrix Must Be Square.");
        }

        const n = this.n;
        const terms = 100;

        const realY = new Float64Array(n);
        const imagY = new Float64Array(n);

        let realCurrent = new Float64Array(realX);
        let imagCurrent = new Float64Array(imagX);

        for (let idx = 0; idx < n; idx++) {
            realY[idx] = realCurrent[idx];
            imagY[idx] = imagCurrent[idx];
        }

        // build series terms: for k>=1
        for (let k = 1; k <= terms; k++) {
            const [hReal, hImag] = this.vectorMultiplyComplex([realCurrent, imagCurrent]);

            const scale = imagCoefficient / k;
            for (let idx = 0; idx < n; idx++) {
                const a = hReal[idx];
                const b = hImag[idx];

                const sReal = b * scale;
                const sImag = -a * scale;

                realCurrent[idx] = sReal;
                imagCurrent[idx] = sImag;

                realY[idx] += realCurrent[idx];
                imagY[idx] += imagCurrent[idx];
            }
        }
        console.log(Math.max(...realCurrent))

        return [realY, imagY];
    }

    lanczosReduce(state) {
        // Return d x n orthonormal basis matrix Q and alphas and betas
        const n = this.n;

        const Qre = new Array(d);
        const Qim = new Array(d);
        for (let i = 0; i < d; i++) {
            Qre[i] = new Float64Array(n);
            Qim[i] = new Float64Array(n);
        }
        
        const beta0 = norm(state);
        for (let i = 0; i < n; i++) {
            Qre[0][i] = state[0][i] / beta0
            Qim[0][i] = state[1][i] / beta0
        }

        const alphas = new Float64Array(d);
        const betas = new Float64Array(d);
        for (let i = 0; i < d; i++) {
            let w = this.vectorMultiplyComplex([Qre[i], Qim[i]]);
            // MGS twice
            alphas[i] = dotProduct([Qre[i], Qim[i]], w)[0];
            oneMGS(w, Qre, Qim, 1);
            alphas[i] += dotProduct([Qre[i], Qim[i]], w)[0];
            betas[i] = oneMGS(w, Qre, Qim);
            if (i === d - 1) {
                break;
            }
            for (let idx = 0; idx < n; idx++) {
                Qre[i+1][idx] = w[0][idx];
                Qim[i+1][idx] = w[1][idx];
            }
        }

        // Construct tridiagonal matrix T. All elements should be real, as the Hamiltonian is Hermitian. Q and T are dense now as they are less fat
        const T = [];
        for (let i = 0; i < d; i++) {
            T.push((new Array(d)).fill(0));
        }
        for (let i = 0; i < d; i++) {
            T[i][i] = alphas[i];
            if (i < d - 1) {
                T[i][i+1] = betas[i];
                T[i+1][i] = betas[i];
            }
        }

        return [[Qre, Qim], T];
    }
}


function gaussian(x0, p0, standard_dev) {
    const re = new Float64Array(size);
    const im = new Float64Array(size);
    for(let i = 0; i < size; i++) {
        let standardDistance = (i-x0)/standard_dev;
        let density = Math.exp(-0.5*standardDistance*standardDistance);
        re[i] = density * Math.cos(p0 * (i-x0));
        im[i] = density * Math.sin(p0 * (i-x0)); 
    }
    const state = [re, im];
    normalize(state);
    return state;
}

function createHamiltonian(m, potential) {
    const diagValues = Float64Array.from(potential);
    const sideValues = new Float64Array(size-1);
    const range = new Int32Array(size);

    for(let i = 0; i < size; i++) {
        diagValues[i] += -2 * (-0.5/m)
        range[i] = i
    }
    for(let i = 0; i < size-1; i++) {
        sideValues[i] = 1 * (-0.5/m)
    }

    const lowerRange = range.slice(0,size-1);
    const upperRange = range.slice(1,size);

    const Vals = [...diagValues, ...sideValues, ...sideValues];
    const Rows = [...range, ...lowerRange, ...upperRange];
    const Cols = [...range, ...upperRange, ...lowerRange]

    return new SparseMatrix([size, size], Rows, Cols, Vals);
}

function evolve(state0, stateHamiltonian, time) {
    const newState = stateHamiltonian.expMultiply(state0, time);
    normalize(newState);
    return newState;
}


function InitializeLanczos(state0, m, potential) {
    const hamiltonian = createHamiltonian(m, potential);
    const [Q, T] = hamiltonian.lanczosReduce(state0);

    // Eigendecompose T, only have to exp the eigenvalues
    const eig = math.eigs(T);
    const lambdas = eig.values.map(v => (v && typeof v.re === 'number') ? v.re : +v);
    const S = eig.eigenvectors.map(obj => obj.vector);

    return [state0, Q, lambdas, S];
}

function evolveLanczos(initializer, time) {
    const [state0, Q, lambdas, S] = initializer;
    const d = lambdas.length;
    
    const re = new Float64Array(d); 
    const im = new Float64Array(d);

    // Compute y = e^itA = Q e^itT Qt x0 = Q e^itT e1 |x| = |x| Q S e^itL S^t e1
    for (let idx = 0; idx < d; idx++) {
        re[idx] = S[idx][0]; // St e1
        im[idx] = re[idx] * -1 * Math.sin(lambdas[idx] * time); 
        re[idx] = re[idx] * Math.cos(lambdas[idx] * time); 
    }
    // e^itL St e1

    const re2 = new Float64Array(d);
    const im2 = new Float64Array(d);
    for (let i = 0; i < d; i++) {
        let realSum = 0.0;
        let imagSum = 0.0;
        for (let idx = 0; idx < d; idx++) {
            realSum += S[idx][i] * re[idx];
            imagSum += S[idx][i] * im[idx];
        }
        re2[i] = realSum;
        im2[i] = imagSum;
    }
    // S e^itL St e1

    const Yre = new Float64Array(size);
    const Yim = new Float64Array(size);
    const normx0 = norm(state0);

    for (let i = 0; i < size; i++) {
        let realSum = 0.0;
        let imagSum = 0.0;
        for (let idx = 0; idx < d; idx++) {
            realSum += Q[0][idx][i] * re2[idx] - Q[1][idx][i] * im2[idx];
            imagSum += Q[0][idx][i] * im2[idx] + Q[1][idx][i] * re2[idx];
        }
        Yre[i] = realSum * normx0;
        Yim[i] = imagSum * normx0;
    }
    // |x| Q e^itT e1

    return [Yre, Yim];
}