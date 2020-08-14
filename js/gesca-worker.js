importScripts("jstat.min.js");
importScripts("https://cdnjs.cloudflare.com/ajax/libs/mathjs/7.1.0/math.min.js");
importScripts("underscore-min.js");
let parsedData;
let summary;
/**
 * Create an array with {num} elements with {val} value
 * @param num: number of elements
 * @param val: value of each element
 * @returns {[]}
 */
function create1Darray(num, val) {
    let ret = [];
    for (let i = 0; i < num; i++) {
        ret.push(val)
    }
    return ret;
}

function parseData(data, model) {
    function create2DarrayOfValue(row, col, val) {
        let arr = [];
        if (row === 1) {
            for (let j = 0; j < col; j++) {
                arr.push(val)
            }
        } else {
            for (let i = 0; i < row; i++) {
                let _row = [];
                for (let j = 0; j < col; j++) {
                    _row.push(val)
                }
                arr.push(_row)
            }
        }

        return arr;
    }

    let ret = {};
    let temp = {};
    let colnames = Object.keys(data[0]);
    let splitLines = model.split('\n');
    for (let t = 0; t < splitLines.length; t++) {
        if (splitLines[t].startsWith("#") || splitLines[t].length === 0) {
            splitLines.splice(t, 1);
            t--;
        }
    }
    temp.measurement = {};
    for (let t = 0; t < splitLines.length; t++) {
        if (splitLines[t].indexOf("=~") >= 0) {
            let sp = splitLines[t].split("=~");
            if (temp.measurement.hasOwnProperty(sp[0])) {
                console.log(`Duplicate latent variable ${sp[0]}`)
            } else {
                temp.measurement[sp[0]] = sp[1];
                splitLines.splice(t, 1);
                t--;
            }
        }
    }


    //Extract latent variable
    temp.latent = Object.keys(temp.measurement);
    temp.indicators = [];
    //Extract indicator

    temp.latent.forEach(key => {
        let _indicators = temp.measurement[key].split("+");
        temp.indicators = _.union(temp.indicators, _indicators)
    });
    //Check if indicators do not match columns name
    let indicatorUnmatch = _.difference(temp.indicators, colnames);
    if (indicatorUnmatch.length > 0) {
        if (self.debug) {
            console.log(`Indicators does not match in data ${indicatorUnmatch.join(", ")}`)
        }
    }
    ret.W00 = create2DarrayOfValue(temp.indicators.length, temp.latent.length, 0);
    //Update values in W00
    temp.latent.forEach((value, index) => {
        let _indicators = temp.measurement[value].split("+");
        _indicators.forEach((val2) => {
            let rowIndex = temp.indicators.indexOf(val2);
            ret.W00[rowIndex][index] = 99;
        })
    });

    ret.C00 = math.transpose(ret.W00);


    //Construct B00
    temp.structuralmodel = {};
    for (let t = 0; t < splitLines.length; t++) {
        if (splitLines[t].indexOf("~") >= 0) {
            let sp = splitLines[t].split("~");
            if (temp.structuralmodel.hasOwnProperty(sp[0])) {
                console.log(`Duplicate latent variable ${sp[0]}`)
            } else {
                temp.structuralmodel[sp[0]] = sp[1];
                splitLines.splice(t, 1);
                t--;
            }
        }
    }
    ret.B00 = create2DarrayOfValue(temp.latent.length, temp.latent.length, 0);
    for (let key of Object.keys(temp.structuralmodel)) {
        let c = temp.latent.indexOf(key);
        let r = temp.latent.indexOf(temp.structuralmodel[key]);
        ret.B00[r][c] = 99;
    }

    //Process return data
    ret.data = create2DarrayOfValue(data.length, temp.indicators.length, 0);
    data.forEach(function (d, i) {
        temp.indicators.forEach((val, id) => {
            ret.data[i][id] = +d[val];
        })
    });
    ret.latent = temp.latent;
    ret.indicators = temp.indicators;

    return ret;
}

function cloneArray(source, target, Rstart, Cstart) {
    let dim = jStat(source).dimensions();
    let rows = dim.rows, cols = dim.cols;
    for (let i = 0; i < rows; i++) {
        for (let j = 0; j < cols; j++) {
            target[i + Rstart][j + Cstart] = source[i][j];
        }
    }
    return target;
}

function analyzeModel(params) {
    console.log("Initialization....");
    let z0 = self.parsedData.data; //raw data
    let dim = jStat(z0).dimensions(); // rows x columns
    let nobs = dim.rows; // number of observations
    let nvar = dim.cols; // number of variables
    let ng = params.ng;// number of groups

    let nobs_g = self.create1Darray(ng, 0); // number of obs per group
    for (let g = 0; g < ng; g++) {
        nobs_g[g] = nobs;
    }

    let case_index = jStat.zeros(ng, 2);
    let _from = -1;
    for (let g = 0; g < ng; g++) {
        let _start = _from + 1;
        let _end = _from + nobs_g[g];
        case_index[g][0] = _start;
        case_index[g][1] = _end;
    }
    //
    let nlv = params.loadtype.length; //number of latent variables
    let totalvars = nlv + nvar; // total variables

    let A00 = math.concat(self.parsedData.C00, self.parsedData.B00);
    let _I = math.identity(nvar)._data;
    let V00 = math.concat(_I, self.parsedData.W00);
    let Wi = self.parsedData.W00;
    let Ai = A00;
    //
    // //Calculate index
    let windex0 = [];
    let aindex0 = [];
    let _W00 = math.flatten(math.transpose(self.parsedData.W00));
    _W00.forEach((_d, _i) => {
        if (_d === 99) {
            windex0.push(_i);
        }
    })
    let _A00 = math.flatten(math.transpose(A00));
    _A00.forEach((_d, _i) => {
        if (_d === 99) {
            aindex0.push(_i);
        }
    })


    let W0 = jStat.zeros(ng * nvar, ng * nlv);
    let A0 = jStat.zeros(ng * nlv, totalvars);
    let V0 = jStat.zeros(ng * nvar, ng * totalvars);
    //
    let W = JSON.parse(JSON.stringify(W0));
    let A = JSON.parse(JSON.stringify(A0));
    let V = JSON.parse(JSON.stringify(V0));

    let kk = -1, ss = -1, ll = -1;
    for (let j = 0; j < ng; j++) {
        let k = kk + 1;
        kk += nvar;
        let s = ss + 1;
        ss += nlv;
        let l = ll + 1;
        ll += totalvars;

        W0 = self.cloneArray(self.parsedData.W00, W0, k, s);
        windex0.forEach(d => {
            let r = d % nvar;
            let c = math.floor(d / nvar);
            Wi[r][c] = math.random();
        })
        W = self.cloneArray(Wi, W, k, s);
        A0 = self.cloneArray(A00, A0, s, 0);

        aindex0.forEach(d => {
            let r = d % jStat(A).dimensions().rows;
            let c = math.floor(d / jStat(A).dimensions().rows);
            Ai[r][c] = math.random();
        })


        A = self.cloneArray(Ai, A, s, 0);
        V0 = self.cloneArray(V00, V0, k, l);
        let I = jStat.identity(nvar);
        let _WI = math.concat(I, Wi);
        V = self.cloneArray(_WI, V, k, l);

    }


    let I = jStat.zeros(totalvars * ng, totalvars);
    kk = -1;
    for (let g = 0; g < ng; g++) {
        let k = kk + 1;
        kk += totalvars;
        let _I = jStat.identity(totalvars);
        I = self.cloneArray(_I, I, k, 0);

    }
    let constrainmatrix = self.generateConstraintMatrix(A0);
    let PHT = constrainmatrix.PHT;
    let num_nzct = constrainmatrix.num_nzct;
    let num_const = constrainmatrix.num_const;

    // Start Bootstrap
    let num_nnz_W00 = math.flatten(self.parsedData.W00).filter(d => d !== 0).length;
    let num_nnz_C00 = math.flatten(self.parsedData.C00).filter(d => d !== 0).length;
    let num_nnz_B00 = math.flatten(self.parsedData.B00).filter(d => d !== 0).length;

    let vec_FIT = jStat.zeros(params.nbt, 1);
    let vec_FIT_m = jStat.zeros(params.nbt, 1);
    let vec_FIT_s = jStat.zeros(params.nbt, 1);
    let vec_AFIT = jStat.zeros(params.nbt, 1);
    let vec_GFI = jStat.zeros(params.nbt, 1);
    let vec_SRMR = jStat.zeros(params.nbt, 1);

    let MatW = jStat.zeros(params.nbt, num_nnz_W00 * ng);
    let Matload = jStat.zeros(params.nbt, num_nnz_C00 * ng);
    let Matbeta = jStat.zeros(params.nbt, num_nnz_B00 * ng);
    let Matsmc = jStat.zeros(params.nbt, num_nnz_C00 * ng);
    let MatcorF = jStat.zeros(params.nbt, Math.pow(nlv, 2) * ng);

    let MatTE_S = [];
    let MatID_S = [];
    let MatTE_M = [];
    let MatID_M = [];
    let Z;
    let WR, Cr, Br, samplesizes, FIT, FIT_S, FIT_M, AFIT, GFI, SRMR, R2, AVE, Alpha, rho, LV_MEAN, LV_VAR,
        corr_corres, Dimension, latencorr, TE_S, ID_S, TE_M, ID_M, mW = [], mC = [], mB = [], mSMC, mCF = [], lb, ub;


    for (let b = 0; b <= params.nbt; b++) {
        self.postMessage({job:'iteration',iter:b})
        if (b === 0) {
            if (params.moption > 1) {

            } else {
                Z = self.bootsample(z0, case_index, nvar, nobs_g, ng, b, nobs);
            }
        } else {
            Z = self.bootsample(z0, case_index, nvar, nobs_g, ng, b, nobs);

        }


        if (b === 0) {
            if (params.moption === 3) {
            } else {
                self.RESULT = self.alternativeLeastSquareAlgorithm(Z, W0, A0, W, A, V, I, PHT, nvar, nlv, ng, params.itmax, params.ceps);
            }
        } else {
            self.RESULT = self.alternativeLeastSquareAlgorithm(Z, W0, A0, W, A, V, I, PHT, nvar, nlv, ng, params.itmax, params.ceps);
        }

        let it = self.RESULT.it;
        let imp = self.RESULT.imp;
        let Gamma = self.RESULT.Gamma;
        let Psi = self.RESULT.Psi;
        let f0 = self.RESULT.f0;
        A = self.RESULT.A;
        let corF = math.multiply(math.transpose(Gamma), Gamma);
        let CR = math.transpose(jStat.col(A, jStat.arange(nvar)));
        let BR = math.transpose(jStat.col(A, jStat.arange(nvar, totalvars)));
        let DF = nobs * nvar;

        let npw = math.flatten(W0).filter(d => d === 99).length;
        let dpht = math.diag(PHT);
        let cnzt;
        if (num_nzct === 0) {
            cnzt = dpht.filter(d => d === 1).length;
        } else {
            cnzt = num_nzct + dpht.filter(d => d === 1).length;
        }
        let NPAR = cnzt + npw;

        let Fit = 1 - f0 / (math.sum(math.diag(math.multiply(math.transpose(Psi), Psi))));
        let dif_m = math.subtract(jStat.col(Psi, jStat.arange(nvar)), math.multiply(Gamma, math.transpose(CR)));
        let dif_s = math.subtract(jStat.col(Psi, jStat.arange(nvar, totalvars)), math.multiply(Gamma, math.transpose(BR)));
        let Fit_m = 1 - math.sum(math.diag(math.multiply(math.transpose(dif_m), dif_m))) / math.sum(math.diag(math.multiply(math.transpose(Z), Z)));
        let Fit_s = 1 - math.sum(math.diag(math.multiply(math.transpose(dif_s), dif_s))) / math.sum(math.diag(corF));

        let Afit = 1 - ((1 - Fit) * (DF) / (DF - NPAR));

        self.modelfit = self.modelFit(Z, self.RESULT.W, A, nvar, nlv, ng, case_index);


        let Gfi = self.modelfit.GFI;
        let Srmr = self.modelfit.SRMR;
        let COR_RES = self.modelfit.COR_RES;

        let total_s = jStat.zeros(ng * nlv, nlv);
        let indirect_s = jStat.zeros(ng * nlv, nlv);
        let total_m = jStat.zeros(ng * nlv, nvar);
        let indirect_m = jStat.zeros(ng * nlv, nvar);
        var k = kk = -1;
        for (let g = 0; g < ng; g++) {
            k = kk + 1;
            kk += nlv;
            let B = jStat.col(BR, jStat.arange(k, kk + 1));
            let C = jStat.col(CR, jStat.arange(k, kk + 1))
            self.effects = self.estimatedEffects(B, C);
            let te_s = self.effects.te_s;
            let te_m = self.effects.te_m;
            let ie_s = self.effects.ie_s;
            let ie_m = self.effects.ie_m;
            total_s = math.subset(total_s, math.index(math.range(k, kk + 1), math.range(0, nlv)), te_s);
            indirect_s = math.subset(indirect_s, math.index(math.range(k, kk + 1), math.range(0, nlv)), ie_s);
            total_m = math.subset(total_m, math.index(math.range(k, kk + 1), math.range(0, nvar)), te_m);
            indirect_m = math.subset(indirect_m, math.index(math.range(k, kk + 1), math.range(0, nvar)), ie_m);

        }
        if (b === 0) {
            if (params.moption === 2) {
                console.log("Missing value")
            } else if (self.moption === 3) {

            }

            if (it <= params.itmax) {
                if (imp <= params.ceps) {
                    console.log(`The algorithm converged in ${it} iterations (convergence criterion = ${params.ceps})`)
                } else {
                    console.log(`The algorithm FAILED to convergence in ${it} iterations (convergence criterion = ${params.ceps})`)
                }
            }
            // WR = JSON.parse(JSON.stringify(W));
            WR = W;
            Cr = CR;
            Br = BR;
            samplesizes = nobs_g;
            FIT = Fit;
            FIT_M = Fit_m;
            FIT_S = Fit_s;
            AFIT = Afit;
            GFI = Gfi;
            SRMR = Srmr;

            R2 = jStat.zeros(ng, nlv);
            AVE = jStat.zeros(ng, nlv);
            Alpha = jStat.zeros(ng, nlv);
            rho = jStat.zeros(ng, nlv);
            Dimension = jStat.zeros(ng, nlv);
            let lvmean = jStat.zeros(ng, nlv);
            let lvvar = jStat.zeros(ng, nlv);
            corr_corres = jStat.zeros(ng * nvar, nvar);

            let ss = -1;
            let kk = -1;
            let z0_g;
            for (let g = 0; g < ng; g++) {
                let s = ss + 1;
                let k = kk + 1;
                ss = ss + nlv;
                kk += nvar;
                if (self.moption === 3) {

                } else {
                    z0_g = jStat.row(z0, jStat.arange(case_index[g][0], case_index[g][1] + 1))
                }
                let W_g = math.subset(W, math.index(math.range(k, kk + 1), math.range(s, ss + 1)));
                let CF_g = math.subset(corF, math.index(math.range(s, ss + 1), math.range(s, ss + 1)));
                let B = math.transpose(jStat.col(BR, jStat.arange(s, ss + 1)));
                let stdL = jStat.col(CR, jStat.arange(s, ss + 1));
                for (let j = 0; j < nlv; j++) {
                    let _val = math.multiply(math.transpose(jStat.col(B, j)), jStat.col(CF_g, j));
                    R2[g][j] = _val[0][0];
                    let zind = [];
                    math.flatten(math.transpose(jStat.col(self.parsedData.W00, j))).forEach((d, i) => {
                        if (d !== 0) zind.push(i);
                    })
                    let nnzload = zind.length;
                    if (nnzload > 0) {
                        let sumload = jStat.sum(math.flatten(jStat.pow(math.subset(stdL, math.index(zind, j)), 2)));
                        let sumload_rho1 = math.pow(jStat.sum(math.flatten(math.subset(stdL, math.index(zind, j)))), 2);
                        let sumload_rho2 = math.sum(math.subtract(1, jStat.pow(math.flatten(math.subset(stdL, math.index(zind, j))), 2)));
                        AVE[g][j] = sumload / nnzload;
                        rho[g][j] = sumload_rho1 / (sumload_rho1 + sumload_rho2);
                    }
                    let nzj = zind.length;
                    if (nzj > 1) {
                        let zsubset = jStat.col(z0_g, zind);
                        Alpha[g][j] = self.cronbachalpha(zsubset);
                        let eigs = math.eigs(self.correlationMat(zsubset)).values;
                        Dimension[g][j] = eigs.filter(d => d > 1).length;
                    } else {
                        Alpha[g][j] = 1;
                    }

                }
                //Calculate latent scores
                let lvscore_g = self.lvscore(z0_g, W_g);
                jStat(lvscore_g).mean().map((d, i) => {
                    lvmean[g][i] = d;
                })
                jStat(lvscore_g).variance().map((d, i) => {
                    lvvar[g][i] = d;
                })
                let sample_corr = self.correlationMat(z0_g);
                //Update Corr_res
                for (let r = k; r < kk + 1; r++) {
                    for (let c = 0; c < nvar; c++) {
                        if (r === c) continue;
                        if (r < c) {
                            corr_corres[r][c] = COR_RES[r][c];
                        } else {
                            corr_corres[r][c] = sample_corr[r][c];
                        }
                    }
                }
            }
            LV_MEAN = lvmean;
            LV_VAR = lvvar;

            let  _mW =[];
            math.flatten(math.transpose(W)).forEach(d => {
                if (d !== 0) _mW.push(d);
            });
            mW.push(_mW);

            let  _mC =[];
            math.flatten(math.transpose(Cr)).forEach(d => {
                if (d !== 0) _mC.push(d);
            });
            mC.push(_mC);
            mSMC = jStat.pow(mC,2);
            let _mCF =[];
            math.flatten(math.transpose(corF)).forEach(d => {
                if (d !== 0) _mCF.push(d);
            });
            mCF.push(_mCF);

            let  _mB =[];
            math.flatten(math.transpose(Br)).forEach(d => {
                if (d !== 0) _mB.push(d);
            });
            mB.push(_mB);

            latencorr = corF;
            TE_S = total_s;
            ID_S = indirect_s;
            TE_M = total_m;
            ID_M = indirect_m;

            ////

        }
        else { //Bootstrap sample solution

            let vecw = [];
            math.flatten(math.transpose(W)).forEach(d => {
                if (d !== 0) vecw.push(d);
            });
            let vecload = [];
            math.flatten(math.transpose(CR)).forEach(d => {
                if (d !== 0) vecload.push(d);
            });
            let vecbeta = [];
            math.flatten(math.transpose(BR)).forEach(d => {
                if (d !== 0) vecbeta.push(d);
            });
            let veccorF = [];
            math.flatten(math.transpose(corF)).forEach(d => {
                if (d !== 0) veccorF.push(d);
            });

            vec_FIT[b - 1] = Fit;
            vec_FIT_m[b - 1] = Fit_m;
            vec_FIT_s[b - 1] = Fit_s;
            vec_AFIT[b - 1] = Afit;
            vec_GFI[b - 1] = Gfi;
            vec_SRMR[b - 1] = Srmr;

            vecw.map((d, i) => {
                MatW[b - 1][i] = d;
            })
            vecload.map((d, i) => {
                Matload[b - 1][i] = d;
            })
            vecbeta.map((d, i) => {
                Matbeta[b - 1][i] = d;
            })
            vecload.map((d, i) => {
                Matsmc[b - 1][i] = d * d;
            });
            veccorF.map((d, i) => {
                MatcorF[b - 1][i] = d;
            });
            let _temptt_s = [];
            math.flatten(math.transpose(total_s)).forEach(d => {
                if (d !== 0) _temptt_s.push(d);
            })
            MatTE_S.push(_temptt_s);

            let _indirect_s = [];
            math.flatten(math.transpose(indirect_s)).forEach(d => {
                if (d !== 0) _indirect_s.push(d);
            })
            MatID_S.push(_indirect_s);

            let _total_m = [];
            math.flatten(math.transpose(total_m)).forEach(d => {
                if (d !== 0) _total_m.push(d);
            })
            MatTE_M.push(_total_m);

            let _indirect_m = [];
            math.flatten(math.transpose(indirect_m)).forEach(d => {
                if (d !== 0) _indirect_m.push(d);
            })
            MatID_M.push(_indirect_m);

        }



    }
    lb = math.ceil(params.nbt * 0.025);
    ub = math.ceil(params.nbt * 0.975);
    let sortFIT = vec_FIT.sort((a, b) => a - b);
    let sortFIT_m = vec_FIT_m.sort((a, b) => a - b);
    let sortFIT_s = vec_FIT_s.sort((a, b) => a - b);
    let sortAFIT = vec_AFIT.sort((a, b) => a - b);
    let sortGFI = vec_GFI.sort((a, b) => a - b);
    let sortSRMR = vec_SRMR.sort((a, b) => a - b);
    let sortw = math.transpose(math.transpose(MatW).map(d => d.sort()));
    let sortload = math.transpose(math.transpose(Matload).map(d => d.sort()));
    let sortbeta = math.transpose(math.transpose(Matbeta).map(d => d.sort()));
    let sortsmc = math.transpose(math.transpose(Matsmc).map(d => d.sort()));
    let sortcorF = math.transpose(math.transpose(MatcorF).map(d => d.sort()));
    let sortte_s = math.transpose(math.transpose(MatTE_S).map(d => d.sort()));
    let sortid_s = math.transpose(math.transpose(MatID_S).map(d => d.sort()));
    let sortte_m = math.transpose(math.transpose(MatTE_M).map(d => d.sort()));
    let sortid_m = math.transpose(math.transpose(MatID_M).map(d => d.sort()));

    return {
        mW,
        mC,
        mB,
        MatW,
        Matload,
        Matbeta,
        Matsmc,
        MatcorF,
        sortAFIT,
        sortbeta,
        sortcorF,
        sortFIT,
        sortFIT_m,
        sortFIT_s,
        sortGFI,
        sortid_m,
        sortid_s,
        sortload,
        sortw,
        sortSRMR,
        sortsmc,
        sortte_m,
        sortte_s,
        WR,
        Cr,
        Br,
        samplesizes,
        FIT,
        FIT_M,
        FIT_S,
        AFIT,
        GFI,
        SRMR,
        R2,
        AVE,
        Alpha,
        rho,
        LV_MEAN,
        LV_VAR,
        corr_corres,
        Dimension,
        latencorr,
        lb,
        ub,
        vec_SRMR,
        vec_GFI,
        vec_AFIT,
        vec_FIT_s,
        vec_FIT_m,
        vec_FIT
    }
}
function lvscore(X, W, option = 1) {
    let nobs = jStat.rows(X);
    let nlv = jStat.cols(W);
    let ctz = math.subtract(X, math.multiply(math.ones(nobs, 1), [math.mean(X, 0)]))._data;
    let covz = math.divide(math.multiply(math.transpose(ctz), ctz), nobs);
    let Dstz = math.sqrt(math.diag(math.diag(covz)));
    let UstdW = math.multiply(math.inv(Dstz), W);
    let lvscore;
    if (option === 1) {
        let sumUstdW = jStat(UstdW).sum();
        let rUstdW = JSON.parse(JSON.stringify(UstdW));
        for (let j = 0; j < nlv; j++) {
            let _val = math.flatten(math.divide(math.subset(UstdW, math.index(math.range(0, jStat.rows(UstdW)), j)), sumUstdW[j]));
            _val.map((d, i) => {
                rUstdW[i][j] = d;
            })
        }
        lvscore = math.multiply(X, rUstdW);

    } else if (option === 2) {
        lvscore = math.multiply(X, UstdW);
    }
    return lvscore;

}

function cronbachalpha(X) {
    let nvar = jStat.cols(X);
    let crX = self.correlationMat(X);
    let rows = jStat.rows(crX);
    let cols = jStat.cols(crX);
    let sumcorr = 0;
    for (let r = 0; r < rows; r++) {
        for (let c = r + 1; c < cols; c++) {
            sumcorr += crX[r][c];
        }
    }
    return nvar * sumcorr / (0.5 * nvar * (nvar - 1) + (nvar - 1) * sumcorr)


}
function estimatedEffects(B,C){
    //Calculate the total and indirect effects in measurement and structural model
    let nlv = jStat.rows(B);
    let te_s = math.transpose(math.subtract(math.inv(math.subtract(math.identity(nlv), B)), math.identity(nlv)));
    let ie_s = math.subtract(te_s, math.transpose(B));
    let te_m = math.transpose(math.transpose(math.multiply(math.inv(math.transpose(math.subtract(math.identity(nlv), B))), math.transpose(C))))
    let ie_m = math.subtract(te_m, math.transpose(C));
    return {te_s: te_s._data, ie_s: ie_s._data, te_m: te_m._data, ie_m: ie_m._data}
}
function correlationMat(data) {
    let nrow = jStat.rows(data);
    let C = math.subtract(math.identity(nrow), math.divide(math.ones(nrow, nrow), nrow));
    let D = math.diag(jStat(data).stdev());
    let Xs = math.multiply(C, data, math.inv(D))
    return math.divide(math.multiply(math.transpose(Xs), Xs), nrow)._data;
}
function modelFit(Z0, W, T, nvar, nlv, ng, case_index) {
    let GFI, SRMR, COR_RES;
    let obs = jStat.rows(Z0);

    let gfi_1 = 0;
    let gfi_2 = 0;
    let srmr_1 = 0;
    let kk = -1;
    let ss = -1;
    for (let g = 0; g < ng; g++) {
        let k = kk + 1;
        kk += nvar;
        let s = ss + 1;
        ss += nlv;
        let zz = math.subset(Z0, math.index(math.range(case_index[g][0], case_index[g][1] + 1), math.range(k, kk + 1)));
        let w = math.subset(W, math.index(math.range(k, kk + 1), math.range(s, ss + 1)));
        let v = math.concat(math.identity(nvar), w)._data;
        let t = jStat.row(T, jStat.arange(s, ss + 1));
        let omega = math.subtract(v, math.multiply(w, t));
        let ee = math.multiply(zz, omega);
        let samcov = math.divide(math.multiply(math.transpose(math.subtract(zz, math.divide(math.multiply(math.ones(obs, obs), zz), obs))), math.subtract(zz, math.divide(math.multiply(math.ones(obs, obs), zz), obs))), (obs - 1))._data;
        let samcor = self.correlationMat(zz);
        let tp_precov = math.multiply(math.inv(math.multiply(omega, math.transpose(omega))), math.multiply(omega, math.diag(jStat(ee).variance()), math.transpose(omega)))
        let precov = math.transpose(math.multiply(math.inv(math.transpose(math.multiply(omega, math.transpose(omega)))), math.transpose(tp_precov)))
        let COV_RES = math.subtract(samcov, precov);
        let prerij = JSON.parse(JSON.stringify(precov));
        for (let i = 0; i < nvar; i++) {
            for (let j = 0; j < nvar; j++) {
                prerij[i][j] = precov[i][j] / math.sqrt(precov[i][i] * precov[j][j]);
            }
        }
        let srmr = 0;
        for (let i = 0; i < nvar; i++) {
            for (let j = 0; j < nvar; j++) {
                if (j > i) {
                    let corr_residual = math.pow(samcor[i][j] - prerij[i][j], 2);
                    srmr = srmr + corr_residual;
                }
            }
        }

        srmr_1 += srmr;
        gfi_1 += math.sum(math.diag(math.pow(COV_RES, 2)));
        gfi_2 += math.sum(math.diag(math.pow(samcov, 2)));

        COR_RES = math.subtract(samcor, prerij);

    }
    let nvar_tot = ng * nvar;
    let srmr_2 = nvar_tot * (nvar_tot + 1) / 2;
    SRMR = math.sqrt(srmr_1 / srmr_2);
    GFI = 1 - (gfi_1 / gfi_2);

    return {GFI, SRMR, COR_RES}
}
function alternativeLeastSquareAlgorithm(Z, W0, A0, W, A, V, I, PHT, nvar, nlv, ng, itmax, ceps) {

    let numberOfColumnsA = math.matrix(A).size()[1];
    let numberOfRowsA = math.matrix(A).size()[0];
    let aindex0 = [];


    math.flatten(math.transpose(A0)).forEach((d, i) => {
        if (d >= 1) aindex0.push(i)
    });
    let Psi = math.multiply(Z, V, I);
    let Gamma = math.multiply(Z, W);
    let it = 0;
    let imp = 100000;
    let f0 = math.pow(10, 10);
    while (it <= itmax && imp > ceps) {
        it += 1;

        //Update A
        for (let t = 0; t < numberOfColumnsA; t++) {
            if (math.sum(math.flatten(jStat.col(A0, t))) !== 0) {
                let H1 = jStat.identity(numberOfColumnsA);
                H1[t][t] = 0;
                let aindex = [];
                math.flatten(jStat.col(A0, t)).forEach((d, i) => {
                    if (d >= 1) aindex.push(i);
                })
                if (aindex.length > 0) {
                    let a = jStat.col(A, t);
                    aindex.map(d => {
                        a[d][0] = 0;
                    })
                    let e = jStat.zeros(1, numberOfColumnsA);
                    e[0][t] = 1;
                    let Y = math.subtract(Psi, math.multiply(Gamma, math.add(math.multiply(A, H1), math.multiply(a, e))))
                    let X = jStat.col(Gamma, aindex);
                    let _val = math.multiply(math.inv(math.multiply(math.transpose(X), X)), math.multiply(math.transpose(X), Y, math.transpose(e)));
                    aindex.forEach((ai, i) => {
                        A[ai][t] = _val[0][i];
                    })
                }
            }
        }

        let vecA = self.extractValFromIndex(A, aindex0);
        let _tempAA = math.flatten(math.transpose(A));
        aindex0.forEach((d, i) => {
            _tempAA[d] = vecA[i];
        })
        A = math.transpose(math.reshape(_tempAA, [numberOfColumnsA, numberOfRowsA]));


        let kk = -1;
        for (let g = 0; g < ng; g++) {
            let k = kk + 1;
            kk = kk + nlv;
            let p = (g + 1) * nvar + g * nlv;
            let s = -1;
            for (let j = k; j <= kk; j++) {
                s = s + 1;
                let windex = [];
                math.flatten(jStat.col(W0, j)).forEach((d, i) => {
                    if (d === 99) windex.push(i);
                })
                let w = jStat.col(W, j);
                windex.map(d => {
                    w[d][0] = 0;
                });

                let beta = [math.subtract(jStat.row(I, p + s), jStat.row(A, j))];//1D array
                let H2 = jStat.identity(ng * nlv);
                H2[j][j] = 0;
                let H3 = jStat.identity(numberOfColumnsA * ng);
                H3[p + s][p + s] = 0;
                let AA1 = math.multiply(W, H2, A);
                let AA2 = math.multiply(V, H3, I);
                let AA3 = math.multiply(w, beta);
                let Delta = math.subtract(math.subtract(AA1, AA2), AA3);
                let Zp = jStat.col(Z, windex);
                let temptheta = math.multiply(math.multiply(math.inv(math.multiply(math.transpose(Zp), Zp)), math.transpose(Zp)), math.multiply(Z, Delta), math.transpose(beta))
                let theta = math.transpose(math.multiply(math.inv(math.transpose(math.multiply(beta, math.transpose(beta)))), math.transpose(temptheta)))
                let zw = math.multiply(Zp, theta);
                theta = math.multiply(theta, math.inv(math.sqrt(math.multiply(math.transpose(zw), zw))));
                windex.map((d, i) => {
                    W[d][j] = theta[i][0];
                    V[d][p + s] = theta[i][0];
                })
            }
        }
        Gamma = math.multiply(Z, W);
        Psi = math.multiply(Z, V);
        let dif = math.subtract(Psi, math.multiply(Gamma, A));
        let f = math.trace(math.multiply(math.transpose(dif), dif));
        imp = f0 - f;
        f0 = f;
    }
    return {W, A, Psi, Gamma, f0, it, imp}
}
function extractValFromIndex(data, IndexArr) {
    let ret = [];
    let A00 = math.flatten(math.transpose(data));
    IndexArr.forEach(d => {
        ret.push(A00[d])
    })
    return ret;
}
function extractFromArray(data, from, to) {
    let ret = [];
    for (let i = from; i <= to; i++) {
        ret.push(data[i]);
    }
    return ret;
}

function bootsample(z0, case_index, nvar, number_of_observation_per_group, ng, b, number_of_observations) {
    let Z = jStat.zeros(number_of_observations, ng * nvar);
    let bz01;
    let kk = -1;
    for (let g = 0; g < ng; g++) {
        let k = kk + 1;
        kk += nvar;
        if (b === 0) {
            bz01 = self.extractFromArray(z0, case_index[g][0], case_index[g][1])
        } else {
            let rrb = math.ceil(math.multiply(math.random([1, number_of_observations]), number_of_observations - 1))[0];
            bz01 = jStat.zeros(number_of_observation_per_group[g], nvar);
            let z0_g = jStat.row(z0, jStat.arange(case_index[g][0], case_index[g][1] + 1));
            for (let i = 0; i < number_of_observation_per_group[g]; i++) {
                bz01[i] = JSON.parse(JSON.stringify(z0_g[rrb[i]]))

                // math.subset(bz01,math.index(i, math.range(0,nvar)), math.subset(z0_g,math.index(rrb[0][i], math.range(0, nvar))) );
            }
        }
        let _z0mean = [jStat(z0).mean()];
        let _identity = jStat.ones(number_of_observation_per_group[g], 1);
        let _product = math.multiply(_identity, _z0mean);
        let ctz = math.subtract(bz01, _product);
        let covz = math.divide(math.multiply(math.transpose(ctz), ctz), number_of_observation_per_group[g]);
        let Dstz = math.sqrt(math.diag(math.diag(covz)));
        let bz0 = math.divide(ctz, Dstz);
        bz0 = math.divide(bz0, math.sqrt(number_of_observation_per_group[g]));
        Z = self.cloneArray(bz0, Z, case_index[g][0], k);


    }
    return Z;
}

function generateConstraintMatrix(A0) {
    let PHT, num_const;
    let dim = jStat(A0).dimensions();
    let vector0 = math.flatten(math.transpose(A0));
    let vect = [];
    let nzct = [];
    vector0.forEach((d, i) => {
        if (d >= 1) {
            vect.push(d)
            nzct.push(i);
        }
    })
    let nzt = nzct.length;
    let nzcst = [];
    let num_nzct = nzcst.length;
    if (num_nzct === 0) {
        PHT = jStat.identity(nzt);
        num_const = 0;
    } else {

    }
    return {PHT, num_nzct, num_const}
}

self.addEventListener('message', ev => {
    let job = ev.data.job;
    let params = ev.data.params;
    switch (job) {
        case 'parseData':
            self.parsedData = self.parseData(ev.data.jdata, ev.data.model);
            self.postMessage({job: '_parsedData',parsedData: self.parsedData});
            break;
        case 'analyzeModel':
            self.summary  = self.analyzeModel(params);
            self.postMessage({job:'summary',ret:self.summary});
            break;
        default:
            console.log("Invalid access");
            self.postMessage('Closing web worker !');
            self.close();

    }

})