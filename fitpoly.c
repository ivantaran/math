
#include "fitpoly.h"

#include "statistic.h"
#include <math.h>
#include <stdlib.h>
#include "matrix.h"

double poly(double t, struct fplyset *s);
void scale(struct fplyset *s);
void fillx(double *buffer, double mean, double std, int n);
void filly(double *x, double *buffer, int degree, int n);
void resids(struct fplyset *s);
void statistic(struct fplyset *s);
int weights(double *r, double *w, double mean, int len);

void scale(struct fplyset *s) {
    s->smean = scale_mean(s->len);
    s->sstd = scale_std(s->smean, s->len);
}

double poly(double t, struct fplyset *s) {
    int i;
    double result = s->k[0];
    
    for (i = 1; i < s->n; ++i)
        result += s->k[i]*pow(t, i);
    
    return result;
}

void fillx(double *buffer, double mean, double std, int n) {
    int i;
    
    for (i = 0; i < n; ++i)
        buffer[i] = ((double)i - mean)/std;
}

void filly(double *x, double *buffer, int degree, int n) {
    int i;
    
    for (i = 0; i < n; ++i) {
        buffer[i] = pow(x[i], degree);
    }
}

void resids(struct fplyset *s) {
    int i;
    double t;
    
    for (i = 0; i < s->len; ++i) {
        t = (i - s->smean)/s->sstd;
        s->fx[i] = poly(t, s);
        s->resids[i] = s->data[i] - s->fx[i];
    }
}

void statistic(struct fplyset *s) {
    s->mean = mean(s->resids, s->len);
    s->rmse = stddev(s->resids, s->mean, s->len);
    s->sse = sse(s->resids, s->len);
}

int fitpoly(struct fplyset *s) {
    int i;
    int result;
    double   **m = (double **)mtx_create(s->n, s->len);
    double **qim = (double **)mtx_create(s->n, s->len);
    double   **q = (double **)mtx_create(s->n, s->n  );
    double  **qi = (double **)mtx_create(s->n, s->n  );
    double *x = (double *)malloc(s->len*sizeof(double));

    scale(s);
    
    fillx(x, s->smean, s->sstd, s->len);
    
    for (i = 0; i < s->n; ++i)
        filly(x, m[i], i, s->len);
    
    mtx_sqr(m, q, s->n, s->len);
    
//    vcr_save(x, s->len, "x.txt");
//    mtx_save(m, s->n, s->len, "m.txt");
//    mtx_save(q, s->n, s->n, "q.txt");
    
    result = mtx_hinv(q, qi, s->n);

    if (result) {
        mtx_mul(qi, m, qim, s->n, s->n, s->len);
        mvr_mul(qim, s->data, s->k, s->n, s->len);
    }
    
    mtx_free(m  , s->n);
    mtx_free(q  , s->n);
    mtx_free(qi , s->n);
    mtx_free(qim, s->n);
    free(x);
    
    resids(s);
    statistic(s);
    
    return result;
}

int weights(double *r, double *w, double mn, int len) {
    int i;
    int err = 0;
//    int mxi;
    double m;
    mn = mean(r, len);
    m = mad(r, mn, len)*6.0;
    
    for (i = 0; i < len; ++i) {
        if (fabs(r[i]) > m) {
            w[i] = 0;
            err++;
        }
        else {
            if (w[i] > 0) {
                w[i] = pow(1 - pow(r[i]/m, 2), 2);
            }
        }
    }
    
//    mxi = vcr_maxi(r, len);
//    
//    if (r[mxi] >= m) {
//        w[mxi] = 0;
//        printf("i = %d\n", mxi);
//        return 1;
//    }
    
    return err;
}

int rfitpoly(struct fplyset *s) {
    int i, j;
    int result;
    int err = 0;
    int n = 0;
    
    double   **m = (double **)mtx_create(s->len, s->n);
    double   **mw = (double **)mtx_create(s->n, s->len);
    double **qim = (double **)mtx_create(s->n, s->len);
    double   **q = (double **)mtx_create(s->n, s->n  );
    double  **qi = (double **)mtx_create(s->n, s->n  );
    double *x = (double *)malloc(s->len*sizeof(double));
    double *w = (double *)malloc(s->len*sizeof(double));
    double *r = (double *)malloc(s->len*sizeof(double));

    vcr_fill(w, 1.0, s->len);
    scale(s);
    fillx(x, s->smean, s->sstd, s->len);
    
    for (i = 0; i < s->n; ++i)
        filly(x, mw[i], i, s->len);
    
    mtx_trp(mw, m, s->n, s->len);
    
    for (j = 0; j < 6; ++j) {
        vrm_mul(w, m, mw, s->len, s->n);

        mtx_mul(mw, m, q, s->n, s->len, s->n);

        result = mtx_hinv(q, qi, s->n);
        
        if (result) {
            mtx_mul(qi, mw, qim, s->n, s->n, s->len);
            mvr_mul(qim, s->data, s->k, s->n, s->len);
        }
        
        resids(s);
        statistic(s);
        vcr_mul(s->resids, w, r, s->len);
//        vcr_abs(r, r, s->len);
//        n = weights(r, w, s->mean, s->len);
        n = weights(r, w, s->mean, s->len);
//        if (n == err) {
//            break;
//        }
//        else {
//            err = n;
////            vcr_save(w, s->len, "w.txt");
//        }
    }
    
    mtx_free(mw  , s->n);
    mtx_free(m  , s->len);
    mtx_free(q  , s->n);
    mtx_free(qi , s->n);
    mtx_free(qim, s->n);
    free(x);
    free(w);
    free(r);
    
    return result;
}

double rloess_point(double* data, int len, int* fail) {
    int n;
    double result = 0;
    double p[3];
    
    len += len % 2 - 1;
    
    double *fx = (double *)malloc(len*sizeof(double));
    double *resids = (double *)malloc(len*sizeof(double));
    
    struct fplyset sp;

    sp.data = data;
    sp.len = len;
    sp.n = 3;
    sp.fx = fx;
    sp.k = p;
    sp.resids = resids;
    
    n = rfitpoly(&sp);
    
    if (!n && (fail != 0)) {
        *fail = 1;
    }
    
    result = fx[len/2];
    
    free(resids);
    free(fx);
    
    return result;
}

int rloess(double* data, double* dout, int len, int width) {
    int i, n;
    double p[3];
    
    width += width % 2 - 1;
    
    double *fx = (double *)malloc(width*sizeof(double));
    double *resids = (double *)malloc(width*sizeof(double));
    double *din = data;
    
    struct fplyset sp;

    sp.data = din;
    sp.len = width;
    sp.n = 3;
    sp.fx = fx;
    sp.k = p;
    sp.resids = resids;
    
    n = rfitpoly(&sp);
    if (!n) return -3;
    
    for (i = 0; i < width/2 + 1; ++i)
        dout[i] = fx[i];
    
    for (i = width/2 + 1; i < len - width/2; ++i) {
        sp.data++;
        n = rfitpoly(&sp);
        if (!n) return -3;
        dout[i] = sp.k[0];//;//fx[width/2];
    }

    for (i = width/2 + 1; i < width; ++i)
        dout[len + i - width] = fx[i];
    
//    printf("sse p = %f\n", sp.sse);
//    printf("mean p = %f\n", sp.mean);
//    printf("rmse p = %f\n", sp.rmse);
//    printf("mad p = %f\n", mad(sp.resids, sp.mean, sp.len));
//    vcr_save(fx, width, "fx.txt");
    
    free(resids);
    free(fx);
    return 0;
}
