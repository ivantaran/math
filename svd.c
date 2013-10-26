
#include "svd.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "matrix.h"

#define FALSE 0
#define TRUE  1

void mvr_gcl(double **m, double *v, int h, int index);
void mvr_pcl(double **m, double *v, int h, int index);
void mvr_grw(double **m, double *v, int w, int index);
void mvr_prw(double **m, double *v, int w, int index);
void vcr_hsd(double *v, double *up, double *b, int h, int l);
void vcr_hst(double *v, double up, double b, double *c, int h, int l);
void gvd(double v1, double v2, double *c, double *s);
void gva(double *v1, double *v2, double *c, double *s);
void gvt(double *z1, double *z2, double c, double s);
void mtx_2dtod1(double **m, double *d, double *e, int w, int k);
void mtx_2dtod2(double **b, double *d, double *e, int wb, int k, int l);
void mtx_2dtod3(double **m, double **b, double *d, double *e, int w, int wb, int k, int l);


/* СЃС‚РѕР»Р±РµС† РјР°С‚СЂРёС†С‹ РІ РІРµРєС‚РѕСЂ */
void mvr_gcl(double **m, double *v, int h, int index) {
    int i;
    
    for (i = 0; i < h; ++i)
        v[i] = m[i][index];
}

/* РІРµРєС‚РѕСЂ РІ СЃС‚РѕР»Р±РµС† РјР°С‚СЂРёС†С‹ */
void mvr_pcl(double **m, double *v, int h, int index) {
    int i;
    
    for (i = 0; i < h; ++i)
        m[i][index] = v[i];
}

/* СЃС‚СЂРѕРєСѓ РјР°С‚СЂРёС†С‹ РІ РІРµРєС‚РѕСЂ */
void mvr_grw(double **m, double *v, int w, int index) {
    int i;
    
    for (i = 0; i < w; ++i)
        v[i] = m[index][i];
}

/* РІРµРєС‚РѕСЂ РІ СЃС‚СЂРѕРєСѓ РјР°С‚СЂРёС†С‹ */
void mvr_prw(double **m, double *v, int w, int index) {
    int i;
    
    for (i = 0; i < w; ++i)
        m[index][i] = v[i];
}

/* РћРїСЂРµРґРµР»РµРЅРёРµ РїР°СЂР°РјРµС‚СЂРѕРІ РїСЂРµРѕР±СЂР°Р·РѕРІР°РЅРёСЏ РҐР°СѓСЃС…РѕР»РґРµСЂР° */
void vcr_hsd(double *v, double *up, double *b, int h, int l) {
    int i;
    double c1 = 0; 
    double sd = 0;
    double c = fabs(v[l]);
    double p = 0;
    
    for (i = l + 1; i < h; ++i) {
        c = fmax(fabs(v[i]), c);
    }
    
    if (c > 0) {
        c1 = 1/c;
        for (i = l; i < h; ++i) {
            sd += v[i]*c1*v[i]*c1;
        }
        p = sd;
        p = c*sqrt(fabs(p));
        if (v[l] > 0) p = -p;
        *up = v[l] - p;
        *b = 1/(p*(*up));
    }
}

/* РџСЂРµРѕР±СЂР°Р·РѕРІР°РЅРёРµ РҐР°СѓСЃС…РѕР»РґРµСЂР° РІРµРєС‚РѕСЂР° */
void vcr_hst(double *v, double up, double b, double *c, int h, int l) {
    int i;
    double dup = up;
    double s = c[l]*dup;
    
    for (i = l + 1; i < h; ++i) {
        s += c[i]*v[i];
    }
    
    s *= b;
    c[l] += s*dup;
    
    for (i = l + 1; i < h; ++i) {
        c[i] += s*v[i];
    }
}

void mtx_to2d(double **m, double **b, double *d, double *e, int h, int w, int wb) {
    int i, j, k;
    double up, bb;
    
    double *v = malloc(w*sizeof(double));
    double *s = malloc(w*sizeof(double));
    double *ups = malloc(w*sizeof(double));
    double *bbs = malloc(w*sizeof(double));
    
    memset(ups, 0, w*sizeof(double));
    memset(bbs, 0, w*sizeof(double));
    memset(v, 0, w*sizeof(double));
    memset(s, 0, w*sizeof(double));
    
    for (i = 0; i < w; ++i) {
        if ((i < w - 1) || (h  > w)) {   // !
            mvr_gcl(m, v, h, i);
            vcr_hsd(v, &up, &bb, h, i);
            for (j = i; j < w; ++j) {
                mvr_gcl(m, s, h, j);
                vcr_hst(v, up, bb, s, h, i);
                mvr_pcl(m, s, h, j);
            }
            
            for (k = 0; k < wb; ++k) {
                mvr_gcl(b, s, h, k);
                vcr_hst(v, up, bb, s, h, i);
                mvr_pcl(b, s, h, k);
            }
        }
        
        if (i < w - 2) {
            mvr_grw(m, v, w, i);
            vcr_hsd(v, &up, &bb, w, i + 1);
            ups[i] = up;
            bbs[i] = bb;
            for (j = i; j < h; ++j) {
                mvr_grw(m, s, w, j);
                vcr_hst(v, up, bb, s, w, i + 1);
                if (j == i)
                    for (k = i + 2; k < w; ++k)
                        s[k] = v[k];
                mvr_prw(m, s, w, j);
            }
        }
    }
    if (w > 1)
        for (i = 1; i < w; ++i) {
            d[i] = m[i][i];
            e[i] = m[i - 1][i];
        }
    
    d[0] = m[0][0];
    e[0] = 0;
    for (i = w - 1; i >= 0; --i) {
        if (i < w - 1) mvr_grw(m, v, w, i);
        for (k = 0; k < w; ++k)
            m[i][k] = 0;
        m[i][i] = 1;
        if (i < w - 2)
            for (k = i; k < w; ++k) {
                mvr_gcl(m, s, h, k);
                vcr_hst(v, ups[i], bbs[i], s, w, i + 1);
                mvr_pcl(m, s, h, k);
            }
    }
    
    free(v);
    free(s);
    free(ups);
    free(bbs);
}

void gvd(double v1, double v2, double *c, double *s) {
    double a1 = fabs(v1);
    double a2 = fabs(v2);
    double w, q;
    
    if (a1 > a2) {
        w = v2/v1;
        q = sqrt(1 + w*w);
        *c = 1/q;
        if (v1 < 0) *c = -(*c);
        *s = (*c)*w;
    }
    else
        if (v2 != 0) {
            w = v1/v2;
            q = sqrt(1 + w*w);
            *s = 1/q;
            if (v2 < 0) *s = -(*s);
            *c = (*s)*w;
        }
        else {
            *c = 1;
            *s = 0;
        }
}

void gva(double *v1, double *v2, double *c, double *s) {
    double a1 = fabs(*v1);
    double a2 = fabs(*v2);
    double w, q;
    
    if (a1 > a2) {
        w = (*v2)/(*v1);
        q = sqrt(1 + w*w);
        *c = 1/q;
        if (*v1 < 0) *c = -(*c);
        *s = (*c)*w;
        *v1 = a1*q;
        *v2 = 0;
    }
    else
        if (*v2 != 0) {
            w = (*v1)/(*v2);
            q = sqrt(1 + w*w);
            *s = 1/q;
            if (*v2 < 0) *s = -(*s);
            *c = (*s)*w;
            *v1 = a2*q;
            *v2 = 0;
        }
        else {
            *c = 1;
            *s = 0;
        }
}

void gvt(double *z1, double *z2, double c, double s) {
    double w = (*z1)*c + (*z2)*s;
    *z2 = -(*z1)*s + (*z2)*c;
    *z1 = w;
}

void mtx_2dtod1(double **m, double *d, double *e, int w, int k) {
    int i, j;
    double cs, sn, h;
    
    for (i = k - 1; k >= 0; --k) {
        if (i == k - 1)
            gva(&d[i], &e[i + 1], &cs, &sn);
        else
            gva(&d[i], &h, &cs, &sn);
        
        if (i > 0) {
            h = 0;
            gvt(&e[i], &h, cs, sn);
        }
        
        for (j = 0; j < w; ++j)
            gvt(&m[j][i], &m[j][k], cs, sn);
    }
}

void mtx_2dtod2(double **b, double *d, double *e, int wb, int k, int l) {
    int i, j;
    double cs, sn, h;

    for (i = l; i < k + 1; ++i) {
        if (i == l)
            gva(&d[i], &e[i], &cs, &sn);
        else
            gva(&d[i], &h, &cs, &sn);
        
        if (i < k) {
            h = 0;
            gvt(&e[i + 1], &h, cs, sn);
        }
        
        for (j = 0; j < wb; ++j)
            gvt(&cs, &sn, b[i][j], b[l - 1][j]);
    }
}

void mtx_2dtod3(double **m, double **b, double *d, double *e, int w, int wb, int k, int l) {
    int i, j;
    double cs, sn, h, f, g, t;
    
    f = ((d[k - 1] - d[k])*(d[k - 1] + d[k]) +
         (e[k - 1] - e[k])*(e[k - 1] + e[k]))/(2*e[k]*d[k - 1]);
    
    if (fabs(f) > 1e10) 
        g = fabs(f);
    else
        g = sqrt(1 + f*f);
    
    if (f >= 0)
        t = f + g;
    else
        t = f - g;
    f = ((d[l] - d[k])*(d[l] + d[k]) + e[k]*(d[k - 1]/t - e[k]))/d[l];
    
    for (i = l; i < k; ++i) {
        if (i == l)
            gvd(f, e[i + 1], &cs, &sn);
        else
            gva(&e[i], &h, &cs, &sn);
        gvt(&d[i], &e[i + 1], cs, sn);
        h = 0;
        gvt(&h, &d[i + 1], cs, sn);
        
        for (j = 0; j < w; ++j)
            gvt(&m[j][i], &m[j][i + 1], cs, sn);
        
        gva(&d[i], &h, &cs, &sn);
        gvt(&e[i + 1], &d[i + 1], cs, sn);
        
        if (i < k - 1) {
            h = 0;
            gvt(&h, &e[i + 2], cs, sn);
        }
        
        for (j = 0; j < wb; ++j)
            gvt(&b[i][j], &b[i + 1][j], cs, sn);
    }
}

int mtx_2dtod(double **m, double **b, double *d, double *e, int h, int w, int wb) {
    int i, j, k, l, ll;
    double bmx;
    int niter = 0;
    int niterm = 10*w;
    int elzero;
    int result;
    int ok = TRUE;
    
    bmx = d[0];

    if (w > 0) {
        for (i = 1; i < w; ++i)
            bmx = fmax(fabs(d[i]) + fabs(e[i]), bmx);
    }
    
    for (k = w - 1; k >= 0; --k) {
//        niterm = 10*w;  //?
        do {
            result = FALSE;
            if (k != 0) {
                if ((bmx + d[k]) - bmx == 0)
//                if (d[k] == 0)
                    mtx_2dtod1(m, d, e, w, k);

                for (ll = k; ll >= 0; --ll) {
                    l = ll;
                    if (l == 0) {
                        elzero = FALSE;
                        break;
                    }
                    else 
                        if ((bmx - e[l]) - bmx == 0) {
//                        if (e[l] == 0) {
                            elzero = TRUE;
                            break;
                        }
                        else 
                            if ((bmx + d[l - 1]) - bmx == 0)
//                            if (d[l - 1] == 0)
                                elzero = FALSE;
                }

                if ((l > 0) && !elzero)
                    mtx_2dtod2(b, d, e, wb, k, l);

                if (l != k) {
                    result = TRUE;
                    mtx_2dtod3(m, b, d, e, w, wb, k, l);
                    niter++;
                    if (niter > niterm) {
//                        printf("%d\n", niter);
                        ok = FALSE;
                        result = FALSE;
                    }
                }
            }
        } while (result);
        
        if (d[k] < 0) {
            d[k] = -d[k];
            for (j = 0; j < w; ++j)
                m[j][k] = -m[j][k];
        }
    }
    
    return ok;
}

void mtx_dsort(double **m, double **b, double *d, int h, int w, int wb) {
    int i, j, k, index;
    double t;
    
    if (w < 2) return;
    
    index = 1;
    do {
        if (d[index] > d[index - 1]) {
            for (i = 1; i < w; ++i) {
                t = d[i - 1];
                k = i - 1;
                
                for (j = i; j < w; ++j)
                    if (t < d[j]) {
                        t = d[j];
                        k = j;
                    }
                
                if (k != i - 1) {
                    d[k] = d[i - 1];
                    d[i - 1] = t;
                    
                    for (j = 0; j < w; ++j) {
                        t = m[j][k];
                        m[j][k] = m[j][i - 1];
                        m[j][i - 1] = t;
                    }
                    
                    for (j = 0; j < wb; ++j) {
                        t = b[k][j];
                        b[k][j] = b[i - 1][j];
                        b[i - 1][j] = t;
                    }
                }
            }
            index = 1;
        }
        else 
            index++;
    }
    while (index < w);
}

void mtx_svd_analyse(double **m, double **b, double *d, double **x, double *r, int h, int w, int wb, int *kn, double frac) {
    int i, j, k, kk;
    double eps = 1e-15;
    double sinmax = 0;
    double sinmin, s1;
    
    frac = fabs(frac);
    if (frac < eps) frac = eps;
    
    for (i = 0; i < w; ++i)
        sinmax = fmax(sinmax, d[i]);
    
    sinmin = sinmax*frac;
    
    kk = w;
    for (i = 0; i < w; ++i) 
        if (d[i] <= sinmin) {
            kk = i;
            break;
        }
    
    for (i = 0; i < h; ++i)
        if (i < kk) {
            s1 = 1/d[i];
            
            for (j = 0; j < wb; ++j)
                b[i][j] *= s1;
        }
        else
            for (j = 0; j < wb; ++j) {
                if (i == kk)
                    r[j] = b[i][j]*b[i][j];
                else
                    r[j] += b[i][j]*b[i][j];
                
                if (i < w) b[i][j] = 0;
            }
    
    for (i = 0; i < w; ++i)
        for (j = 0; j < wb; ++j) {
            x[i][j] = 0;
            for (k = 0; k < w; ++k)
                x[i][j] += m[i][k]*b[k][j];
        }
    *kn = kk;
}

void mtx_svm(double **m, double **b, double *d, double lam, double *x1, double *x2, int h, int w, int *kn, double frac) {
    int i, k, kk;
    double eps = 1e-15;
    double sinmax = 0;
    double sinmin, g, den1, den2;
    double lam2 = lam*lam;
    double lamp = lam*0.1;
    double lamp2 = lamp*lamp;

    double *p1 = malloc(w*sizeof(double));
    double *p2 = malloc(w*sizeof(double));
    
    frac = fabs(frac);
    if (frac < eps) frac = eps;
    
    for (i = 0; i < w; ++i)
        sinmax = fmax(sinmax, d[i]);
    
    sinmin = sinmax*frac;
    
    kk = w;
    for (i = 0; i < w; ++i) 
        if (d[i] <= sinmin) {
            kk = i;
            break;
        }
    
    for (i = 0; i < h; ++i) {
        g = b[i][0];
        if (i < kk) {
            den1 = 1/(d[i]*d[i] + lam2);
            den2 = 1/(d[i]*d[i] + lamp2);
            p1[i] = g*d[i]*den1;
            p2[i] = g*d[i]*den2;
        }
        else {
            if (i < w) {
                p1[i] = 0;
                p2[i] = 0;
            }
        }
    }
    
    for (i = 0; i < w; ++i) {
        x1[i] = 0;
        x2[i] = 0;
        
        for (k = 0; k < w; ++k) {
            x1[i] += m[i][k]*p1[k];
            x2[i] += m[i][k]*p2[k];
        }
    }
    
    *kn = kk;
    
    free(p1);
    free(p2);
}

int mtx_mar(double **m, double **b, double lam, double *x1, double *x2, int h, int w, int *kn, double frac) {
    int result;
   
    double *d = malloc(w*sizeof(double));
    double *e = malloc(w*sizeof(double));


    mtx_to2d(m, b, d, e, h, w, 1);
    
    result = mtx_2dtod(m, b, d, e, h, w, 1);
    
    mtx_dsort(m, b, d, h, w, 1);

//    vcr_print(d, w);

    mtx_svm(m, b, d, lam, x1, x2, h, w, kn, frac);
    
    free(d);
    free(e);
    
    return result;
}

int mtx_svd_solver(double **m, double **b, double **x, double *r, int h, int w, int wb, int *kn, double frac) {
    int result;
   
    double *d = malloc(w*sizeof(double));
    double *e = malloc(w*sizeof(double));

    mtx_to2d(m, b, d, e, h, w, wb);
    
    result = mtx_2dtod(m, b, d, e, h, w, wb);
    
//    mtx_print(m, h, w);
//    mtx_print(b, h, wb);
//    vcr_print(d, w);
//    vcr_print(e, w);
    
    mtx_dsort(m, b, d, h, w, wb);
    
//    mtx_print(m, h, w);
//    mtx_print(b, h, wb);
//    vcr_print(d, w);
//    vcr_print(e, w);
    
    mtx_svd_analyse(m, b, d, x, r, h, w, wb, kn, frac);
    
    free(d);
    free(e);
    
    return result;
}
