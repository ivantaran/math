/* 
 * File:   svd.h
 * Author: Ivan
 *
 * Created on 27 Декабрь 2012 г., 10:29
 */

#ifndef SVD_H
#define	SVD_H

#ifdef	__cplusplus
extern "C" {
#endif

    void mtx_to2d(double **m, double **b, double *d, double *e, int h, int w, int wb);
    int mtx_2dtod(double **m, double **b, double *d, double *e, int h, int w, int wb);
    void mtx_dsort(double **m, double **b, double *d, int h, int w, int wb);
    void mtx_svd_analyse(double **m, double **b, double *d, double **x, double *r, int h, int w, int wb, int *kn, double frac);
    int mtx_svd_solver(double **m, double **b, double **x, double *r, int h, int w, int wb, int *kn, double frac);
    void mtx_svm(double **m, double **b, double *d, double lam, double *x1, double *x2, int h, int w, int *kn, double frac);
    int mtx_mar(double **m, double **b, double lam, double *x1, double *x2, int h, int w, int *kn, double frac);
    
#ifdef	__cplusplus
}
#endif

#endif	/* SVD_H */

