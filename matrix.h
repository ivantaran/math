/* 
 * File:   matrix.h
 * Author: yanco
 *
 * Created on 16 Декабрь 2012 г., 20:16
 */

#ifndef MATRIX_H
#define	MATRIX_H

#ifdef	__cplusplus
extern "C" {
#endif


    double **mtx_create(int h, int w);
    void mtx_free(double **m, int h);
    void mtx_zfill(double **m, int h, int w);
    void mtx_add(double **m1, double **m2, double **mr, int h, int w);
    void mtx_mul(double **m1, double **m2, double **mr, int h1, int w1h2, int w2);
    void vcr_fill(double *m, double value, int h);
    void vcr_add(double *m1, double *m2, double *mr, int h);
    void vcr_dlt(double *m1, double *m2, double *mr, int h);
    void vcr_mul(double *m1, double *m2, double *mr, int h);
    void vcr_sqr(double *v, double **m, int h);
    void vcr_abs(double *v, double *vr, int len);
    int vcr_maxi(double *v, int len);
    void vrc_mul(double *m, double c, double *mr, int h);
    void vcr_cpy(double *m, double *mr, int h);
    void vrm_tmul(double *v, double **m, double **mr, int h, int w);
    void vrm_mul(double *v, double **m, double **mr, int h, int w);
//    void mtxDlt(double *m1, double *m2, double *mr, int h);
    void mtx_sqr(double **m, double **mr, int h, int w);
    int mtx_ihol(double **m, double **mr, int h);
    int mtx_hol(double **m, double **mr, int h);
    void mtx_jsqr(double **m, double **mr, int h, int w);
    void mtx_trp(double **m, double **mr, int h, int w);
    //void mtxTrp(double **m, double **mr, int h, int w);
    int mtx_hinv(double **m, double **mr, int h);
    void mtx_save(double **m, int h, int w, char *file);
    void mtx_print(double **m, int h, int w);
    void vcr_print(double *m, int h);
    void vcr_save(double *m, int h, char *file);
    int vcr_load_len(char *file);
    int vcr_load(double *buffer, int len, char *file);
    int mtx_load(double **buffer, int h, int w, char *file);

#ifdef	__cplusplus
}
#endif

#endif	/* MATRIX_H */

