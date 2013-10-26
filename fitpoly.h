/* 
 * File:   fitpoly.h
 * Author: Ivan
 *
 * Created on 21 Декабрь 2012 г., 10:46
 */

#ifndef FITPOLY_H
#define	FITPOLY_H

#ifdef	__cplusplus
extern "C" {
#endif

    struct fplyset {
        int n;
        int len;
        double smean;
        double sstd;
        double sse;
        double rmse;
        double mean;
        double *k;
        double *data;
        double *resids;
        double *fx;
    };

    int fitpoly(struct fplyset *s);
    int rfitpoly(struct fplyset *s);
    double rloess_point(double* data, int len, int* fail);
    int rloess(double* data, double* dout, int len, int width);

#ifdef	__cplusplus
}
#endif

#endif	/* FITPOLY_H */

