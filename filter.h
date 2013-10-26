/* 
 * File:   fir.h
 * Author: Ivan
 *
 * Created on 21 Декабрь 2012 г., 12:23
 */

#ifndef FILTER_H
#define	FILTER_H

#ifdef	__cplusplus
extern "C" {
#endif


    struct firset {
        int order;
        int len;
        double mean;
        double rmse;
        double sse;
        double *k;
        double *result;
        double *data;
    };


    struct sthset {
        int order;
        int len;
        double mean;
        double rmse;
        double sse;
        double *result;
        double *data;
    };

    void fir(struct firset *s);
    void smooth(struct sthset *s);
    
#ifdef	__cplusplus
}
#endif

#endif	/* FILTER_H */

