/* 
 * File:   statistic.h
 * Author: yanco
 *
 * Created on 23 Декабрь 2012 г., 18:16
 */

#ifndef STATISTIC_H
#define	STATISTIC_H

#ifdef	__cplusplus
extern "C" {
#endif


    double scale_mean(int len);
    double scale_std(double mean, int len);
    double mean(double *data, int len);
    double stddev(double *data, double mean, int len);
    double sse(double *data, int len);
    double mad(double *data, double mean, int len);
    
#ifdef	__cplusplus
}
#endif

#endif	/* STATISTIC_H */

