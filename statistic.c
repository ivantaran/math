
#include "statistic.h"
#include <math.h>

double scale_mean(int len) {
    return ((double)len - 1)/2;
}

double scale_std(double mean, int len) {
    int i;
    double std = 0;
    double dq = 0;
    
    for (i = 0; i < len; ++i) {
        dq = i - mean;
        dq *= dq;
        std += dq;
    }
    std = sqrt(std/(len - 1));

    return std;
}

double mean(double *data, int len) {
    int i;
    double value = 0;
    
    for (i = 0; i < len; ++i) {
        value += data[i];
    }
    
    return value/len;
}

double stddev(double *data, double mean, int len) {
    int i;
    double std = 0;
    double dq = 0;
    
    for (i = 0; i < len; ++i) {
        dq = data[i] - mean;
        dq *= dq;
        std += dq;
    }
    std = sqrt(std/(len - 1));

    return std;
}

double mad(double *data, double mean, int len) {
    int i;
    double md = 0;
    
    for (i = 0; i < len; ++i) {
        md += fabs(data[i] - mean);
    }
    md /= len;

    return md;
}

double sse(double *data, int len) {
    int i;
    double value = 0;
    
    for (i = 0; i < len; ++i) {
        value += data[i]*data[i];
    }
    
    return value;
}