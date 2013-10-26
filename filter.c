
#include "filter.h"
#include "statistic.h"

void statistic_f(struct firset *s);
void statistic_s(struct sthset *s);


void statistic_f(struct firset *s) {
    s->mean = mean(s->result, s->len);
    s->rmse = stddev(s->result, s->mean, s->len);
    s->sse = sse(s->result, s->len);
}

void statistic_s(struct sthset *s) {
    s->mean = mean(s->result, s->len);
    s->rmse = stddev(s->result, s->mean, s->len);
    s->sse = sse(s->result, s->len);
}

void fir(struct firset *s) {
    int i, j, k;
    
    for (j = 0; j < s->len; ++j) {
        s->result[j] = 0;
        k = (j - s->order < 0) ? j : s->order - 1;
        for (i = (j - s->order < 0) ? 0 : j - s->order + 1; i < j + 1; ++i) {
            s->result[j] += s->k[k]*s->data[i];
            k--;
        }
    }
    statistic_f(s);
}

void smooth(struct sthset *s) {
    int i;
    double value = s->data[0];
    s->order = s->order + s->order%2 - 1;
    s->result[0] = value;
    
    for(i = 1; i < s->order/2 + 1; ++i) {
        value += s->data[2*i - 1] + s->data[2*i];
        s->result[i] = value/(2*i + 1);
    }
    
    for (i = s->order; i < s->len; ++i) {
        value += s->data[i] - s->data[i - s->order];
        s->result[i - s->order/2] = value/s->order;
    }
    
    for(i = s->len - s->order/2; i < s->len; ++i) {
        value -= s->data[2*i - s->len - 1] + s->data[2*i - s->len];
        s->result[i] = value/(2*(s->len - i) - 1);
    }
    statistic_s(s);
}
