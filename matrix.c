//---------------------------------------------------------------------------


#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


#define BOOL int
#define TRUE 1
#define FALSE 0
//---------------------------------------------------------------------------

//выделение памяти под матрицу
double **mtx_create(int h, int w) {
    int i;
  double **res;
  res = (double **)malloc(h*sizeof(double *));
  for (i = 0; i < h; i++) 
      res[i] = (double *)malloc(w*sizeof(double));
  return res;
}
//---------------------------------------------------------------------------

//double *mtxCreate(int h) {
//  double *res;
//  res = malloc(h*sizeof(double));
//  return res;
//}
//---------------------------------------------------------------------------
//удаление матрицы
void mtx_free(double **m, int h) {
    int i;
    for (i = 0; i < h; i++) 
        free(m[i]);
    free(m);
}
//---------------------------------------------------------------------------
//копирование матрицы
//void mtxCpy(double **m, double **mr, int h, int w) {
//    int i, j;
//    for (i = 0; i < h; i++)
//        for (j = 0; j < w; j++) mr[i][j] = m[i][j];
//}

//void mtxCpy(double *m, double *mr, int h) {
//    int i;
//    for (i = 0; i < h; i++) 
//        mr[i] = m[i];
//}

//---------------------------------------------------------------------------


//умножение матриц

void mtx_mul(double **m1, double **m2, double **mr, int h1, int w1h2, int w2) {
    int i, j, k;
    for (j = 0; j < w2; j++)
        for (i = 0; i < h1; i++) {
            mr[i][j] = 0;
            for (k = 0; k < w1h2; k++) mr[i][j] += m1[i][k]*m2[k][j];
        }
}

//void mtxMul(double **m, double c, double **mr, int h, int w) {
//    int i, j;
//  for (i = 0; i < h; i++)
//    for (j = 0; j < w; j++) mr[i][j] = m[i][j]*c;
//}
//
//void mtxMul(double *v, double c, double *vr, int h) {
//    int i;
//  for (i = 0; i < h; i++) vr[i] = v[i]*c;
//}
//

void mvr_mul(double **m, double *v, double *vr, int h, int w) {
    int i, k;
    for (i = 0; i < h; i++) {
        vr[i] = 0.0;
        for (k = 0; k < w; k++) vr[i] += m[i][k]*v[k];
    }
}

//---------------------------------------------------------------------------
//сложение матриц
void mtx_add(double **m1, double **m2, double **mr, int h, int w) {
    int i, j;
    for (i = 0; i < h; i++)
        for (j = 0; j < w; j++) 
            mr[i][j] = m1[i][j] + m2[i][j];
}

void vcr_add(double *m1, double *m2, double *mr, int h) {
    int i;
    for (i = 0; i < h; i++) 
        mr[i] = m1[i] + m2[i];
}

void vcr_fill(double *m, double value, int h) {
    int i;
    for (i = 0; i < h; ++i)
        m[i] = value;
}

//---------------------------------------------------------------------------
//вычитание матриц
//void mtxDlt(double **m1, double **m2, double **mr, int h, int w) {
//    int i, j;
//    for (i = 0; i < h; i++)
//        for (j = 0; j < w; j++) mr[i][j] = m1[i][j] - m2[i][j];
//}
//
void vcr_dlt(double *m1, double *m2, double *mr, int h) {
    int i;
    for (i = 0; i < h; i++) 
        mr[i] = m1[i] - m2[i];
}

void vcr_cpy(double *m, double *mr, int h) {
    memcpy(mr, m, h*sizeof(double));
}

void vcr_mul(double *m1, double *m2, double *mr, int h) {
    int i;
    for (i = 0; i < h; i++) 
        mr[i] = m1[i]*m2[i];
}

void vcr_sqr(double *v, double **m, int h) {
    int i, j;
    for (j = 0; j < h; ++j)
        for (i = 0; i < h; ++i) {
            m[i][j] = v[i]*v[j];
        }
}

int vcr_maxi(double *v, int len) {
    int i;
    int m = 0;
    
    for (i = 1; i < len; ++i) {
        if (v[m] < v[i]) {
            m = i;
        }
    }
    
    return m;
}

void vcr_abs(double *v, double *vr, int len) {
    int i;
    
    for (i = 1; i < len; ++i) {
        vr[i] = fabs(v[i]);
    }
}

void vrc_mul(double *m, double c, double *mr, int h) {
    int i;
    for (i = 0; i < h; i++) 
        mr[i] = m[i]*c;
}

void vrm_mul(double *v, double **m, double **mr, int h, int w) {
    int i, j;
    for (i = 0; i < h; i++)
        for (j = 0; j < w; j++) 
            mr[j][i] = m[i][j]*v[i];
}

void vrm_tmul(double *v, double **m, double **mr, int h, int w) {
    int i, j;
    for (i = 0; i < h; i++)
        for (j = 0; j < w; j++) 
            mr[i][j] = m[i][j]*v[j];
}

//---------------------------------------------------------------------------
void mtx_zfill(double **m, int h, int w) {
    int i;
    for (i = 0; i < h; ++i)
        memset(m[i], 0, w*sizeof(double));
}

//квадрат транспонированной матрицы
void mtx_sqr(double **m, double **mr, int h, int w) {
    int i, j, k;
    double **tmp = mtx_create(h, h);
    mtx_zfill(mr, h, h);
  
    for (j = 0; j < w; j++) {
        for (i = 0; i < h; i++)
            for (k = 0; k < h; k++) 
                tmp[k][i] = m[k][j]*m[i][j];
        
        mtx_add(mr, tmp, mr, h, h);
    }
  mtx_free(tmp, h);
}

void mtx_jsqr(double **m, double **mr, int h, int w) {
    int i, j, k;
    double **tmp = mtx_create(h, h);
    mtx_zfill(mr, h, h);
  
    for (j = 0; j < w; j++) {
        for (i = 0; i < h; i++)
            for (k = 0; k < h; k++) 
                tmp[k][i] = m[k][j]*m[i][j];
        
        mtx_add(mr, tmp, mr, h, h);
    }
  mtx_free(tmp, h);
}

//---------------------------------------------------------------------------
//транспонирование
void mtx_trp(double **m, double **mr, int h, int w) {
    int i, j;
    for (j = 0; j < w; j++)
        for (i = 0; i < h; i++) 
            mr[j][i] = m[i][j];
}

//---------------------------------------------------------------------------
//разложение Холецкого

int mtx_hol(double **m, double **mr, int h) {
    int l, j, k;
    double tmp;

    mtx_zfill(mr, h, h);

    for (k = 0; k < h; k++)
        for (j = k; j < h; j++) {
            tmp = 0.0;
            for (l = 0; l < k; l++) tmp += mr[k][l]*mr[j][l];
            mr[j][k] = m[j][k] - tmp;
            if (mr[j][k] <= 0 && j == k) return FALSE;
            if (j == k) mr[j][k] = sqrt(mr[j][k]);
            else mr[j][k] = mr[j][k]/mr[k][k];
        }
    return TRUE;
}

int mtx_ihol(double **m, double **mr, int h) {
    int l, j, k;
    double tmp;

    mtx_zfill(mr, h, h);

    for (k = 0; k < h; k++)
        for (j = k; j < h; j++) {
            tmp = 0.0;
            for (l = 0; l < k; l++) tmp += mr[k][l]*mr[j][l];
            mr[j][k] = m[j][k] - tmp;
            if (mr[j][k] <= 0 && j == k) {
//                return FALSE;
                mr[j][k] = 1;
                puts("incomplete cholesky\n");
            }
            if (j == k) mr[j][k] = sqrt(mr[j][k]);
            else mr[j][k] = mr[j][k]/mr[k][k];
        }
    return TRUE;
}

//---------------------------------------------------------------------------
//обращение матрицы через Холецкого
int mtx_hinv(double **m, double **mr, int h) {
    int i, j, k, l;
    double **mh = mtx_create(h, h);
  
    if (!mtx_ihol(m, mh, h)) {
          mtx_free(mh, h);
          return FALSE;
    }
  
    for (i = 0; i < h; i++) {
        for (l = i; l < h; l++) {
            if (l == i) mr[l][h - 1] = 1/mh[l][l];
            else {
                mr[l][h - 1] = 0;
                for (k = i; k < l; k++) mr[l][h - 1] -= mh[l][k]*mr[k][h - 1];
                mr[l][h - 1] = mr[l][h - 1]/mh[l][l];
            }
        }
        for (l = h - 1; l >= i; l--) {
            if (i == h - 1) mr[l][i] = mr[l][h-1]/mh[l][l];
            else {
                mr[l][i] = mr[l][h-1];
                for (k = h - 1; k >= l + 1; k--) mr[l][i] -= mh[k][l]*mr[k][i];
                mr[l][i] = mr[l][i]/mh[l][l];
            }
        }
    }
    
    for (j = 0; j < h; j++)
        for (i = 0; i < j; i++) mr[i][j] = mr[j][i];
    
    mtx_free(mh, h);
    return TRUE;
}

//---------------------------------------------------------------------------
//сохранение матрицы в текстовом виде в файл
void mtx_save(double **m, int h, int w, char *file) {
    int i, j;
    
    FILE *f = fopen(file, "wt");

    if (f != 0)
        for (i = 0; i < h; i++) {
            for (j = 0; j < w; j++) 
                fprintf(f, "%12.12f ", m[i][j]);
            fprintf(f, "\n");
        }
    
    fclose(f);
}

void mtx_print(double **m, int h, int w) {
    int i, j;
    
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) 
            printf("%16.16lf ", m[i][j]);
        printf("\n");
    }
    puts("");
}

void vcr_print(double *m, int h) {
    int i;
    
    for (i = 0; i < h; i++) 
        printf("%16.16g ", m[i]);
    puts("\n");
}

void vcr_save(double *m, int h, char *file) {
    int i;
    
    FILE *f = fopen(file, "wt");
    
    if (f != 0)
        for (i = 0; i < h; i++) {
            fprintf(f, "%16.16lf\n", m[i]);
        }
    
    fclose(f);
}

int mtx_load(double **buffer, int h, int w, char *file) {
    int i = 0;
    int j = 0;
    int n;
    
    FILE *f = fopen(file, "rt");
    
    if (f != 0)
        for (j = 0; j < h; ++j) {
            for (i = 0; i < w; ++i) {
                n = fscanf(f, "%lf", &buffer[j][i]);
                if (n < 1) break;
            }
            if (n < 1) break;
        }
    
    fclose(f);
    return j*i;
}

int vcr_load(double *buffer, int len, char *file) {
    int i = 0;
    int n;
    FILE *f = fopen(file, "rt");
    
    if (f != 0)
        for (i = 0; i < len; ++i) {
            n = fscanf(f, "%lf", &buffer[i]);
            if (n < 1) break;
        }
    
    fclose(f);
    return i;
}

int vcr_load_len(char *file) {
    int i = -1;
    int n;
    double tmp;
    
    FILE *f = fopen(file, "rt");
    
    if (f != 0)
        do {
            n = fscanf(f, "%lf", &tmp);
            i++;
        } while (n == 1);
    
    fclose(f);
    return i;
}
