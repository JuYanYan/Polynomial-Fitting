#include "ployfitting.h"
#include <stdlib.h>
// 变量交换值
void swapValue(RealNumber *a, RealNumber *b)
{
    *a = *a + *b;
    *b = *a - *b;
    *a = *a - *b;
}
// 载入数据到矩阵里面
void LoadData_Into_Matrix(RealNumber *m, const RealNumber *v, int anNum, int inNum)
{
    int        i, j, k;
    RealNumber t;
    for (i = 0; i < anNum; i++)
    {
        for (j = 0; j < inNum; j++)
        {
            t = 1.0;
            for (k = 0; k < i; k++)
                t = t * v[j];
            *BUFF2D(m, i, j, inNum) = t;
        }
    }
}
// 计算正规方程组的系数矩阵
void calcNormalEqCof(RealNumber *a, RealNumber *b, const RealNumber *y, int anNum, int inNum)
{
    int         i, j, k;
    RealNumber  t;
    for (i = 0; i < anNum; i++)
    {
        for (j = 0; j < anNum; j++)
        {
            t = 0.0;
            for (k = 0; k < inNum; k++)
                t = t + (*BUFF2D(a, i, k, inNum) * (*BUFF2D(a, j, k, inNum)));
            *BUFF2D(b, i, j, anNum + 1) = t;
        }
        t = 0.0;
        for (k = 0; k < inNum; k++)
        {
            t = t + (y[k] * (*BUFF2D(a, i, k, inNum)));
            *BUFF2D(b, i, anNum, anNum + 1) = t;                    // 矩阵增广
        }
    }
}
// 列主元LU分解
void DirectLU(RealNumber *a, RealNumber *x, int n)
{
    int          i, r, k, j;
    int          n2;
    RealNumber  *s, *t;
    RealNumber   m;
    n2 = n + 1;
    s = (RealNumber*)malloc(sizeof(RealNumber) * n);
    t = (RealNumber*)malloc(sizeof(RealNumber) * n);
    for (r = 0; r < n; r++)
    {
        m = 0;
        j = r;
        for (i = r; i < n; i++)                 // 选主元
        {
            s[i] = *BUFF2D(a, i, r, n2);
            for (k = 0; k < r; k++)
                s[i] = s[i] - ((*BUFF2D(a, i, k, n2)) * (*BUFF2D(a, k, r, n2)));
            s[i] = s[i] > 0 ? s[i] : -s[i];     // 绝对值
            if (s[i] > m)
            {
                j = i;
                m = s[i];
            }
        }
        if (j != r)                             // 主元的位置所在行j不是r的话，就调换
        {
            for (i = 0; i < n2; i++)
                swapValue(a + r * n2 + i, a + j * n2 + i);
        }
        for (i = r; i < n2; i++)                // 计算第r行的元素
        {
            for (k = 0; k < r; k++)
            {
                *BUFF2D(a, r, i, n2) = *BUFF2D(a, r, i, n2) - ((*BUFF2D(a, r, k, n2)) * (*BUFF2D(a, k, i, n2)));
            }
        }
        for (i = r + 1; i < n2; i++)            // 计算第r列的元素
        {
            for (k = 0; k < r; k++)
            {
                *BUFF2D(a, i, r, n2) = *BUFF2D(a, i, r, n2) - (*BUFF2D(a, i, k, n2)) * (*BUFF2D(a, k, r, n2));
            }
            *BUFF2D(a, i, r, n2) = *BUFF2D(a, i, r, n2) / (*BUFF2D(a, r, r, n2));
        }
    }
    for (i = 0; i < n; i++)
    {
        t[i] = *BUFF2D(a, i, n, n2);
    }
    for (i = n - 1; i >= 0; i--)                // 回代法计算出最后的解
    {
        for (r = n - 1; r > i; r--)
        {
            t[i] = t[i] - *BUFF2D(a, i, r, n2) * x[r];
        }
        x[i] = t[i] / *BUFF2D(a, i, i, n2);
    }
    free(s);
    free(t);
}
