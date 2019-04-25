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
            *(m + i*inNum + j) = t;
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
                t = t + (*(a + i*inNum + k)*(*(a + j*inNum + k)));
            *(b + i*(anNum + 1) + j) = t;
        }
        t = 0.0;
        for (k = 0; k < inNum; k++)
        {
            t = t + (y[k] * (*(a + i*inNum + k)));
            *(b + i*(anNum + 1) + anNum) = t;
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
            s[i] = *(a + i*n2 + r);
            for (k = 0; k < r; k++)
                s[i] = s[i] - ((*(a + i*n2 + k)) * (*(a + k*n2 + r)));
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
                swapValue(a + r*n2 + i, a + j*n2 + i);
        }
        for (i = r; i < n2; i++)                // 计算第r行的元素
        {
            for (k = 0; k < r; k++)
            {
                *(a + r*n2 + i) = *(a + r*n2 + i) - ((*(a + r*n2 + k)) * (*(a + k*n2 + i)));
            }
        }
        for (i = r + 1; i < n2; i++)            // 计算第r列的元素
        {
            for (k = 0; k < r; k++)
            {
                *(a + i*n2 + r) = *(a + i*n2 + r) - (*(a + i*n2 + k)) * (*(a + k*n2 + r));
            }
            *(a + i*n2 + r) = *(a + i*n2 + r) / (*(a + r*n2 + r));
        }
    }
    for (i = 0; i < n; i++)
    {
        t[i] = *(a + i*n2 + n);
    }
    for (i = n - 1; i >= 0; i--)                // 回代法计算出最后的解
    {
        for (r = n - 1; r > i; r--)
        {
            t[i] = t[i] - *(a + i*n2 + r) * x[r];
        }
        x[i] = t[i] / *(a + i*n2 + i);
    }
    free(s);
    free(t);
}

