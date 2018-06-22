# Polynomial-Fitting
A Polynomial fitting which uses the least square method.And it is believable.
# How to use
You can use it like:
```
bool PloyFitting(RealNumber *anOut, const RealNumber *x, const RealNumber *y, int xNum, int outNum)
{
    RealNumber  *a, *b;
    a = new RealNumber[outNum*xNum];
    b = new RealNumber[(outNum + 1)*(outNum + 1)];

    LoadData_Into_Matrix(a, x, outNum, xNum);

    calcNormalEqCof(a, b, y, outNum, xNum);

    DirectLU(b, anOut, outNum);

    delete[] a;
    delete[] b;
    return true;
}
```
Or.....You can do it that you want to.
