#pragma once
typedef double RealNumber;
void             LoadData_Into_Matrix(RealNumber *m, const RealNumber *v, int anNum, int inNum);
void             calcNormalEqCof(RealNumber *a, RealNumber *b, const RealNumber *y, int anNum, int inNum);
void             DirectLU(RealNumber *a, RealNumber *x, int n);
