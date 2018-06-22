#pragma once
typedef double RealNumber;
void             LoadData_Into_Matrix(RealNumber *m, const RealNumber *v, unsigned int anNum, unsigned int inNum);
void             calcNormalEqCof(RealNumber *a, RealNumber *b, const RealNumber *y, unsigned int anNum, unsigned int inNum);
void             DirectLU(RealNumber *a, RealNumber *x, unsigned int n);
