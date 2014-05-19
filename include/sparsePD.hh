#ifndef SPARSE_HH_
#define SPARSE_HH_

SpMat genSMPD(int n);
MatrixXd gen_index(int n);
VectorXd gen_b(MatrixXd B, MatrixXd F, int n);
MatrixXd vectorToMatrix(MatrixXd B, VectorXd u, int n);
MatrixXd exPoissonLDLT();
#endif