/*
 * edp_06.cpp
 *
 *  Created on: Dec 29, 2012
 *      Author: ataias
 *
 *      Testes sobre salvar matrizes da Eigen em arquivos binários
 */

#ifndef IOSTREAM_H
#define IOSTREAM_H
#include<iostream>
#endif

#ifndef DENSE_H
#define DENSE_H
#include<Dense>
#endif

#ifndef EIGEN_H
#define EIGEN_H
#include<Eigen>
#endif

#ifndef FSTREAM_H
#define FSTREAM_H
#include<fstream>
#endif

using namespace std;
using namespace Eigen;

void write(Matrix3d a, Matrix3d c){
	ofstream ofs("matriz", ios::binary);
	ofs.write((char *)&a, sizeof(a));
	ofs.write((char *)&c, sizeof(c));
	ofs.close();
}

void read(Matrix3d *b, Matrix3d *d){
	ifstream ifs("matriz", ios::binary);
	ifs.read((char*)b, sizeof(*b));
	ifs.read((char*)d, sizeof(*d));
	ifs.close();
}
int main()
{
	Matrix3d a=Matrix3d::Random(3,3);
	Matrix3d c=Matrix3d::Random(3,3);
	cout << a << endl << endl;
	cout << c << endl << endl;
	write(a,c);
	Matrix3d b,d;
	read(&b,&d);
	cout << b << endl << endl;
	cout << d << endl;
	return(0);
}
