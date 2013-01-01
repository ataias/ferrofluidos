/*
 * edp_06.cpp
 *
 *  Created on: Dec 29, 2012
 *      Author: ataias
 *
 *      Testes sobre salvar matrizes da Eigen em arquivos binários
 */

#include<stdheader.hpp>

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
	cout << sizeof(a)/9 << " " << sizeof(double) << endl;
	write(a,c);
	Matrix3d b,d;
	read(&b,&d);
	cout << b << endl << endl;
	cout << d << endl;
	return(0);
}
