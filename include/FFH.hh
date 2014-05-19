/**
 * @file stdheader.hpp
 * @author Ataias Pereira Reis
 * Created on: Dec 29, 2012
 * Este arquivo contém arquivos cabeçalho básicos para inserir em qualquer
 * arquivo que faça parte do projeto. Isso evita o eclipse ficar reclamando
 * da ausência de funções.
 */

#ifndef FFH_HH_
#define FFH_HH_

	/*------------C and C++ libraries---------------*/
	#ifndef IOSTREAM_H
	#define IOSTREAM_H
	#include<iostream>
	#endif /* IOSTREAM_H */
	/*----------------------------------------------*/

	/*---------------Eigen Libraries --------------*/
	#ifndef DENSE_H
	#define DENSE_H
	#include<Eigen/Dense>
	#endif

	#ifndef SPARSE_H
	#define SPARSE_H
	#include<Eigen/Sparse>
	#endif
	/*---------------------------------------------*/

  /** Os seguintes typedef's foram criados com o intuito de serem usados
  * na seção de matrizes esparsas.
  * */
	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	/*---------------------------------------------*/

	using namespace Eigen;
	using namespace std;
#endif /* FFH_HH_ */
