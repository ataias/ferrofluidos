/**
 * @file stdheader.hpp
 * @author Ataias Pereira Reis
 * Created on: Dec 29, 2012
 * Este arquivo contém arquivos cabeçalho básicos para inserir em qualquer
 * arquivo que faça parte do projeto. Isso evita o eclipse ficar reclamando
 * da ausência de funções.
 */

#ifndef STDHEADER_HPP_
#define STDHEADER_HPP_

	/*------------C and C++ libraries---------------*/
	#ifndef IOSTREAM_H
	#define IOSTREAM_H
	#include<iostream>
	#endif /* IOSTREAM_H */

	#ifndef STRING_H
	#define STRING_H
	#include<string>
	#endif/* STRING_H */

	#ifndef CSTDLIB_H
	#define CSTDLIB_H
	#include<cstdlib>
	#endif /* CSTDLIB_H */

	#ifndef CTIME_H
	#define CTIME_H
	#include<ctime>
	#endif /* CTIME_H */

	/*To deal with files*/
	#ifndef FSTREAM_H
	#define FSTREAM_H
	#include<fstream>
	#endif

	#ifndef VECTOR_H
	#define VECTOR_H
	#include<vector>
	#endif
	/*----------------------------------------------*/


	/*---------------BOOST LIBRARIES----------------*/
	#ifndef LEXICAL_CAST_H
	#define LEXICAL_CAST_H
	#include <boost/lexical_cast.hpp>
	#endif
	/*----------------------------------------------*/

	/*---------------Eigen Libraries --------------*/
	#ifndef DENSE_H
	#define DENSE_H
	#include<eigen3/Eigen/Dense>
	#endif

	#ifndef SPARSE_H
	#define SPARSE_H
	#include<eigen3/Eigen/Sparse>
	#endif

	/*---------------------------------------------*/

	/*----Python-----------------------------------*/
	#include<python2.7/Python.h>
	/*---------------------------------------------*/

	#include <../include/NavierStokes.hpp>
	#include <../include/Poisson.hpp>

	#define PENTADIAGONAL 5
	#define READ 0
	#define WRITE 1
	#define BINARY 0
	#define ASCII 1
	#define SUCCESS 1

  /** Os seguintes typedef's foram criados com o intuito de serem usados
  * na seção de matrizes esparsas.
  * */
	typedef Eigen::SparseMatrix<double> SpMat;
	typedef Eigen::Triplet<double> T;
	/*---------------------------------------------*/
#endif /* STDHEADER_HPP_ */
