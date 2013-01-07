/**
 * @file main.cpp
 * @author  Ataias Pereira Reis <ataiasreis@gmail.com>
 * @version 2.0.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * Este arquivo é o programa principal que faz chamadas de todas as outras
 * funções deste projeto que lida com equações diferenciais parciais.
 * O intuito de uso é no projeto de iniciação científica de alunos do Vortex
 * de iniciação científica na área de fluidos.
 * Ele faz uso da biblioteca Eigen e implementa alguns algoritmos para resolver
 * equações diferenciais parciais, incluindo a equação de Navier-Stokes em malhas
 * quadradas, a equaçaõ é a seguinte:
 * \f$\rho\left( \frac{\partial \textbf{v}}{\partial t}+
 \textbf{v}\cdot\nabla\textbf{v}\right)=
 -\nabla p+\mu\nabla^2\textbf{v}+\textbf{f}\f$
 */

#include<stdheader.hpp>


int main(){
	std::cout << "Testing Poisson Solver";
	return(0);
}
