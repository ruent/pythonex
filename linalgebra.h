 #pragma once

#include <vector>
#include <algorithm> //for max_element
#include <cassert>
 
std::vector<double> solveEqn(std::vector<std::vector<double>> &Jac, std::vector<double> &Pri)
{
	 
	unsigned int n = Jac.size();
	unsigned int m = Jac[0].size();
	//py::print("inside lin algebra, num of entries for Jac is: ", n, " and Jac[0] has num of entries: ", m);
	 
	//this code is for square systems so
	//assert(n == m);

	//construct the augmented matrix by extending Jac: Jac = [Jac:Pri] 
	for (unsigned int i = 0; i < n; ++i)
	{
		Jac[i].push_back(Pri[i]);
	}
	 

	//find maximum element of the absolute values in a row, for each row
	std::vector<double> max_in_a_row_all;
	double val_max = 0.0;// *std::max_element(Jac[0].begin(), Jac[0].end());
	double val_min = 0.0;// *std::max_element(Jac[0].begin(), Jac[0].end());

	for (unsigned int i = 0; i < n; ++i)
	{
		val_max = *std::max_element(Jac[i].begin(), Jac[i].end()); //max_element returns a pointer so you dereference it.
		val_min = *std::min_element(Jac[i].begin(), Jac[i].end());
		max_in_a_row_all.push_back(std::max<double>(abs(val_max), abs(val_min)));
		assert(max_in_a_row_all[i] > 0); //if not then there is no unique soln; you might write an error class here
		//if (max_in_a_row_all)
		//{raise() }
	}
	/*
	py::print("augmented jacobian");
	py::print(Jac[0][0], Jac[0][1], Jac[0][2], Jac[0][3]);
	py::print(Jac[1][0], Jac[1][1], Jac[1][2], Jac[1][3]);
	py::print(Jac[2][0], Jac[2][1], Jac[2][2], Jac[2][3]);
	py::print("maximums in each row", max_in_a_row_all);
	py::print("can you confirm?");
	*/

	//first iterate over columns: 
	for (unsigned int j = 0; j < m - 1; ++j)
	{
		//find the index of the maximum of abs(element) in a column, for each column
		double relative_col_max = 0;
		double relative_col_max_prev = 0;
		int index_of_max = j;
		for (unsigned int i = j; i < n; ++i)
		{
			if (abs(Jac[i][j]) / max_in_a_row_all[i] > relative_col_max)
			{
				relative_col_max = abs(Jac[i][j]);
				if (relative_col_max_prev < relative_col_max)
				{
					index_of_max = i;
				}
				relative_col_max_prev = relative_col_max;
			}
		}
	
		std::vector<double> temp(Jac[j]);
	
		Jac[j] = Jac[index_of_max];
		Jac[index_of_max] = temp;
		double xxx = 0;
		for (unsigned int i = j + 1; i < n; ++i)
		{
			xxx = Jac[i][j] / Jac[j][j];
			for (unsigned int k = 0; k < m + 1; ++k)
			{
				Jac[i][k] = Jac[i][k] - xxx * Jac[j][k];
			}
		}

	}
	/*
	py::print("augmented jacobian after row switch (turned off now) and Gaussian elimination");
	py::print(Jac[0][0], Jac[0][1], Jac[0][2], Jac[0][3]);
	py::print(Jac[1][0], Jac[1][1], Jac[1][2], Jac[1][3]);
	py::print(Jac[2][0], Jac[2][1], Jac[2][2], Jac[2][3]);
	*/ 
	std::vector<double> solution(n);
	//solution.resize(Pri.size());

	solution[n - 1] = Jac[n - 1][n] / Jac[n - 1][n - 1];
	for (unsigned int i = 2; i <= n; ++i)
	{
		double aux = 0;
		for (unsigned int j = n - i + 1; j < n; ++j)
		{
			aux += Jac[n - i][j] * solution[j];
		}
		solution[n - i] = (Jac[n - i][n] - aux) / Jac[n - i][n - i];
	}
	  
	return solution;
	//return Pri;


}
