
#pragma once
#include <vector>
#include <string>
#include <algorithm> //for max_element
#include <cassert>
#include "curve.h"
#include "instruments.h"
#include "linalgebra.h"


void solveNewtonRaphson(Curve& C, FixingInstruments &fixingInstruments, SwapInstruments& swapInstruments, std::string interpMethod)
{
	py::print("solving now");
	unsigned int numOfFixings = fixingInstruments.size();
	unsigned int numOfSwaps = swapInstruments.size();
	//for (auto& p : C.points_values)	py::print("points and values", p.first, " ... " , p.second);

	py::print("number of fixing instruments is: ", numOfFixings);
	py::print("number of swap instruments is: ", numOfSwaps);
	std::vector<double> Pri(numOfSwaps);
	std::vector<std::vector<double>> Jac(numOfSwaps);

	bool end_criteria;

	for (int i = 0; i < 5; ++i)
	{
		py::print("iteration number: ", i);
		end_criteria = true;

		for (unsigned int m = 0; m < numOfSwaps; ++m)
		{

			InstrumentPG<SwapInstrument> instrumentPG(swapInstruments[m]);
			std::tuple<double, std::vector<double>> aux;
			py::print("iteration number: ", i, "instrument number", m);
			//INTERPOLATION METHOD CHOICE

			aux = instrumentPG.pGLinZerosSingleCurveFutures(C.points_values, numOfFixings);
			/*
			py::print("aux first coordinate",-std::get<0>(aux));
			//py::print("aux second coordinate", std::get<1>(aux));

			Pri[m] = -std::get<0>(aux);
			Jac[m] =  std::get<1>(aux);
			end_criteria &= abs(std::get<0>(aux)) < 0.00001;
			py::print("price", Pri[m]);
			*/
		}
		/*
		py::print("C values before!");
		py::print(C.points_values);
		std::vector<double> soln = solveEqn(Jac, Pri);

		int j = 0;
		unsigned int jj = 0;
		for (auto k = C.points_values.begin(), e = C.points_values.end(); k != e; ++k) //these are pointers, so dont worry about not-referencing!
		//for (auto & k : C.points_values)
		{
			if (jj < numOfFixings) { ++jj; continue; }
			//py::print("inside solver, after algebra", k->second);
			k->second += soln[j];
			j += 1;

		}
		py::print("C values after!");
		py::print(C.points_values);

		if (end_criteria == true)
		{
			py::print("... accuracy achieved in", i, " iterations.");
			break;
		}
		*/
	}
	if (end_criteria != true) py::print("... maximum iterations!");


}
#pragma once
