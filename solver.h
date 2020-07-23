 
#pragma once
#include <vector>
#include <string>
#include <algorithm> //for max_element
#include <cassert>
#include "curve.h"
#include "instruments.h"
#include "linalgebra.h"
 
//template<typename TInterp>
class Solver
{
public:
	Solver() {}

	//if Cross Cur Swap: libor=Domestic Float, OISCurve is domestic OIS, C is the foreign projection to solve for
	//if Tenor basis swap: libor is known forward curve (say 3M), OIScurve is irrelavant (could be any curve), and C is the unknown forward curve (say 6M), foreignOIS is the common discount factor
	//the naming is for a Cross Currency swap so dont be fooled by it: for tenor basis swap it is just the common discount factor on both legs
	void NewtonRaphson(Curve& C, Curve& OIS, Curve& domLibor, Curve& domOIS, FixingInstruments &fixingInstruments, SwapInstruments& swapInstruments,  int numOfIter = 1000, double spotFxRate = 0.0, double tol = 0.000001)
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
		
		for (int i = 0; i < numOfIter; ++i)
		{
			//py::print("iteration number: ", i );
			end_criteria = true;
			GoodInterp phi = GoodInterp(C.points_values); //pass this to Cubic PG
			//py::print("an aux interpolation in solver: ", phi.Interp(44000));
			 
			for (unsigned int m = 0; m < numOfSwaps; ++m)
			{
				InstrumentPG<SwapInstrument> instrumentPG(swapInstruments[m]);
				bool finalReturn = false;
				if ((swapInstruments[m].paymType == "flofloCCS" || swapInstruments[m].paymType == "flofloCCSCubS") && swapInstruments[m].numOfPeriods != 1) finalReturn = true;

				//py::print("calling PG with swap number: ", m);
				auto aux = instrumentPG.defaultPG(C.points_values, numOfFixings, finalReturn, spotFxRate, C.curveSpotDate, OIS, domLibor, domOIS, phi);
				//auto aux = instrumentPG.pGCurveCCSCubS(C.points_values, numOfFixings, finalReturn, spotFxRate, C.curveSpotDate, OIS, domLibor, domOIS, phi);
				
				
				//py::print("before iterations", C.points_values);
				//py::print( numOfFixings, finalReturn, spotFxRate);
				//py::print(  C.curveSpotDate, domOIS.curveSpotDate, domLibor.curveSpotDate, foreignOIS.curveSpotDate);
				Pri[m] = -std::get<0>(aux);
				Jac[m] = std::get<1>(aux);
				//py::print("gradient in the solver for swap number: ", m);
				if (m== numOfSwaps-1 && i == numOfIter -1)	py::print(Jac[m]);
				end_criteria &= abs(std::get<0>(aux)) < tol;
				   
			}
			 
			 
			std::vector<double> soln = solveEqn(Jac, Pri); // J(x) y = -F for solving the root, i.e. F(x) = 0

			int j = 0;
			unsigned int jj = 0;
			for (auto k = C.points_values.begin(), e = C.points_values.end(); k != e; ++k) //these are pointers, so dont worry about not-referencing!//for (auto & k : C.points_values)
			{
				if (jj < numOfFixings) { ++jj; continue; }
				
				k->second += soln[j];
				//py::print("xx inside solver, after algebra", k->second, soln[j]);
				if (k->second < 0.1) 
				{
					py::print("opps, discount rate became negative:", k->second);
					k->second = 1.09; //discount can go negative! but then in log this will be trouble...
				}				
				if (k->second > 1.3) {
					k->second = 1.01; //discount can go too high
					py::print("opps, discount rate went too high!");
				}
				 
				j += 1;

			}
			//py::print("after iterations",C.points_values);
			if (end_criteria == true)
			{
				py::print("... accuracy achieved in", i, " iterations.");
				py::print("last price vector was: ", Pri);
				break;
			}
			 
		}
		if (end_criteria != true) {
			py::print("... maximum iterations!");
			py::print("last price vector was: ", Pri);
			//py::print("last grad vector modified in the LinAlgebra module is: ", Jac);
		}
    }

	void SteepestDescent(Curve& C, Curve& OIS, Curve& domLibor, Curve& domOIS, FixingInstruments &fixingInstruments, SwapInstruments& swapInstruments, int numOfIter = 1000, double spotFxRate = 0.0, double tol = 0.000001, bool pDetail = false)
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

		

		for (int i = 0; i < numOfIter; ++i)
		{
			//py::print("iteration number: ", i );
			end_criteria = true;
			double gA = 0; //g is sum_i price(instrumnet_i)^2
			double gB = 0; //g is sum_i price(instrumnet_i)^2
			double gC = 0; //g is sum_i price(instrumnet_i)^2
			double gZero = 0; //g is sum_i price(instrumnet_i)^2
			std::vector<double> gGrad(numOfSwaps);//gradient of g, so this is 2*Jac^T * Pri
			double gGradNorm;

			std::vector<double> curvePointsV(numOfSwaps);

			GoodInterp phi = GoodInterp(C.points_values); //pass this to Cubic PG

			for (unsigned int m = 0; m < numOfSwaps; ++m)
			{
				InstrumentPG<SwapInstrument> instrumentPG(swapInstruments[m]);
				bool finalReturn = false;
				if ((swapInstruments[m].paymType == "flofloCCS" || swapInstruments[m].paymType == "flofloCCSCubS") && swapInstruments[m].numOfPeriods != 1) finalReturn = true;

				
				auto aux = instrumentPG.defaultPG(C.points_values, numOfFixings, finalReturn, spotFxRate, C.curveSpotDate, OIS, domLibor, domOIS, phi);
				//py::print("before iterations", C.points_values);
				//py::print( numOfFixings, finalReturn, spotFxRate);
				//py::print(  C.curveSpotDate, domOIS.curveSpotDate, domLibor.curveSpotDate, foreignOIS.curveSpotDate);
				Pri[m] = std::get<0>(aux);
				Jac[m] = std::get<1>(aux);
				//py::print("gradient in the solver, swap: ", m);
				//py::print(Jac[m]);
				end_criteria &= abs(std::get<0>(aux)) < 0.000001;

				//get g and gGrad
				gA += Pri[m] * Pri[m];

			}
			if (pDetail == true) py::print("gA", gA);
			py::print("gA", gA);
			//get gradinet of g
			for (unsigned int m = 0; m < numOfSwaps; ++m)
			{
				for (unsigned int mm = 0; mm < numOfSwaps; ++mm)
				{
					gGrad[m] += 2*Pri[m] * Jac[mm][m];
				}

				gGradNorm += gGrad[m] * gGrad[m];
			}

			gGradNorm = sqrt(gGradNorm);
			if (gGradNorm == 0) 
			{
				end_criteria = true;
				break;
			}

			//normalize gradient vector
			for (unsigned int m = 0; m < numOfSwaps; ++m)
			{
				gGrad[m] = gGrad[m] / gGradNorm;
			}

			double val_max = *std::max_element(gGrad.begin(), gGrad.end());
			double val_min = *std::min_element(gGrad.begin(), gGrad.end());
			double bound = std::max<double>(abs(val_max), abs(val_min));
			double cpmin = 2.0;			
			for (auto k = C.points_values.begin(), e = C.points_values.end(); k != e; ++k)
			{
				if (k->second < cpmin) cpmin = k->second;
			}
			py::print("max in grad and min curve point and bound on alpha (iteration)",i, bound, cpmin, cpmin/bound);
			if (pDetail == true) py::print("gradient norm : ", gGradNorm);
			//py::print("gradient", gGrad);
			double alphaA = 0;
			double alphaB;
			double alphaC = std::min<double>(cpmin/bound/1.1, 0.001);

			//calculate gC
			//for this first get new curve pointValues
			int j = 0;
			unsigned int jj = 0;
			for (auto k = C.points_values.begin(), e = C.points_values.end(); k != e; ++k) //these are pointers, so dont worry about not-referencing!//for (auto & k : C.points_values)
			{
				if (jj < numOfFixings) { ++jj; continue; }
				//py::print("inside solver, after algebra", k->second);
				curvePointsV[j] = k -> second; //SAVE CURVE POINT VALUES HERE
				if (pDetail == true)py::print("1- k->second, alphaC, gGrad[j]");
				if (pDetail == true)py::print(k->second, alphaC, gGrad[j]);
				k->second += -alphaC * gGrad[j];
				if (pDetail == true)py::print("1- k->second", k->second);
				j += 1;
			}
			//then put this back in instrumentPG
			GoodInterp phiNew = GoodInterp(C.points_values); //pass this to Cubic PG
			for (unsigned int m = 0; m < numOfSwaps; ++m)
			{
				InstrumentPG<SwapInstrument> instrumentPG(swapInstruments[m]);
				bool finalReturn = false;
				if ((swapInstruments[m].paymType == "flofloCCS" || swapInstruments[m].paymType == "flofloCCSCubS") && swapInstruments[m].numOfPeriods != 1) finalReturn = true;
				auto aux = instrumentPG.defaultPG(C.points_values, numOfFixings, finalReturn, spotFxRate, C.curveSpotDate, OIS, domLibor, domOIS, phiNew);
				Pri[m] = std::get<0>(aux);
				if (pDetail == true)py::print("inside solver m, pri[m]", m, Pri[m]);
				gC += Pri[m] * Pri[m];
			}
			double gC_old = gC;
			if (pDetail == true)py::print("gC", gC);
			py::print("alphaC", alphaC);
			while (gC >= gA)
			{
				gC = 0;
				alphaC = alphaC / 2;
				//calculate gC
				//for this first get new curve pointValues
				j = 0;
				jj = 0;
				for (auto k = C.points_values.begin(), e = C.points_values.end(); k != e; ++k) //these are pointers, so dont worry about not-referencing!//for (auto & k : C.points_values)
				{
					if (jj < numOfFixings) { ++jj; continue; }
					if (pDetail == true) py::print("2 inside solver, after algebra", k->second);
					//while (curvePointsV[j] - alphaC * gGrad[j] < 0.4)
					//	alphaC = alphaC / 2;
					k->second = curvePointsV[j]-alphaC * gGrad[j];
					if (pDetail == true)py::print("2 inside solver, after algebra", k->second);
					j += 1;
				}
				//then put this back in instrumentPG
				GoodInterp phiNewNew = GoodInterp(C.points_values); //pass this to Cubic PG
				for (unsigned int m = 0; m < numOfSwaps; ++m)
				{
					InstrumentPG<SwapInstrument> instrumentPG(swapInstruments[m]);
					bool finalReturn = false;
					if ((swapInstruments[m].paymType == "flofloCCS" || swapInstruments[m].paymType == "flofloCCSCubS") && swapInstruments[m].numOfPeriods != 1) finalReturn = true;
					auto aux = instrumentPG.defaultPG(C.points_values, numOfFixings, finalReturn, spotFxRate, C.curveSpotDate, OIS, domLibor, domOIS, phiNewNew);
					Pri[m] = std::get<0>(aux);

					gC += Pri[m] * Pri[m];
					if (pDetail == true)py::print("inside the loop gC", gC);

				}
				py::print("INSDIE WHILE alphaC", alphaC);
				py::print("INSDIE WHILE GC, a decrease of amount:", gC -gA);
				if (alphaC < tol/100000000000000)
				{
					py::print("exit alphaC", alphaC);
					end_criteria = true;
					py::print("... no improvement in g3 over g1! (alphaC)");
					break;
				}
			}
			py::print("gC", gC);

			alphaB = alphaC / 2;
			j = 0;
			jj = 0;
			for (auto k = C.points_values.begin(), e = C.points_values.end(); k != e; ++k) //these are pointers, so dont worry about not-referencing!//for (auto & k : C.points_values)
			{
				if (jj < numOfFixings) { ++jj; continue; }
				if (pDetail == true)py::print("3- inside solver, after algebra", k->second);
				k->second = curvePointsV[j]-alphaB * gGrad[j];
				if (pDetail == true)py::print("3- inside solver, after algebra", k->second);
				j += 1;
			}
			//then put this back in instrumentPG
			GoodInterp phiNewNewNew = GoodInterp(C.points_values); //pass this to Cubic PG
			for (unsigned int m = 0; m < numOfSwaps; ++m)
			{
				InstrumentPG<SwapInstrument> instrumentPG(swapInstruments[m]);
				bool finalReturn = false;
				if ((swapInstruments[m].paymType == "flofloCCS" || swapInstruments[m].paymType == "flofloCCSCubS") && swapInstruments[m].numOfPeriods != 1) finalReturn = true;
				auto aux = instrumentPG.defaultPG(C.points_values, numOfFixings, finalReturn, spotFxRate, C.curveSpotDate, OIS, domLibor, domOIS, phiNewNewNew);
				Pri[m] = std::get<0>(aux);
				Jac[m] = std::get<1>(aux);

				gB += Pri[m] * Pri[m];
			}

			//see p. 659, Burden-Faires
			double hA = (gB - gA) / alphaB;
			double hB = (gC - gB) / (alphaC-alphaB);
			double hC = (hB - hA) / alphaC;
			double alphaZero = (alphaB - hA / hC) / 2.0;
			if (pDetail == true)py::print("gA,gB,gC", gA, gB, gC);
			if (pDetail == true)py::print("alphaA,alphaB,alphaC", alphaA, alphaB, alphaC);
			if (pDetail == true)py::print("hA,hB,hC, alphaZero", hA,hB,hC, alphaZero);
			j = 0;
			jj = 0;
			for (auto k = C.points_values.begin(), e = C.points_values.end(); k != e; ++k) //these are pointers, so dont worry about not-referencing!//for (auto & k : C.points_values)
			{
				if (jj < numOfFixings) { ++jj; continue; }
				if (pDetail == true)py::print("4 inside solver, after algebra k->second, curvePointsV[j]", k->second, curvePointsV[j]);
				k->second = curvePointsV[j] -alphaZero * gGrad[j];
				if (pDetail == true)py::print("4 inside solver, after algebra", k->second);
				j += 1;
			}
			//then put this back in instrumentPG
			GoodInterp phiNewNewNewNew = GoodInterp(C.points_values); //pass this to Cubic PG
			for (unsigned int m = 0; m < numOfSwaps; ++m)
			{
				InstrumentPG<SwapInstrument> instrumentPG(swapInstruments[m]);
				bool finalReturn = false;
				if ((swapInstruments[m].paymType == "flofloCCS" || swapInstruments[m].paymType == "flofloCCSCubS") && swapInstruments[m].numOfPeriods != 1) finalReturn = true;
				auto aux = instrumentPG.defaultPG(C.points_values, numOfFixings, finalReturn, spotFxRate, C.curveSpotDate, OIS, domLibor, domOIS, phiNewNewNewNew);
				Pri[m] = std::get<0>(aux);
				Jac[m] = std::get<1>(aux);

				gZero += Pri[m] * Pri[m];
				
			}
			if (pDetail == true)py::print("gZero, gC", gZero, gC);
			double alpha;
			double g;
			if (gZero < gC)
			{
				alpha = alphaZero;
				g = gZero;
			}
			else
			{
				alpha = alphaC;
				g = gC;
			}

			j = 0;
			jj = 0;
			for (auto k = C.points_values.begin(), e = C.points_values.end(); k != e; ++k) //these are pointers, so dont worry about not-referencing!//for (auto & k : C.points_values)
			{
				if (jj < numOfFixings) { ++jj; continue; }
				if (pDetail == true)py::print("5- inside solver, after algebra", k->second);
				k->second = curvePointsV[j] - alpha * gGrad[j];
				if (pDetail == true)py::print("5- inside solver, after algebra", k->second);
				j += 1;
			}

			 
			
			if (abs(g-gA)<tol)
			{
				end_criteria = true;
				py::print("... accuracy achieved in", i, " iterations!");
				py::print("last price vector was: ", Pri);
				break;
			}
			if (pDetail == true)py::print("i: gA vs g:", i, gA, g);
		}
		
		if (end_criteria != true) {
			if (pDetail == true)py::print("... maximum iterations!");
			if (pDetail == true)py::print("last price vector was: ", Pri);
			//py::print("last grad vector modified in the LinAlgebra module is: ", Jac);
		}
		if (pDetail == true)py::print("after after  the loop");
	}

};

