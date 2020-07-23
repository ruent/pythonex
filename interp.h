#pragma once
#include "linalgebra.h"

class GoodInterp
{
	//For log - linear interpolation between points(t_0, d_0) and (t_1, d_1), for a = (t_1 - t) / (t_1 - t_0)

public:
	const std::map<int, double, std::less<int>>& pointsValues;
	
	std::vector<std::vector<double>> X; // K * d'' = x WARNING X IS A COLUMN MATRIX! So TRASNPOSE BEFORE you do multiplication as X * d
	std::vector<double> h; //double since this will be used in division
	std::vector<int> timeG;
	unsigned int n;
	//GoodInterp() {}

	GoodInterp(const std::map<int, double, std::less<int>>& pointsValues) : pointsValues(pointsValues)
	{
		//fix the size
		n = pointsValues.size();

		//h.resize(n - 1);// size the h vector
		//timeG.resize(n);// size the timeGrid vector
		std::vector<double> auxPointsValues;// (n);
		unsigned int i = 0;
		int a;
		for (auto k = pointsValues.begin(), e = pointsValues.end(); k != e; ++k)
		{
			if (i > 0)
			{
				//h[i - 1] = k->first - a;
				h.push_back(k->first - a);
			}
			a = k->first;
			//auxPointsValues[i] = k->second;
			auxPointsValues.push_back(k->second);
			//timeG[i] = k->first;
			timeG.push_back(k->first);
			++i;
			
		}
		//find the X matrix
		//X.resize(n);
		 
		for (unsigned int i = 0; i < n; ++i)
		{
			std::vector<double> aux(n);
			for (unsigned int j=0;j<n;++j)
			{
				aux[j] = (j == i) ? 1 : 0;
			}
			
			//X[i] = finddDotDot(aux);
			X.push_back(finddDotDot(aux));
		}
		 
		//py::print("inside GoodInterp");
		//py::print("X[n - 1][n - 1]",X[n - 1][n - 1]);
		//verify that X works
		/*
		std::vector<double> soln = finddDotDot(auxPointsValues);
		py::print("direct soln: ", soln);
		std::vector<double> solnFromX(n);
		for (int i = 0; i < n; ++i)
		{
			double aux = 0;
			for (int j = 0; j < n; ++j)
			{
				aux += X[j][i] * auxPointsValues[j];
			}
			solnFromX[i] = aux;
		}
		
		py::print("soln from X: ", solnFromX);
		//test ends
		*/


	};

	std::vector<double> Interp(int date) const
	//void Interp(int date)
	{
		int lowerIndex = 0;
		for (unsigned int i = 0; i < n; ++i)
		{
			if (date > timeG[i]) lowerIndex = i; //NOT >=, PLEASE!!!!!!!!
			else break;
		}

		//std::vector<double> aux;// (n);
		//for (unsigned int j = 0; j < n; ++j)
		//{
		//	aux.push_back( 0);
		//}
		//return aux;
		//py::print("lowerIndex",lowerIndex);
		 
		std::vector<double> rowFromA(n);// a row of A
		for (unsigned int j = 0; j < n; ++j)
		{
			if (j == lowerIndex)
			{
				rowFromA[j] = (timeG[lowerIndex + 1] - date) / h[lowerIndex];
				//py::print("n, j, (timeG[lowerIndex + 1] - date) / h[lowerIndex]", n,j, rowFromA[j], h[lowerIndex]);
			}
			else {
				if (j == lowerIndex + 1)
				{
					rowFromA[j] = (date - timeG[lowerIndex]) / h[lowerIndex];
					//py::print("n, j,  (date - timeG[lowerIndex]) / h[lowerIndex], h[lowerIndex]", n, j, rowFromA[j], h[lowerIndex]);
				}
				else rowFromA[j] = 0;
			}
		}
		
		//py::print("h", h);
		//py::print("timegrid", timeG);
		//py::print("rowFromA", rowFromA);

		std::vector<double> rowFromB(n);// a row of A
		for (unsigned int j = 0; j < n; ++j)
		{
			double aux = (timeG[lowerIndex + 1] - date) / h[lowerIndex];
			if (j == lowerIndex) rowFromB[j] = h[lowerIndex] * h[lowerIndex] * (aux - aux * aux*aux) / 6.0;
			else rowFromB[j] = 0;
		}
		//py::print("rowFromB", rowFromB);
		//py::print("X",X);
		
		//calculate rowFromB *X
		std::vector<double> rowFromBX(n);// a row of A

		for (unsigned int j = 0; j < n; ++j)
		{
			double aux = 0;
			for (unsigned int i = 0; i < n; ++i)
			{
				aux += rowFromB[i] * X[j][i];
			}
			rowFromBX[j] = aux;
		}
		//py::print("rowFromBX", rowFromBX);
		
		std::vector<double> rowFromC(n);// a row of A
		for (unsigned int j = 0; j < n; ++j)
		{
			rowFromC[j] = 0;
		}
		double aux = (date - timeG[lowerIndex]) / h[lowerIndex];
		rowFromC[lowerIndex + 1] = h[lowerIndex] * h[lowerIndex] * (aux - aux * aux*aux) / 6.0;

		//py::print("rowFromC", rowFromC);
		 
		//calculate rowFromC *X
		std::vector<double> rowFromCX(n);// a row of A
		for (unsigned int j = 0; j < n; ++j)
		{
			double aux = 0;
			for (unsigned int i = 0; i < n; ++i)
			{
				aux += rowFromC[i] * X[j][i];
			}
			rowFromCX[j] = aux;
		}
		 
		std::vector<double> gradient(n);// 
		for (unsigned int j = 0; j < n; ++j)
		{
			gradient[j] = rowFromA[j] - rowFromBX[j] - rowFromCX[j];
		}

		//py::print("gradient", gradient);
		return gradient;
		 
	}

	std::vector<double> finddDotDot(const std::vector<double>&  auxPointsValues) const //k is the index of the perturbed curve point
	{
		//py::print("what the fuck!!!!!!!!!!!!!!!!!!");
		//py::print(pointsValues.size()); //WHY THE FUCK I Cannot even read this from here?? Because the reference has only information of the initial memory slot's address??
		 
		//fix the sizes
		std::vector<double> dDelta(n-1); //dDelta =[d_1 -d_0, d_2 - d_1 ...]
		 
		double b;
		for (unsigned int k=0; k<n;++k)
		{
			if (k > 0) 
			{
				dDelta[k - 1] = auxPointsValues[k] - b;			
			}
			b = auxPointsValues[k];
		}
				 
		std::vector<std::vector<double>> K(n); // K * d'' = x //row matrix, n-rows, n-columns

		std::vector<double> aux(n); 
		for (unsigned int j = 0; j < n; ++j)//construct first row of the matrix K
		{
			aux[j] = j==0 ? 1 : 0; //d_0'' = 0 
		}
		K[0]=aux;
		//py::print(K);
		
		//main part of K
		for (unsigned int i = 1; i < n-1; ++i) //for rows 1,2,..., n-2
		{
			for (unsigned int j = 0; j < n; ++j)//reset aux to zero
			{
				aux[j] = 0;  
			}
			aux[i - 1] = h[i - 1];//h_i-1 d_i-1''...
			aux[i] = 2 * (h[i - 1] + h[i]);
			aux[i + 1] = h[i];
			K[i] = aux;
		}
		
		 
		for (unsigned int j = 0; j < n; ++j)//construct last row of the matrix K
		{
			aux[j] = j == (n - 1) ? 1: 0; //d_0'' = 0 
		}
		K[n-1] = aux;
		//py::print(K);
		//construct x vector
		aux[0] = 0; //for d_0'' =0
		aux[n-1] = 0; //for d_n-1'' =0
		for (unsigned int i = 1; i < n - 1; ++i) //for rows 1,2,..., n-2
		{
			aux[i] = 6*(dDelta[i]/h[i]- dDelta[i-1] / h[i-1]);//h_i-1 d_i-1''...		
		}
		//now solve: K d'' = aux (x)
		//Note: solver modifies K so it is not reusable if you need them as they were
		std::vector<double> dDotDot = solveEqn(K, aux);
		//py::print(dDotDot);
		return dDotDot;
		  
	}
};

	
 

class LogLinInterp
{
	//For log - linear interpolation between points(t_0, d_0) and (t_1, d_1), for a = (t_1 - t) / (t_1 - t_0)

	int tA; //tA \leq tB
	int tB;
	double dA; 
	double dB;
	
	int t; //interpolation point
	int curveSpot;
public:
	double a; //calculated
	double d; //calculated
	double zeroR;

	LogLinInterp(int tA, int tB, double dA, double dB, int t, int curveSpot): tA(tA), tB(tB), dA(dA), dB(dB),  t(t), curveSpot(curveSpot)
	{
		a = (tB - t) / double(tB - tA);
		d = exp(a*log(dA) + (1 - a) * log(dB));
		zeroR = -log(dB / dA) / (tB - tA); //this should use year frac!! it is not impacting the results since i dont do extrapolation...
	}

	LogLinInterp(int tA, int tB, double dA, double dB, int t, int curveSpot, double a, double d ): tA(tA), tB(tB), dA(dA), dB(dB), t(t), curveSpot(curveSpot), a(a), d(d)
	{
		;
	}

	double sensA() const//sensitivity with respect to point A
	{
		return a * d / dA;
	}
	double sensB() const //sensitivity with respect to point A
	{
		return (1-a) * d / dB;
	}
};

class LinZInterp
{
	//For log - linear interpolation between points(t_0, d_0) and (t_1, d_1), for a = (t_1 - t) / (t_1 - t_0)

	int tA; //tA \leq tB
	int tB;
	double dA;
	double dB;

	int t; //interpolation point
	int curveSpot;
public:
	double a; //calculated
	double d; //calculated
	double zeroR =0.0;

	LinZInterp(int tA, int tB, double dA, double dB, int t, int curveSpot) : 
		tA(tA), tB(tB), dA(dA), dB(dB), t(t), curveSpot(curveSpot)
	{
		a = (tB - t) / double(tB - tA);
		//py::print("a", a, "tA",tA, "tB", tB, "t", t, "curveSpot", curveSpot);
		//py::print("dA", dA, "dB", dB);
		double rA;
		if (tA == curveSpot) 	rA = 0.0;
		else rA = -log(dA) / double(tA - curveSpot+1) / 360.0; 

		double rB = -log(dB) / double(tB - curveSpot+1)/360.0;
		double r = a * rA + (1 - a) * rB;
		d = exp(-r * (t - curveSpot)/360.0);
		//py::print("rA", rA, "rB", rB, "r", r, "(t - curveSpot)/360.0", (t - curveSpot) / 360.0,"d", d);
	}

	LinZInterp(int tA, int tB, double dA, double dB, int t, int curveSpot, double a, double d) : 
		tA(tA), tB(tB), dA(dA), dB(dB), t(t), curveSpot(curveSpot), a(a), d(d)
	{
		;
	}

	double sensA() const//sensitivity with respect to point A
	{
		//d = e^{-r (t-t0)}
		//d' = - d r' (t-t0);
		//r' = a rA' + (1-a) rB', etc...
		//py::print("d * (t - curveSpot) * a / dA / double(tA- curveSpot)", d, (t - curveSpot), a, dA, double(tA - curveSpot));
		return d * (t - curveSpot) * a / dA / double(tA- curveSpot+1) ;

	}
	double sensB() const //sensitivity with respect to point A
	{
		//py::print("d * (t - curveSpot) * a / dB / double(tB- curveSpot)", d, (t - curveSpot), a, dB, double(tB - curveSpot));
		return d * (t - curveSpot) * (1-a) / dB / double(tB - curveSpot+1);
	}
};