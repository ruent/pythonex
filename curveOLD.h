#pragma once
#include <string> 
#include <vector>
#include <map>
#include "dates.h"
#include "instruments.h"
namespace py = pybind11;

class Curve
{
public:

	//std::vector<SwapInstrument> swapInstruments;
	std::map<int, double, std::less<int>> points_values; //curvePoints;curveValues;

	//InstrumentPGs priceGradients;

	explicit Curve(std::vector<SwapInstrument>& swapInstruments) : swapInstruments(swapInstruments) //curvePoints(curvePoints),
	{
		//for each instrument create price and derivative functors
		 
		for (auto irument : swapInstruments)
		{
			points_values.insert(std::pair<int, double>(irument.getmaturityPoint(), 0.5));
			py::print("before init: matpoint:", irument.getmaturityPoint(), "numofperiods", irument.getnumPeriods());
			InstrumentPG<SwapInstrument> priceGradient( irument );
			py::print("after init: matpoint:", priceGradient.instrument.getmaturityPoint(), "numofperiods", priceGradient.instrument.getnumPeriods());
			priceGradients.push_back(priceGradient);
			//std::map<int, double, std::less<int>> aux = { {728, 0.5}, {1092, 0.5} };
			//priceGradient(aux);
			
		}
		py::print("... curve initiated.");
		//pass above data to the solver, pass as constant reference
	}


	double dRate(int date) const
	{
		return 0;
	}

	std::map<int, double, std::less<int>> getPointValues()
	{
		return points_values;
	}

	std::vector<SwapInstrument> getCurveInstruments()
	{
		return swapInstruments;
	}
};
