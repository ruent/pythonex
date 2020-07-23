#pragma once
#pragma once
#include <string> 
#include <vector>
#include <map>
#include "dates.h"

#include "dates.h"
namespace py = pybind11;

//#include "instrumentpg.h"
//#include "instruments.h"
typedef std::vector<SwapInstrument> SwapInstruments;
typedef std::vector<FixingInstrument> FixingInstruments;
#include "interp.h" //only for cubic spline

class Curve
{
public:

	std::map<int, double, std::less<int>> points_values; //curvePoints;curveValues;
	int curveSpotDate;
	//GoodInterp goodInterp;

	explicit Curve(FixingInstruments &fixingInstruments, SwapInstruments &swapInstruments)
	{
		for (auto irument : fixingInstruments)
		{
			
			points_values.insert(std::pair<int, double>(irument.curvePoint, irument.curvePointV)); //0.5 is an initial value for the discount rate
			//what is the size after this?
		}
		for (auto irument : swapInstruments)
		{
			points_values.insert(std::pair<int, double>(irument.curvePoint, irument.curvePointV)); //0.5 is an initial value for the discount rate
			//what is the size after this?
		}
		//GoodInterp aux(points_values);
		//goodInterp = aux;
		curveSpotDate = points_values.begin()->first;
		py::print(" ... curve initialized.");
	}

	
	double dRateLogL(int date) const
	{
		double a = 0;
		double d = 0;
		int foundfirstday = 0;
		double foundfirstdisc = 0.0;
		double foundfirstzerorate = 0.0;


		int an_unfortunate_index = 0;
		for (auto k = points_values.begin(), e = points_values.end(); k != e; ++k)
		{
			an_unfortunate_index += 1;
			if (foundfirstday != 0)
			{
				if (date <= k->first)
				{
					a = (k->first - date) / double(k->first - foundfirstday);
					d = exp(a*log(foundfirstdisc) + (1 - a) * log(k->second));
					assert(count > 0);
					//py::print("an upper lower bound located for quarter:", i);
					break;
				}
				else
				{
					if (an_unfortunate_index == points_values.size())
					{
						//discount is assumed to be this one:
						d = foundfirstdisc * exp(-foundfirstzerorate * (date - foundfirstday));
						break;
					}
					else
					{
						foundfirstzerorate = -log(k->second / foundfirstdisc) / (k->first - foundfirstday);
						foundfirstday = k->first;
						foundfirstdisc = k->second;
						//py::print("lower bound should be updated:", i);
					}
				}
			}

			if (date <= k->first && foundfirstday == 0)
			{
				a = (k->first - date) / (double(k->first) - double(curveSpotDate));
				d = exp((1 - a) * log(k->second));//k->second = exp(k_first * f), f = -log(k_second)/k_first, d = 
				//py::print("smaller than the smallest:");
				break;
			}
			if (date > k->first && foundfirstday == 0)
			{
				foundfirstday = k->first;
				foundfirstdisc = k->second;
				foundfirstzerorate = 1; //in ++iter, this will be previous periods fixed inst. fwd rate; You can delete this line all together!
				//py::print("a candidate for lower bound:", i);
			}

		}
		return d;
	}

	double dRateCubS(int date) 
	{
		GoodInterp goodInterp(points_values);
		std::vector<double> gradient = goodInterp.Interp(date);
		double aux = 0;
		int i = 0;
		for (auto k = points_values.begin(), e = points_values.end(); k != e; ++k)
		{
			aux += k->second * gradient[i];
			++i;
		}
		return aux;
	}
	 

	double fRateLogL(int date, std::string tenor, bool immShift, std::string dayCountConv, bool isArithAver = false)
	{
		double x = dRateLogL(date);
		DayCountCalculator dayCCFloat = DayCountCalculator(dayCountConv);
		double yFrac;

		int fwdEndDate;
		if (isArithAver == false)
		{
			if (tenor == "3M")
			{	
				fwdEndDate = shiftByThreeMonths(date, immShift);
				yFrac = dayCCFloat.yFrac(date, fwdEndDate);
				return (x / dRateLogL(fwdEndDate) - 1)  / yFrac; //annualized fwd rate: multiply by 360/Act
			}
			else 
				if (tenor == "6M")
			    {
					fwdEndDate = shiftBySixMonths(date);
					yFrac = dayCCFloat.yFrac(date, fwdEndDate);
					return (x / dRateLogL(fwdEndDate) - 1) / yFrac; //annualized fwd rate: multiply by 360/Act
			    }
				else
				{
					if (tenor == "1M")
					{
						fwdEndDate = shiftByOneMonth(date);
						yFrac = dayCCFloat.yFrac(date, fwdEndDate);
						return (x / dRateLogL(fwdEndDate) - 1) / yFrac; //annualized fwd rate: multiply by 360/Act
					}
					else
					{
						py::print("Opps! I don't have that tenor!");
						return -1.0;
					}
				}
		}
		else
		{
			fwdEndDate = shiftByThreeMonths(date, immShift);
			yFrac = dayCCFloat.yFrac(date, fwdEndDate);
			return log(x / dRateLogL(fwdEndDate)) / yFrac; //annualized fwd rate: multiply by 360/Act
		}
			

	}

	double fRateCubS(int date, std::string tenor, bool immShift, std::string dayCountConv, bool isArithAver = false)
	{
		double x = dRateCubS(date);
		DayCountCalculator dayCCFloat = DayCountCalculator(dayCountConv);
		double yFrac;

		int fwdEndDate;
		if (isArithAver == false)
		{
			if (tenor == "3M")
			{
				fwdEndDate = shiftByThreeMonths(date, immShift);
				yFrac = dayCCFloat.yFrac(date, fwdEndDate);
				return (x / dRateCubS(fwdEndDate) - 1) / yFrac; //annualized fwd rate: multiply by 360/Act
			}
			else
				if (tenor == "6M")
				{
					fwdEndDate = shiftBySixMonths(date);
					yFrac = dayCCFloat.yFrac(date, fwdEndDate);
					return (x / dRateCubS(fwdEndDate) - 1) / yFrac; //annualized fwd rate: multiply by 360/Act
				}
				else
				{
					if (tenor == "1M")
					{
						fwdEndDate = shiftByOneMonth(date);
						yFrac = dayCCFloat.yFrac(date, fwdEndDate);
						return (x / dRateCubS(fwdEndDate) - 1) / yFrac; //annualized fwd rate: multiply by 360/Act
					}
					else
					{
						py::print("Opps! I don't have that tenor!");
						return -1.0;
					}
				}
		}
		else
		{
			fwdEndDate = shiftByThreeMonths(date, immShift);
			yFrac = dayCCFloat.yFrac(date, fwdEndDate);
			return log(x / dRateCubS(fwdEndDate)) / yFrac; //annualized fwd rate: multiply by 360/Act
		}


	}

	 
};





