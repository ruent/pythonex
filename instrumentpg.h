#pragma once
#include <cassert>
namespace py = pybind11;
#include "dates.h"
#include "instruments.h"
#include "interp.h"
#include "instrsens.h"
#include "curve.h"
#include "pricer.h"

//forward declaration
//class Curve;

/*
DEBUG LOG:
spreadFloatSensitivty at -1 was put wihoiur curent-cutoff>0
*/

//template <typename T, typename TInterp>
template <typename T >
class InstrumentPG
{
public:
	T &instrument;
	InstrumentPG(T &instrument) : instrument(instrument) 
	{ //py::print("yeah, instrumentPG created!"); 
	}; 

	std::tuple<double, std::vector<double>>  pGCurveFutures(std::map<int, double, std::less<int>>& points_values, int numOfFixings, int curveSpot, Curve &oisCurve) const
	{
		double pvFixedLeg = 0.0;
		double pvFloatLeg = 0.0;

		//initialize day count calculators
		DayCountCalculator dayCCFixed = DayCountCalculator(instrument.dayCFixed);
		DayCountCalculator dayCCFloat = DayCountCalculator(instrument.dayCFloat);
		double yFracFixed = 0;
		double yFracFloat = 0;

		//find pay date
		int paydate = instrument.spotDate;
		int auxpaydate = paydate;
		//py::print("shifting for IMM dates!");
		paydate = shiftByThreeMonths(paydate, true);
		yFracFixed = dayCCFixed.yFrac(auxpaydate, paydate);
		yFracFloat = dayCCFloat.yFrac(auxpaydate, paydate);

		//initialize counters etc.
		double a = 0; //for future reference, not used now!
		double d = 0;

		int foundfirst_cpdate = 0;
		double foundfirst_cpvalue = 0.0;
		double foundFirstZeroRate = 0.0;
		int count = -1; //counts the curve points smaller than the paydate	
		int an_unfortunate_index = 0;

		std::vector<double> Gradient(points_values.size() - numOfFixings); // this should have as many entries as points_values, except FIXINGS

		InstrSens floatSens = InstrSens(instrument);

		//iterate over curve points
		for (auto k = points_values.begin(), e = points_values.end(); k != e; ++k)
		{
			//ignore the fixing dates
			//increase count for proper gradient allocation
			double gradientPrev = 0;
			double gradientNow = 0;

			an_unfortunate_index += 1;
			if (an_unfortunate_index < numOfFixings)
			{
				//count += 1;
				foundfirst_cpdate = k->first;
				foundfirst_cpvalue = k->second;
				continue;
			}

			if (foundfirst_cpdate != 0)
			{
				if (paydate <= k->first)
				{
					//float pay = N(d_prevcurve - d)
					//fixed pay = - N fixedR yfracFixed d

					//always zero: a = (k->first - paydate) / double(k->first - foundfirst_cpdate);
					//d = exp(a*log(foundfirst_cpvalue) + (1 - a) * log(k->second));
					d = k->second;
					if (count > 0) Gradient[count - 1] += floatSens.floatSensA(foundfirst_cpvalue, k->second, oisCurve.dRateLogL(paydate)); // instrument.notional; //float paym sens due to d_prev
					//py::print("count", count);
					//py::print("gradientNow", gradientNow);
					Gradient[count] += floatSens.floatSensB(foundfirst_cpvalue, k->second, oisCurve.dRateLogL(paydate));// -instrument.notional; //float paym sens due to current point discount
					Gradient[count] += -instrument.fixedRate *yFracFixed*floatSens.fixedSensB();// -instrument.fixedRate *yFracFixed * instrument.notional; //fix paym sens due to current point discount
					//py::print("an upper lower bound located for quarter:", i);
					break;

				}
				else
				{
					if (an_unfortunate_index == points_values.size())
					{
						//futures will never be the end of the curve, but if it is, you will be here.
						//discount is assumed to be this one:
						d = foundfirst_cpvalue * exp(-foundFirstZeroRate * (paydate - foundfirst_cpdate));
						Gradient[count - 1] += instrument.notional *  (1 - exp(-foundFirstZeroRate * (paydate - foundfirst_cpdate))); //due to float payment
						Gradient[count - 1] += -instrument.fixedRate*yFracFixed * instrument.notional *exp(-foundFirstZeroRate * (paydate - foundfirst_cpdate)); //beyond last curve point
						break;
					}
					else
					{
						foundFirstZeroRate = -log(k->second / foundfirst_cpvalue) / (k->first - foundfirst_cpdate);
						foundfirst_cpdate = k->first;
						foundfirst_cpvalue = k->second;
						count += 1;
						//py::print("lower bound should be updated:", i);
					}
				}
			}

			if (paydate <= k->first && foundfirst_cpdate == 0)
			{
				py::print("you should never come here");
				//d is just the curve point (the first curve point)
				//float pay = N(1-d)
				//fixed pay = - N fixedR yfracFixed d
				break;
			}

			if (paydate > k->first && foundfirst_cpdate == 0)
			{
				//this is point of entry 
				foundFirstZeroRate = -log(k->second / foundfirst_cpvalue) / (k->first - foundfirst_cpdate);
				foundfirst_cpdate = k->first;
				foundfirst_cpvalue = k->second;
				//py::print("Gradient size", Gradient.size());
				//py::print(Gradient[0]);
				count += 1;
			}
			//Gradient[count] = gradientNow; //this is weird! why do I have to do this?

		}
		pvFixedLeg += instrument.fixedRate * yFracFixed * floatSens.fixedPV(d, oisCurve.dRateLogL(paydate));// instrument.notional*d;
		//py::print("foundfirst_cpvalue", foundfirst_cpvalue);
		//py::print("d", d);
		pvFloatLeg += (foundfirst_cpvalue / d - 1) * instrument.notional *oisCurve.dRateLogL(paydate); //no time frac here, be careful in future implementations
		//py::print("c++ pvFixedLeg: ", pvFixedLeg, "pvFloatLeg: ", pvFloatLeg, "yearF: ", yFracFixed);
		//py::print("price: ", -pvFixedLeg + pvFloatLeg, "Gradient!", Gradient);
		std::tuple<double, std::vector<double>> result(-pvFixedLeg + pvFloatLeg, Gradient);
		return result;
	}//

	//this one takes a discount curve and calculates another, as in USD OIS and USD LIBOR 
	std::tuple<double, std::vector<double>>  pGCurveSwap(std::map<int, double, std::less<int>>& points_values, int numOfFixings, int curveSpot, Curve &oisCurve) const
	{
		
		
		//initialize day count calculators
		DayCountCalculator dayCCFixed = DayCountCalculator(instrument.dayCFixed);
		DayCountCalculator dayCCFloat = DayCountCalculator(instrument.dayCFloat);
		double yFracFixed = 0;
		double yFracFloat = 0;
		double pvFixedLeg = 0.0;
		double pvFloatLeg = 0.0;
		
		

		assert(instrument.fixedTenor >= instrument.floatTenor);
		int alpha = instrument.fixedTenor / instrument.floatTenor; // 1 if 3M-3M or 2 if 6M-3M etc.
		int paydate = instrument.spotDate;
		
		double d_prev = 1; //initial discount //I use this in fwd rate calc; I don't know if dprev/d sacrifices any precision in 3month rate
		int paydate_prev = instrument.spotDate;
		int count_prev = 0; //this will be the least upper bound for the previous payment date
		double sensA_prev = 0;
		double sensB_prev = 0;

		std::vector<double> Gradient(points_values.size() - numOfFixings);
		
		for (int i = 1, e = instrument.numOfPeriods; i <= e; ++i)//i is used for parity; parity is for 3M-6M setup
		{
			//shift paydate while keeping a prev for the next period
			int auxpaydate = paydate;
			if (instrument.floatTenor == 3 || instrument.fixedTenor == 3)
				paydate = shiftByThreeMonths(paydate, false);
			else {
				if (instrument.floatTenor == 6 || instrument.fixedTenor == 6)
				{
					paydate = shiftBySixMonths(paydate);
				}
				else
				{
					if (instrument.floatTenor == 1 || instrument.fixedTenor == 1)
					{
						paydate = shiftByOneMonth(paydate);
					}
					else py::print("this tenor type is not coded for PriceJacobian!");
				}
			}

			yFracFixed = dayCCFixed.yFrac(paydate_prev, paydate);
			yFracFloat = dayCCFloat.yFrac(auxpaydate, paydate);
			paydate_prev = auxpaydate;

			InstrSens floatSens = InstrSens(instrument);

			double a = 0;
			double d = 0;
			int foundfirst_cpdate = 0;
			double foundfirst_cpvalue = 0.0;
			double foundFirstZeroRate = 0.0;
			int count = -1; //should be the index of the least upper bound point: 
							//e.g. if: last fixingdate < paydate < first non-fixing date, count = 0
			int gradientCut = numOfFixings - 1;
			int an_unfortunate_index = 0;
			for (auto k = points_values.begin(), e = points_values.end(); k != e; ++k)
			{
				an_unfortunate_index += 1;


				if (foundfirst_cpdate != 0)
				{
					LogLinInterp phi = LogLinInterp(foundfirst_cpdate, k->first, foundfirst_cpvalue, k->second, paydate, curveSpot);

					if (paydate <= k->first)
					{
						a = phi.a;
						d = phi.d;
						//Executive summary: below to lines for risk in N(dprev-d) due to d_prev; d_prev changes with curve[count_old] and curve[count_old-1]
						//different from linear interp you also need to bring two curve points with yourself here
						if (count >= gradientCut)
						{
							if (count_prev > gradientCut) {
								Gradient[count_prev - gradientCut - 1] += floatSens.floatSensA(foundfirst_cpvalue, k->second, oisCurve.dRateLogL(paydate)) *sensA_prev;
							}
							if (count_prev >= gradientCut) {
								Gradient[count_prev - gradientCut] += floatSens.floatSensA(foundfirst_cpvalue, k->second, oisCurve.dRateLogL(paydate)) *sensB_prev;
							}
							//Executive summary: below two lines for risk in floatPV, N(dprev-d), due to d; d changes with curve[count] and curve[count-1]
							if (count > gradientCut) {
								Gradient[count - gradientCut - 1] += floatSens.floatSensB(foundfirst_cpvalue, k->second, oisCurve.dRateLogL(paydate)) *phi.sensA(); //float paym sens due to current point discount (note: N (d_prev -d) = N( d_prev - a d_k-1 - (1-a)d_k))
							}
							Gradient[count - gradientCut] += floatSens.floatSensB(foundfirst_cpvalue, k->second, oisCurve.dRateLogL(paydate)) * phi.sensB(); //float paym sens due to current point discount
							if (i % alpha == 0 && count > gradientCut) {
								Gradient[count - gradientCut - 1] += -instrument.fixedRate *yFracFixed * floatSens.fixedSensA()* phi.sensA(); //fix paym sens due to previous point discount
							}
							if (i % alpha == 0) Gradient[count - gradientCut] += -instrument.fixedRate *yFracFixed * floatSens.fixedSensB() * phi.sensB(); //fix paym sens due to current point discount
						}
						sensA_prev = phi.sensA();
						sensB_prev = phi.sensB();
						break;
					}
					else
					{
						if (an_unfortunate_index == points_values.size())
						{
							py::print("extrapolation region!");
							py::print("to be coded yet!");
							break;
						}
						else
						{
							foundFirstZeroRate = phi.zeroR;
							foundfirst_cpdate = k->first;
							foundfirst_cpvalue = k->second;
							count += 1;
						}
					}
				}

				if (paydate <= k->first && foundfirst_cpdate == 0)
				{
					//this part is not at work at the moment, a curve will always start from today
					break;
				}
				if (paydate > k->first && foundfirst_cpdate == 0)
				{
					foundFirstZeroRate = 1000;
					foundfirst_cpdate = k->first;
					foundfirst_cpvalue = k->second;
					count += 1;
				}

			}
			 
			//done with loop over curve points, record, calculate
			if (i % alpha == 0) pvFixedLeg += instrument.fixedRate * yFracFixed * floatSens.fixedPV(d, oisCurve.dRateLogL(paydate));
			pvFloatLeg += (floatSens.floatPV(d_prev, d, oisCurve.dRateLogL(paydate)) + instrument.spreadFloat * yFracFloat * floatSens.floatSpreadPV(d, oisCurve.dRateLogL(paydate))); //no time frac here, be careful in future implementations
			d_prev = d; //d_prev is assigned right after an instrument is located on the curve points grid.
						//hence dprev is the latest discount, hence the latest discount for the latest paydate for the current instrument
			count_prev = count;			
		}
	
		std::tuple<double, std::vector<double>> result(-pvFixedLeg + pvFloatLeg, Gradient);
		return result;
	}//

	//find new fwd given disc fixed
	std::tuple<double, std::vector<double>>  pGCurveTBS(std::map<int, double, std::less<int>>& points_values, int numOfFixings, 
		bool finalReturn, double spotFxRate, int curveSpot, Curve &oisCurve, Curve &domLibor, Curve &domOIS) const
	{
	
		//initialize day count calculators
		DayCountCalculator dayCCFixed = DayCountCalculator(instrument.dayCFixed);
		DayCountCalculator dayCCFloat = DayCountCalculator(instrument.dayCFloat);
		double yFracFixed = 0;
		double yFracFloat = 0;
		double pvFixedLeg = 0.0;
		double pvFloatLeg = 0.0;
		int paydate = instrument.spotDate;

		double d_prev = 1; //initial discount //I use this in fwd rate calc; I don't know if dprev/d sacrifices any precision in 3month rate
		int count_prev = 0; //this will be the least upper bound for the previous payment date
		int paydate_prev = instrument.spotDate;
		double sensA_prev = 0;
		double sensB_prev = 0;
		std::vector<double> Gradient(points_values.size() - numOfFixings);

		//this code is only for tenor basis swaps or cross c
		assert(instrument.fixedRate == 0.0);
		//this code assumes fixed (or the unsolved for float) tenor is at least as large as the float tenor
		assert(instrument.fixedTenor >= instrument.floatTenor);
		
		int alpha = instrument.fixedTenor / instrument.floatTenor; // 1 if 3M-3M or 2 if 6M-3M etc.

		//i is used for parity; parity is for 3M-6M setup
		for (int i = 1, e = instrument.numOfPeriods; i <= e; ++i)
		{
			//shift paydate while keeping a prev for the next period
			int auxpaydate = paydate;
			if (instrument.floatTenor == 3 || instrument.fixedTenor == 3)
				paydate = shiftByThreeMonths(paydate, false);
			else {
				if (instrument.floatTenor == 6 || instrument.fixedTenor == 6)
					paydate = shiftBySixMonths(paydate);
				else py::print("this tenor type is not coded for PriceJacobian!");
			}


			yFracFixed = dayCCFixed.yFrac(paydate_prev, paydate);
			yFracFloat = dayCCFloat.yFrac(auxpaydate, paydate);
			paydate_prev = auxpaydate;

			InstrSens floatSens = InstrSens(instrument);

			double a = 0;
			double d = 0;
			int foundfirst_cpdate = 0;
			double foundfirst_cpvalue = 0.0;
			double foundFirstZeroRate = 0.0;
			int count = -1; //should be the index of the least upper bound point: 
							//e.g. if: last fixingdate < paydate < first non-fixing date, count = 0
			int gradientCut = numOfFixings - 1;
			int an_unfortunate_index = 0;
			for (auto k = points_values.begin(), e = points_values.end(); k != e; ++k)
			{
				an_unfortunate_index += 1;


				if (foundfirst_cpdate != 0)
				{
					LogLinInterp phi = LogLinInterp(foundfirst_cpdate, k->first, foundfirst_cpvalue, k->second, paydate, curveSpot);
					

					if (paydate <= k->first)
					{
						//py::print("A",foundfirst_cpdate, k->first, foundfirst_cpvalue, k->second, paydate, curveSpot);
						//py::print("B",curvepoilowd_prev, curvepoiupd_prev, curvepoilow_prev, curvepoiup_prev, paydate_prev, curveSpot, a_prev, d_prev);
						a = phi.a;
						d = phi.d;
						//Executive summary: below to lines for risk in N(dprev-d) due to d_prev; d_prev changes with curve[count_old] and curve[count_old-1]
						//different from linear interp you also need to bring two curve points with yourself here
						if (count >= gradientCut)
						{
							if (count_prev > gradientCut) {
								Gradient[count_prev - gradientCut - 1] += floatSens.floatSensA(foundfirst_cpvalue, k->second, oisCurve.dRateLogL(paydate)) *sensA_prev;
							}
							if (count_prev >= gradientCut) {
								Gradient[count_prev - gradientCut] += floatSens.floatSensA(foundfirst_cpvalue, k->second, oisCurve.dRateLogL(paydate)) *sensB_prev;
							}
							//Executive summary: below two lines for risk in floatPV, N(dprev-d), due to d; d changes with curve[count] and curve[count-1]
							if (count > gradientCut) {
								Gradient[count - gradientCut - 1] += floatSens.floatSensB(foundfirst_cpvalue, k->second, oisCurve.dRateLogL(paydate)) *phi.sensA(); //float paym sens due to current point discount (note: N (d_prev -d) = N( d_prev - a d_k-1 - (1-a)d_k))
								//py::print("3", -instrument.notional *phi.sensA());
							}
							Gradient[count - gradientCut] += floatSens.floatSensB(foundfirst_cpvalue, k->second, oisCurve.dRateLogL(paydate)) * phi.sensB(); //float paym sens due to current point discount

							if (i % alpha == 0 && count > gradientCut) {
								Gradient[count - gradientCut - 1] += -instrument.fixedRate *yFracFixed * floatSens.fixedSensA()* phi.sensA(); //fix paym sens due to previous point discount
								//py::print("4", -instrument.fixedRate *yFracFixed * instrument.notional* phi.sensA());
							}
							if (i % alpha == 0) Gradient[count - gradientCut] += -instrument.fixedRate *yFracFixed * floatSens.fixedSensB() * phi.sensB(); //fix paym sens due to current point discount
						}
						sensA_prev = phi.sensA();
						sensB_prev = phi.sensB();
						//py::print("an upper lower bound located for quarter:", i);
						break;
					}
					else
					{
						if (an_unfortunate_index == points_values.size())
						{
							//discount is assumed to be this one:
							py::print("extrapolation region!");
							py::print("to be coded yet!");
							break;
						}
						else
						{
							foundFirstZeroRate = phi.zeroR;
							foundfirst_cpdate = k->first;
							foundfirst_cpvalue = k->second;
							count += 1;
							//py::print("lower bound should be updated:", i);
						}
					}
					 
				}

				if (paydate <= k->first && foundfirst_cpdate == 0)
				{
					//this part is not at work at the moment, a curve will always start from today
					break;
				}
				if (paydate > k->first && foundfirst_cpdate == 0)
				{
					 
					foundFirstZeroRate = 1000;
					foundfirst_cpdate = k->first;
					foundfirst_cpvalue = k->second;
					count += 1;
					//py::print("found a lower point, and here d_prev is ", d_prev);
				}
			}
			//py::print("done with a pay date (went over all instruments), here I have d_prev, d as ", d_prev, d);
			//py::print("if d came as none check it out! from the prints above...");
			//you are done with loop over curve points, record, calculate and do to next swap paydate
			//fixed leg calculations are not necessary for float-Float
			if (i % alpha == 0) pvFixedLeg += instrument.fixedRate * yFracFixed * floatSens.fixedPV(d, oisCurve.dRateLogL(paydate));
			pvFloatLeg += (floatSens.floatPV(d_prev, d, oisCurve.dRateLogL(paydate)) + instrument.spreadFloat * yFracFloat * floatSens.floatSpreadPV(d, oisCurve.dRateLogL(paydate))); //no time frac here, be careful in future implementations
			d_prev = d; //d_prev is assigned right after an instrument is located on the curve points grid.
						//hence dprev is the latest discount, hence the latest discount for the latest paydate for the current instrument
			count_prev = count;
			
		}

		double pvFloatLegFinExch = finalReturn * instrument.notional * spotFxRate;

		SwapInstrument anInst(instrument.floatTenor, instrument.fixedTenor, "plain", "flofix",
			0.0, instrument.numOfPeriods / alpha, instrument.spotDate, instrument.dayCFixed, instrument.dayCFloat,
			instrument.spreadFloat, instrument.edControl, 0.9);

		double pvdomFloatLeg = instrumentPrice(anInst, domLibor, domOIS, finalReturn);
		//py::print("instrument.numOfPeriods", instrument.numOfPeriods);
		//py::print(instrument.curvePoint, "Eur 6M Leg PV ", pvdomFloatLeg);
		//py::print(instrument.curvePoint, "-pvFixedLeg should be zero ", -pvFixedLeg);
		//py::print(instrument.curvePoint, "pvFloatLeg 3M Euribor ", pvFloatLeg);
		//py::print( "pvFloatLegFinExch should be zero ", pvFloatLegFinExch);
		std::tuple<double, std::vector<double>> result(-pvFixedLeg + pvFloatLeg  -pvdomFloatLeg - pvFloatLegFinExch, Gradient);
		return result;
	}//

	//find new discount given fwdcurve fixed
	std::tuple<double, std::vector<double>>  pGCurveCCS(std::map<int, double, std::less<int>>& points_values, int numOfFixings,
		bool finalReturn, double spotFxRate, int curveSpot, Curve &euribor, Curve &domLibor, Curve &domOIS) const
	{
		//initialize day count calculators
		DayCountCalculator dayCCFixed = DayCountCalculator(instrument.dayCFixed);
		DayCountCalculator dayCCFloat = DayCountCalculator(instrument.dayCFloat);
		double yFracFixed = 0;
		double yFracFloat = 0; double pvFixedLeg = 0.0;
		double pvFloatLeg = 0.0;
		double pvFloatSpreadLeg = 0.0;
		
		int paydate_prev = instrument.spotDate;
		int paydate = instrument.spotDate;

		std::vector<double> Gradient(points_values.size() - numOfFixings);

		//this code is only for tenor basis swaps or cross c
		assert(instrument.fixedRate == 0.0);
		//this code assumes fixed (or the unsolved for float) tenor is at least as large as the float tenor
		assert(instrument.fixedTenor >= instrument.floatTenor);
		int alpha = 1;// instrument.fixedTenor / instrument.floatTenor; // 1 if 3M-3M or 2 if 6M-3M etc.
		  
		 
		double dfinal = 0;
		int indexForFinalExchnage;
		
		for (int i = 1, e = instrument.numOfPeriods; i <= e; ++i)//i is used for parity; parity is for 3M-6M setup
		{
			//shift paydate while keeping a prev for the next period
			int auxpaydate = paydate;
			if (instrument.floatTenor == 3 || instrument.fixedTenor == 3)
			{
				if (instrument.numOfPeriods == 1) 
					paydate = shiftByThreeMonths(paydate, true);
				else 
					paydate = shiftByThreeMonths(paydate, false);
			}
				
			else {
				if (instrument.floatTenor == 6 || instrument.fixedTenor == 6)
					paydate = shiftBySixMonths(paydate);
				else py::print("this tenor type is not coded for PriceJacobian!");
			}

			yFracFixed = dayCCFixed.yFrac(paydate_prev, paydate);
			yFracFloat = dayCCFloat.yFrac(auxpaydate, paydate);
			paydate_prev = auxpaydate;

			InstrSens floatSens = InstrSens(instrument);

			 
			double d = 0;
			int foundfirst_cpdate = 0;
			double foundfirst_cpvalue = 0.0;
			double foundFirstZeroRate = 0.0;
			int count = -1; //should be the index of the least upper bound point: 
							//e.g. if: last fixingdate < paydate < first non-fixing date, count = 0
			int gradientCut = numOfFixings - 1;
			int an_unfortunate_index = 0;
			
			 
			for (auto k = points_values.begin(), e = points_values.end(); k != e; ++k)
			{
				an_unfortunate_index += 1;	 
				if (foundfirst_cpdate != 0)
				{
					LogLinInterp phi = LogLinInterp(foundfirst_cpdate, k->first, foundfirst_cpvalue, k->second, paydate, curveSpot);
					if (paydate <= k->first)
					{
						//py::print("   ");
						//py::print("NEW PAYMENT, see the details:");
						//py::print("A",foundfirst_cpdate, k->first, foundfirst_cpvalue, k->second, paydate, curveSpot);
						//py::print("B",curvepoilowd_prev, curvepoiupd_prev, curvepoilow_prev, curvepoiup_prev, paydate_prev, curveSpot, a_prev, d_prev);
						d = phi.d;
						if (i == instrument.numOfPeriods)
						{
							dfinal = d;
							indexForFinalExchnage = an_unfortunate_index - numOfFixings;
						}
						//Executive summary: below to lines for risk in N(dprev-d) due to d_prev; d_prev changes with curve[count_old] and curve[count_old-1]
						//different from linear interp you also need to bring two curve points with yourself here
						 
						//py::print("euribor.dRateLogL(paydate_prev), euribor.dRateLogL(paydate), d");
						//py::print(euribor.dRateLogL(paydate_prev), euribor.dRateLogL(paydate), d);
						//py::print("a", a);
						//py::print("d", d);
						if (count >= gradientCut)
						{							
							if (count > gradientCut) 
							{
								Gradient[count - gradientCut - 1] += floatSens.floatSensB(euribor.dRateLogL(paydate_prev), euribor.dRateLogL(paydate), d) *phi.sensA(); //float paym sens due to current point discount (note: N (d_prev -d) = N( d_prev - a d_k-1 - (1-a)d_k))
								Gradient[count - gradientCut - 1] += instrument.notional*instrument.spreadFloat *yFracFloat *phi.sensA(); //spread sens
								Gradient[count - gradientCut - 1] += -instrument.notional*instrument.fixedRate *yFracFloat *phi.sensA(); //fixed sens, for futures only
							}
							
							Gradient[count - gradientCut] += floatSens.floatSensB(euribor.dRateLogL(paydate_prev), euribor.dRateLogL(paydate), d) * phi.sensB(); //float paym sens due to current point discount
							Gradient[count - gradientCut] += instrument.notional*instrument.spreadFloat *yFracFloat *phi.sensB(); //spread sens
							Gradient[count - gradientCut] += -instrument.notional*instrument.fixedRate *yFracFloat *phi.sensB(); //fixed sens, for futures only

						}
						//py::print("an upper lower bound located for quarter:", i);
						break;
						 
					}
					else
					{
						if (an_unfortunate_index == points_values.size())
						{
							//discount is assumed to be this one:
							py::print("extrapolation region!");
							py::print("to be coded yet!");
							break;
						}
						else
						{
							foundFirstZeroRate = phi.zeroR;
							foundfirst_cpdate = k->first;
							foundfirst_cpvalue = k->second;
							count += 1;
							//py::print("found a lower point, and here foundfirst_cpdate foundfirst_cpvalue at count ", count, foundfirst_cpdate, foundfirst_cpvalue); 
						}
					}
				}
				 
				if (paydate <= k->first && foundfirst_cpdate == 0)
				{
					//this part is not at work at the moment, a curve will always start from today
					break;
				}
				if (paydate > k->first && foundfirst_cpdate == 0)
				{
					foundFirstZeroRate = 1000;
					foundfirst_cpdate = k->first;
					foundfirst_cpvalue = k->second;
					count += 1;
					//py::print("found a lower point, and here foundfirst_cpdate foundfirst_cpvalue", foundfirst_cpdate,foundfirst_cpvalue);
				}
			}
			 
			//py::print("done with a pay date (went over all instruments), here I have d_prev, d as ", d_prev, d);
			//py::print("if d came as none check it out! from the prints above...");
			//you are done with loop over curve points, record, calculate and do to next swap paydate
			//fixed leg calculations are not necessary for float-Float
			
			 
			pvFloatLeg += floatSens.floatPV(euribor.dRateLogL(paydate_prev), euribor.dRateLogL(paydate), d);  //no time frac here, be careful in future implementations
			pvFloatSpreadLeg  += instrument.spreadFloat * yFracFloat * floatSens.floatSpreadPV(euribor.dRateLogL(paydate), d);
			pvFixedLeg += instrument.fixedRate * yFracFloat * floatSens.fixedPV(euribor.dRateLogL(paydate), d); //this is only for futures instrument
			 
			 
		}
		double pvFloatLegFinExch = 0.0;
		if (finalReturn == true)
		{
			Gradient[indexForFinalExchnage - 1] += instrument.notional; //sensitivity due to the final exchange
			pvFloatLegFinExch = instrument.notional  *dfinal * spotFxRate;
		}
		   

		//this is a colletarelized libor swap with 0.0 fixed rate
		//day count fixed of instrument above is used for the daycount for Libor in the swap below
		SwapInstrument anInst(instrument.floatTenor, instrument.fixedTenor, "", "",
			0.0, instrument.numOfPeriods / alpha, instrument.spotDate, instrument.dayCFixed, instrument.dayCFixed,
			0.0, instrument.edControl, 0.9);

		double pvdomFloatLeg = instrumentPrice(anInst, domLibor, domOIS, finalReturn);
	
		//py::print("Gradient");
		//py::print(Gradient);
		//py::print("instrument.numOfPeriods", instrument.numOfPeriods);
		//py::print(instrument.curvePoint, "pvFixedLeg should be zero unless a futures point", pvFixedLeg);
		//py::print(instrument.curvePoint, "libor 3m Leg PV ", pvdomFloatLeg );
		//py::print(instrument.curvePoint, " 3M EuriborpvFloatLeg  ", pvFloatLeg);
		//py::print(instrument.curvePoint, " 3M Euribor pvFloatSpreadLeg ", pvFloatSpreadLeg);
		//py::print(instrument.curvePoint, " 3M Euribor pvFloatLegFinExch  (dfinal)", pvFloatLegFinExch, dfinal);
		//py::print(instrument.curvePoint, "CCS PV  ", pvFloatLeg + pvFloatSpreadLeg - pvFixedLeg + pvFloatLegFinExch - pvdomFloatLeg);
		//py::print(instrument.curvePoint, "CCS PV euribor vs libor ", pvFloatLeg + pvFloatSpreadLeg - pvFixedLeg + pvFloatLegFinExch, pvdomFloatLeg);
		 
		std::tuple<double, std::vector<double>> result( pvFloatLeg +pvFloatSpreadLeg - pvFixedLeg+ pvFloatLegFinExch - pvdomFloatLeg , Gradient);
		return result;
		
	}//

	std::tuple<double, std::vector<double>>  pGCurveCCSCubS(std::map<int, double, std::less<int>>& points_values, int numOfFixings,
		bool finalReturn, double spotFxRate, int curveSpot, Curve &euribor, Curve &domLibor, Curve &domOIS, const GoodInterp& phi) const
	{
		//initialize day count calculators
		DayCountCalculator dayCCFixed = DayCountCalculator(instrument.dayCFixed);
		DayCountCalculator dayCCFloat = DayCountCalculator(instrument.dayCFloat);
		double yFracFixed = 0;
		double yFracFloat = 0; 
		double pvFixedLeg = 0.0;
		double pvFloatLeg = 0.0;
		double pvFloatSpreadLeg = 0.0;

		int paydate_prev = instrument.spotDate;
		int paydate = instrument.spotDate;

		std::vector<double> Gradient(points_values.size() - numOfFixings);

		//this code is only for tenor basis swaps or cross c
		assert(instrument.fixedRate == 0.0);
		//this code assumes fixed (or the unsolved for float) tenor is at least as large as the float tenor
		assert(instrument.fixedTenor >= instrument.floatTenor);
		int alpha = 1;// instrument.fixedTenor / instrument.floatTenor; // 1 if 3M-3M or 2 if 6M-3M etc.


		double dfinal = 0;
		int indexForFinalExchange =0;
		 
		
		for (int i = 1, e = instrument.numOfPeriods; i <= e; ++i)//i is used for parity; parity is for 3M-6M setup
		{
			//shift paydate while keeping a prev for the next period
			int auxpaydate = paydate;
			if (instrument.floatTenor == 3 || instrument.fixedTenor == 3)
			{
				if (instrument.numOfPeriods == 1)
					paydate = shiftByThreeMonths(paydate, true);
				else
					paydate = shiftByThreeMonths(paydate, false);
			}

			else {
				if (instrument.floatTenor == 6 || instrument.fixedTenor == 6)
					paydate = shiftBySixMonths(paydate);
				else py::print("this tenor type is not coded for PriceJacobian!");
			}
			//py::print("cubic gradient called in PG, for paydate:", paydate);
			std::vector<double> cubicGradient(phi.Interp(paydate));
			//py::print("cubic grad in pGCubicS ", cubicGradient);

			yFracFixed = dayCCFixed.yFrac(paydate_prev, paydate);
			yFracFloat = dayCCFloat.yFrac(auxpaydate, paydate);
			paydate_prev = auxpaydate;

			InstrSens floatSens = InstrSens(instrument);
			double d = 0;
			int count = -1; //should be the index of the least upper bound point: 
							//e.g. if: last fixingdate < paydate < first non-fixing date, count = 0
			int gradientCut = numOfFixings - 1;
			int an_unfortunate_index = 0;

			double euribord_prev =  euribor.dRateLogL(paydate_prev);
			double euribord =   euribor.dRateLogL(paydate);
			//py::print("above two calls to the curve");
			
			//test code, delete!
			//int ii = 0;
			//for (auto k = points_values.begin(), e = points_values.end(); k != e; ++k)
			//{
			//	d += k->second * cubicGradient[i];
			//	++ii;
			//}
			//-test code
			if (d < 0.01)
			{
				//py::print("d<0.01");
				//py::print("d from CubS: ", d);
				//py::print("interpolate on points", points_values);
				//py::print("cubic grad in pGCubicS ", cubicGradient);// [cubicGradient.size() - 1]);
			}

			for (auto k = points_values.begin(), e = points_values.end(); k != e; ++k)
			{
				
				an_unfortunate_index += 1;
				if (count >= 0)
				{
					//py::print("count", count);
					if (paydate <= k->first)
					{
						int u = 0;
						for (auto k = points_values.begin(), e = points_values.end(); k != e; ++k)
						{
							d += k->second *  cubicGradient[u];
							++u;
						}


						if (i == instrument.numOfPeriods)
						{
							dfinal = d;
							indexForFinalExchange = an_unfortunate_index - numOfFixings;
						}
						
						for (unsigned int f = 0; f < points_values.size() - numOfFixings; ++f)
						{
							Gradient[f] += floatSens.floatSensB(euribord_prev, euribord, d) *cubicGradient[gradientCut + 1 + f]; //float paym sens due to current point discount (note: N (d_prev -d) = N( d_prev - a d_k-1 - (1-a)d_k))
							Gradient[f] += instrument.notional*instrument.spreadFloat *yFracFloat *cubicGradient[gradientCut + 1 + f]; //spread sens
							Gradient[f] += -instrument.notional*instrument.fixedRate *yFracFloat *cubicGradient[gradientCut + 1 + f]; //fixed sens, for futures only
							//py::print("cubicGradient[gradientCut + 1 + f]", cubicGradient[gradientCut + 1 + f]);
						}
						//py::print("an upper lower bound located for quarter:", i);
						break;
					}
					else
					{
						if (an_unfortunate_index == points_values.size())
						{
							//discount is assumed to be this one:
							py::print("extrapolation region!");
							py::print("to be coded yet!");
							break;
						}
						else
						{
							count += 1;
							//py::print("found a lower point, and here foundfirst_cpdate foundfirst_cpvalue at count ", count, foundfirst_cpdate, foundfirst_cpvalue); 
						}
					}
				}

				if (paydate <= k->first && count == -1)
				{
					//this part is not at work at the moment, a curve will always start from today
					break;
				}
				if (paydate > k->first && count == -1)
				{
					count += 1;
					//py::print("found a lower point, and here foundfirst_cpdate foundfirst_cpvalue", foundfirst_cpdate,foundfirst_cpvalue);
				}
				 
			}

			//py::print("done with a pay date (went over all instruments), here I have d_prev, d as ", d_prev, d);
			//py::print("if d came as none check it out! from the prints above...");
			//you are done with loop over curve points, record, calculate and do to next swap paydate
			//fixed leg calculations are not necessary for float-Float
			pvFloatLeg += floatSens.floatPV(euribord_prev, euribord, d);  //no time frac here, be careful in future implementations
			pvFloatSpreadLeg += instrument.spreadFloat * yFracFloat * floatSens.floatSpreadPV(euribord, d);
			pvFixedLeg += instrument.fixedRate * yFracFloat * floatSens.fixedPV(euribord, d); //this is only for futures instrument
		}
		
		 
		double pvFloatLegFinExch = 0.0;
		if (finalReturn == true)
		{
			Gradient[indexForFinalExchange -1] += instrument.notional; //sensitivity due to the final exchange
			pvFloatLegFinExch = instrument.notional  *dfinal * spotFxRate;
		}

		 
		//this is a colletarelized libor swap with 0.0 fixed rate
		//day count fixed of instrument above is used for the daycount for Libor in the swap below
		SwapInstrument anInst(instrument.floatTenor, instrument.fixedTenor, "", "",
			0.0, instrument.numOfPeriods / alpha, instrument.spotDate, instrument.dayCFixed, instrument.dayCFixed,
			0.0, instrument.edControl, 0.9);

		double pvdomFloatLeg = instrumentPrice(anInst, domLibor, domOIS, finalReturn);
		 
		//py::print("Gradient");
		//py::print(Gradient);
		//py::print("instrument.numOfPeriods", instrument.numOfPeriods);
		//py::print(instrument.curvePoint, "pvFixedLeg should be zero unless a futures point", pvFixedLeg);
		//py::print(instrument.curvePoint, "libor 3m Leg PV ", pvdomFloatLeg );
		//py::print(instrument.curvePoint, " 3M EuriborpvFloatLeg  ", pvFloatLeg);
		//py::print(instrument.curvePoint, " 3M Euribor pvFloatSpreadLeg ", pvFloatSpreadLeg);
		//py::print(instrument.curvePoint, " 3M Euribor pvFloatLegFinExch  (dfinal)", pvFloatLegFinExch, dfinal);
		//py::print(instrument.curvePoint, "CCS PV  ", pvFloatLeg + pvFloatSpreadLeg - pvFixedLeg + pvFloatLegFinExch - pvdomFloatLeg);
		//py::print(instrument.curvePoint, "CCS PV euribor vs libor ", pvFloatLeg + pvFloatSpreadLeg - pvFixedLeg + pvFloatLegFinExch, pvdomFloatLeg);

		std::tuple<double, std::vector<double>> result(pvFloatLeg + pvFloatSpreadLeg - pvFixedLeg + pvFloatLegFinExch - pvdomFloatLeg, Gradient);
		return result;
		 
	}//


	//THIS GUY SELECTS FROM THE ONES ABOVE
	std::tuple<double, std::vector<double>> defaultPG(std::map<int, double, std::less<int>>& points_values, int numOfFixings, bool finalReturn, double spotFxRate, 
		int curveSpot, Curve &OIS, Curve &domLibor, Curve &domOIS, GoodInterp& phi) const
	{
		//py::print("instrument paym type", instrument.paymType);
		if (instrument.numOfPeriods == 1 && instrument.paymType != "flofloCCS" )
		{
			return pGCurveFutures(points_values, numOfFixings, curveSpot, OIS);
		}
		else
		{
			if(instrument.paymType == "flofix")  return pGCurveSwap(points_values, numOfFixings, curveSpot, OIS);
			else
			{
				if (instrument.paymType == "floflo")  return pGCurveTBS(points_values, numOfFixings, finalReturn, spotFxRate, curveSpot, OIS, domLibor, domOIS);
				else
				{
					if (instrument.paymType == "flofloCCS")  return pGCurveCCS(points_values, numOfFixings, finalReturn, spotFxRate, curveSpot, OIS, domLibor, domOIS);
					else 
					{
						if (instrument.paymType == "flofloCCSCubS")  return pGCurveCCSCubS(points_values, numOfFixings, finalReturn, spotFxRate, curveSpot, OIS, domLibor, domOIS, phi);
						else py::print("not a known type!");
					}
				}
			}
		}
 
	}



 
};

