#pragma once
#include <cassert>
namespace py = pybind11;
#include "dates.h"
#include "instruments.h"
#include "interp.h"
#include "instrsens.h"
#include "curve.h"

//forward declaration
//class Curve;


template <typename T, typename TInterp>
class InstrumentPG
{
public:
	T &instrument;
	InstrumentPG(T &instrument) : instrument(instrument) { //py::print("yeah, instrumentPG created!"); 
	};  
/*
	std::tuple<double, std::vector<double>>  pGLinRatesSingleCurveSMTM(std::map<int, double, std::less<int>>& points_values) const
	{
		double pvFixedLeg = 0.0;
		double pvFloatLeg = 0.0;
		std::vector<double> Gradient(points_values.size()); // this should have as many entries as points_values
		//py::print("xxxxxxxxxxxxxx, points_values size", points_values.size());
		//py::print("getnumPeriods", instrument.getnumPeriods());
		double d_prev = 1; //initial discount 
		int paydate = instrument.spotDate;
		for (int i = 1, e = instrument.numofPeriods; i <= e; ++i)
		{
			paydate = shiftByThreeMonths(paydate);
			//py::print("asset number:", i, "numofperiods ", instrument.numofPeriods, "maturity point:", instrument.maturityPoint);
			//py::print("getnumPeriods", instrument.getnumPeriods());
			//py::print("getmaturityPoint", instrument.getmaturityPoint());
			double a = 0;
			double d = 0;
			int foundfirstday = 0;
			double foundfirstdisc = 0.0;
			int count = 0;
			double a_old = 0.0;

			int count_old = 0;
			int an_unfortunate_index = 0;
			for (auto k = points_values.begin(), e = points_values.end(); k != e; ++k)
			{
				an_unfortunate_index += 1;
				if (foundfirstday != 0)
				{
					if (paydate <= k->first)
					{
						a = (k->first - paydate) / double(k->first - foundfirstday);
						d = a * foundfirstdisc + (1 - a) * k->second;
						assert(count > 0);

						if (count_old > 0)Gradient[count_old - 1] += a_old * instrument.notional; //float paym sens due to previous point discount
						Gradient[count_old] += (1 - a_old) * instrument.notional; //float paym sens due to previous point discount
						Gradient[count - 1] += -a * instrument.notional; //float paym sens due to current point discount (note: N (d_prev -d) = N( d_prev - a d_k-1 - (1-a)d_k))
						Gradient[count] += -(1 - a) * instrument.notional; //float paym sens due to current point discount
						count_old = count;
						a_old = a;
						if (i % 2 == 0) Gradient[count - 1] += -instrument.fixedRate *0.5 * a * instrument.notional; //fix paym sens due to previous point discount
						if (i % 2 == 0) Gradient[count] += -instrument.fixedRate *0.5 * (1 - a)* instrument.notional; //fix paym sens due to current point discount
						//py::print("an upper lower bound located for quarter:", i);
						break;

					}
					else
					{
						if (an_unfortunate_index == points_values.size())
						{
							d = d_prev; //no float payment
							if (i % 2 == 0) Gradient[points_values.size() - 1] += -instrument.fixedRate*0.5 * instrument.notional; //beyond last curve point
							break;
						}
						else
						{
							foundfirstday = k->first;
							foundfirstdisc = k->second;
							//py::print("lower bound should be updated:", i);
						}
					}
				}

				if (paydate <= k->first && foundfirstday == 0)
				{
					a = (k->first - paydate) / (double(k->first) - double(instrument.spotDate));
					d = a * 1 + (1 - a) *k->second;
					Gradient[count] += -(1 - a)* instrument.notional; //due to float payment
					if (i % 2 == 0) Gradient[count] += -instrument.fixedRate *0.5 *(1 - a)* instrument.notional; //due to fixed payment
					//py::print("smaller than the smallest:");
					break;
				}
				if (paydate > k->first && foundfirstday == 0)
				{
					foundfirstday = k->first;
					foundfirstdisc = k->second;
					//py::print("a candidate for lower bound:", i);
				}
				count += 1;
			}

			if (i % 2 == 0)
			{
				pvFixedLeg += instrument.fixedRate * 0.5 * instrument.notional*d;
				//year frac is fixed at 0.5
				//py::print(" fixed rate", instrument.fixedRate, "notio ", instrument.notional, "disc", d);
			}
			pvFloatLeg += (d_prev / d - 1) * instrument.notional *d; //no time frac here, be careful in future implementations
			d_prev = d;
		}
		std::tuple<double, std::vector<double>> result(-pvFixedLeg + pvFloatLeg, Gradient);
		return result;
	}//

	std::tuple<double, std::vector<double>>  pGLinRatesSingleCurveTMTM(std::map<int, double, std::less<int>>& points_values) const
	{
		double pvFixedLeg = 0.0;
		double pvFloatLeg = 0.0;
		std::vector<double> Gradient(points_values.size()); // this should have as many entries as points_values
		//py::print("xxxxxxxxxxxxxx, points_values size", points_values.size());
		//py::print("getnumPeriods", instrument.getnumPeriods());
		double d_prev = 1; //initial discount 

		//initialize day count calculators
		DayCountCalculator dayCCFixed = DayCountCalculator(instrument.dayCFixed);
		DayCountCalculator dayCCFloat = DayCountCalculator(instrument.dayCFloat);
		double yFracFixed = 0;
		double yFracFloat = 0;

		int paydate = instrument.spotDate;
		for (int i = 1, e = instrument.numofPeriods; i <= e; ++i)
		{
			//shift paydate while keeping a prev for the next period
			int auxpaydate = paydate;
			paydate = shiftByThreeMonths(paydate);
			yFracFixed = dayCCFixed.yFrac(auxpaydate, paydate);
			yFracFloat = dayCCFloat.yFrac(auxpaydate, paydate);
			//py::print("asset number:", i, "numofperiods ", instrument.numofPeriods, "maturity point:", instrument.maturityPoint);
			//py::print("getnumPeriods", instrument.getnumPeriods());
			//py::print("getmaturityPoint", instrument.getmaturityPoint());
			double a = 0;
			double d = 0;
			int foundfirstday = 0;
			double foundfirstdisc = 0.0;
			int count = 0;
			double a_old = 0.0;

			int count_old = 0;
			int an_unfortunate_index = 0;


			for (auto k = points_values.begin(), e = points_values.end(); k != e; ++k)
			{
				an_unfortunate_index += 1;
				if (foundfirstday != 0)
				{
					if (paydate <= k->first)
					{
						a = (k->first - paydate) / double(k->first - foundfirstday);
						d = a * foundfirstdisc + (1 - a) * k->second;
						assert(count > 0);

						if (count_old > 0)Gradient[count_old - 1] += a_old * instrument.notional; //float paym sens due to previous point discount
						Gradient[count_old] += (1 - a_old) * instrument.notional; //float paym sens due to previous point discount
						Gradient[count - 1] += -a * instrument.notional; //float paym sens due to current point discount (note: N (d_prev -d) = N( d_prev - a d_k-1 - (1-a)d_k))
						Gradient[count] += -(1 - a) * instrument.notional; //float paym sens due to current point discount
						count_old = count;
						a_old = a;
						Gradient[count - 1] += -instrument.fixedRate *yFracFixed * a * instrument.notional; //fix paym sens due to previous point discount
						Gradient[count] += -instrument.fixedRate *yFracFixed * (1 - a)* instrument.notional; //fix paym sens due to current point discount
						//py::print("an upper lower bound located for quarter:", i);
						break;

					}
					else
					{
						if (an_unfortunate_index == points_values.size())
						{
							d = d_prev; //no float payment
							Gradient[points_values.size() - 1] += -instrument.fixedRate*yFracFixed * instrument.notional; //beyond last curve point
							break;
						}
						else
						{
							foundfirstday = k->first;
							foundfirstdisc = k->second;
							//py::print("lower bound should be updated:", i);
						}
					}
				}

				if (paydate <= k->first && foundfirstday == 0)
				{
					a = (k->first - paydate) / (double(k->first) - double(instrument.spotDate));
					d = a * 1 + (1 - a) *k->second;
					Gradient[count] += -(1 - a)* instrument.notional; //due to float payment
					Gradient[count] += -instrument.fixedRate *yFracFixed *(1 - a)* instrument.notional; //due to fixed payment
					//py::print("smaller than the smallest:");
					break;
				}
				if (paydate > k->first && foundfirstday == 0)
				{
					foundfirstday = k->first;
					foundfirstdisc = k->second;
					//py::print("a candidate for lower bound:", i);
				}
				count += 1;
			}


			pvFixedLeg += instrument.fixedRate * yFracFixed * instrument.notional*d;
			pvFloatLeg += (d_prev / d - 1) * instrument.notional *d; //no time frac here, be careful in future implementations
			d_prev = d;
		}
		std::tuple<double, std::vector<double>> result(-pvFixedLeg + pvFloatLeg, Gradient);
		return result;
	}//
*/

	/*std::tuple<double, std::vector<double>>  pGLinZerosSingleCurve(std::map<int, double, std::less<int>>& points_values, int numOfFixings, int curveSpot) const
	{
		double pvFixedLeg = 0.0;
		double pvFloatLeg = 0.0;

		double d_prev = 1; //initial discount //I use this in fwd rate calc; I don't know if dprev/d sacrifices any precision in 3month rate
		double a_prev = 0.0;
		int count_prev = 0; //this will be the least upper bound for the previous payment date
		double curvepoilow_prev = 1;
		double curvepoiup_prev = 1;
		int curvepoilowd_prev;;
		int curvepoiupd_prev;
		//initialize day count calculators
		DayCountCalculator dayCCFixed = DayCountCalculator(instrument.dayCFixed);
		DayCountCalculator dayCCFloat = DayCountCalculator(instrument.dayCFloat);
		double yFracFixed = 0;
		double yFracFloat = 0;
		int paydate_prev = instrument.spotDate;
		int paydate = instrument.spotDate;

		std::vector<double> Gradient(points_values.size() - numOfFixings);

		//i is used for parity; parity is for 3M-6M setup
		for (int i = 1, e = instrument.numofPeriods; i <= e; ++i)
		{

			//shift paydate while keeping a prev for the next period
			int auxpaydate = paydate;
			if (instrument.tenorType == "3M6M" || instrument.tenorType == "6M3M" || instrument.tenorType == "3M3M")
				paydate = shiftByThreeMonths(paydate, false);
			else {
				if (instrument.tenorType == "12M6M" || instrument.tenorType == "6M12M" || instrument.tenorType == "6M6M")
					paydate = shiftBySixMonths(paydate);
				else py::print("trouble!");
			}

			int alpha = 1;
			if (instrument.tenorType == "3M6M" || instrument.tenorType == "6M3M" || instrument.tenorType == "6M12M" || instrument.tenorType == "12M6M")
				alpha = 2;
					   			 
			yFracFixed = dayCCFixed.yFrac(paydate_prev, paydate);
			yFracFloat = dayCCFloat.yFrac(auxpaydate, paydate);
			paydate_prev = auxpaydate;
			 
			InstrSens floatSens = InstrSens(instrument);
			
			//py::print("asset number:", i, "numofperiods ", instrument.numofPeriods, "maturity point:", instrument.maturityPoint);
			//py::print("getnumPeriods", instrument.getnumPeriods());
			//py::print("getmaturityPoint", instrument.getmaturityPoint());
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
				 //
				//if (an_unfortunate_index < numOfFixings)
				//{
				//	//count += 1;
				//	foundfirst_cpdate = k->first;
				//	foundfirst_cpvalue = k->second;
				//	continue;
				//}
				
				
				//Example: There are 22 fixings, and 12 swaps after them
				//gardientCut = 21
				//first swap first payment has count = 21
				//sensitivity tpo first swap will be the first (0) enry of gradient vectors
				//this happens if paydate is less than firstswap, larger then last fixing
				//and this is when count is at least 21
				//if count=21 there shoud be sens to instrument 22, or swap 0 (0 = count - gradientCut)
			//	if count = 22 there should be sensitivity to inst 22 and 23, i.e. swaps 0 and 1

				//if you price swap 1, then last 4 payments will have count 22
				

				if (foundfirst_cpdate != 0)
				{
					TInterp phi = TInterp(foundfirst_cpdate, k->first, foundfirst_cpvalue, k->second, paydate, curveSpot);
					TInterp phi_prev = TInterp(curvepoilowd_prev, curvepoiupd_prev, curvepoilow_prev, curvepoiup_prev, paydate_prev, curveSpot, a_prev, d_prev);

					//LogLinInterp phi = LogLinInterp(foundfirst_cpdate, k->first, foundfirst_cpvalue, k->second, paydate);
					//LogLinInterp phi_prev = LogLinInterp(curvepoilowd_prev, curvepoiupd_prev, curvepoilow_prev, curvepoiup_prev, a_prev, d_prev);
					//LinZInterp phi = LinZInterp(foundfirst_cpdate, k->first, foundfirst_cpvalue, k->second, paydate, curveSpot);
					//LinZInterp phi_prev = LinZInterp(curvepoilowd_prev, curvepoiupd_prev, curvepoilow_prev, curvepoiup_prev, paydate_prev, curveSpot, a_prev, d_prev);
					
					if (paydate <= k->first)
					{
						
						a = phi.a;
						d = phi.d;
						//Executive summary: below to lines for risk in N(dprev-d) due to d_prev; d_prev changes with curve[count_old] and curve[count_old-1]
						//different from linear interp you also need to bring two curve points with yourself here
						if (count >= gradientCut)
						{
							//py::print("pay period : ", i, " data: ", phi.a, phi.d, phi.sensA(), phi.sensB());
							//py::print("pay period : ", i, " data_prev: ", phi_prev.a, phi_prev.d, phi_prev.sensA(), phi_prev.sensB());
							if (count_prev > gradientCut) {
								Gradient[count_prev - gradientCut - 1] +=  floatSens.floatSensA(curvepoilow_prev, curvepoiup_prev) *phi_prev.sensA();
								//py::print("1", instrument.notional *phi_prev.sensA());
							}
							if (count_prev >= gradientCut) {
								Gradient[count_prev - gradientCut] += floatSens.floatSensA(curvepoilow_prev, curvepoiup_prev) *phi_prev.sensB();
								//py::print("2", instrument.notional *phi_prev.sensB());
							}
							//Executive summary: below two lines for risk in floatPV, N(dprev-d), due to d; d changes with curve[count] and curve[count-1]
							if (count > gradientCut) {
								Gradient[count - gradientCut - 1] += floatSens.floatSensB(foundfirst_cpvalue, k->second) *phi.sensA(); //float paym sens due to current point discount (note: N (d_prev -d) = N( d_prev - a d_k-1 - (1-a)d_k))
								//py::print("3", -instrument.notional *phi.sensA());
							}
							Gradient[count - gradientCut] += floatSens.floatSensB(foundfirst_cpvalue, k->second) * phi.sensB(); //float paym sens due to current point discount
							if (i % alpha == 0 && count > gradientCut) {
								Gradient[count - gradientCut - 1] += -instrument.fixedRate *yFracFixed * instrument.notional* phi.sensA(); //fix paym sens due to previous point discount
								//py::print("4", -instrument.fixedRate *yFracFixed * instrument.notional* phi.sensA());
							}
							if (i % alpha == 0) Gradient[count - gradientCut] += -instrument.fixedRate *yFracFixed * instrument.notional * phi.sensB(); //fix paym sens due to current point discount
						}
						//py::print("an upper lower bound located for quarter:", i);
						break;
					}
					else
					{
						if (an_unfortunate_index == points_values.size())
						{
							//discount is assumed to be this one:
							d = foundfirst_cpvalue * exp(-phi.zeroR * (paydate - foundfirst_cpdate));
							Gradient[count - 1 - gradientCut] += instrument.notional *  (1 - exp(-foundFirstZeroRate * (paydate - foundfirst_cpdate))); //due to float payment
							if (i % alpha == 0) Gradient[count - 1 - gradientCut] += -instrument.fixedRate*yFracFixed * instrument.notional *exp(-foundFirstZeroRate * (paydate - foundfirst_cpdate)); //beyond last curve point
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
					//these two lines were at the wrong place. they cost me a lot!
					curvepoilow_prev = foundfirst_cpvalue;
					curvepoiup_prev = k->second;
					curvepoilowd_prev = foundfirst_cpdate;
					curvepoiupd_prev = k->first;
				}

				if (paydate <= k->first && foundfirst_cpdate == 0)
				{
					//this part is not at work at the moment, a curve will always start from today
					break;
				}
				if (paydate > k->first && foundfirst_cpdate == 0)
				{
					curvepoilow_prev = foundfirst_cpvalue;
					curvepoiup_prev = k->second;
					curvepoilowd_prev = 10000; //this will be needed only for the pass thorugh first curve point; just a placeholder
					curvepoiupd_prev = 10000; //this will be needed only for the pass thorugh first curve point; just a placeholder
					foundFirstZeroRate = 1000;
					foundfirst_cpdate = k->first;
					foundfirst_cpvalue = k->second;
					count += 1;
				}

			}

			//you are done with loop over curve points, record, calculate and do to next swap paydate
			if (i % alpha == 0) pvFixedLeg += instrument.fixedRate * yFracFixed * instrument.notional*d;
			pvFloatLeg += floatSens.floatPV(d_prev, d); //no time frac here, be careful in future implementations
			d_prev = d; //d_prev is assigned right after an instrument is located on the curve points grid.
						//hence dprev is the latest discount, hence the latest discount for the latest paydate for the current instrument
			count_prev = count;
			a_prev = a;
		}

		std::tuple<double, std::vector<double>> result(-pvFixedLeg + pvFloatLeg, Gradient);
		return result;
	}//*/

	std::tuple<double, std::vector<double>>  pGLinZerosSingleCurveTMTM(std::map<int, double, std::less<int>>& points_values, int numOfFixings, int curveSpot) const
	{
		double pvFixedLeg = 0.0;
		double pvFloatLeg = 0.0;

		double d_prev = 1; //initial discount //I use this in fwd rate calc; I don't know if dprev/d sacrifices any precision in 3month rate
		double a_prev = 0.0;
		int count_prev = 0; //this will be the least upper bound for the previous payment date
		double curvepoilow_prev = 0.0;
		double curvepoiup_prev = 0.0;
		//initialize day count calculators
		DayCountCalculator dayCCFixed = DayCountCalculator(instrument.dayCFixed);
		DayCountCalculator dayCCFloat = DayCountCalculator(instrument.dayCFloat);
		double yFracFixed = 0;
		double yFracFloat = 0;
		int paydate_prev = instrument.spotDate;
		int paydate = instrument.spotDate;

		std::vector<double> Gradient(points_values.size() - numOfFixings);
				 
		for (int i = 1, e = instrument.numofPeriods; i <= e; ++i)
		{
			//shift paydate while keeping a prev for the next period
			int auxpaydate = paydate;

			if (instrument.tenorType == "3M3M") {
				py::print("call to duty, 3M3M!, calling 3M shifter");
				paydate = shiftByThreeMonths(paydate, false);
			}
			else{
				if (instrument.tenorType == "6M6M") paydate = shiftBySixMonths(paydate);
				else
				{
					;//py::print("trouble!");
				}
			}

			yFracFixed = dayCCFixed.yFrac(paydate_prev, paydate);
			yFracFloat = dayCCFloat.yFrac(auxpaydate, paydate);
			paydate_prev = auxpaydate;

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
					if (paydate <= k->first)
					{
						a = (k->first - paydate) / double(k->first - foundfirst_cpdate);
						d = exp(a*log(foundfirst_cpvalue) + (1 - a) * log(k->second));

						//Executive summary: below to lines for risk in N(dprev-d) due to d_prev; d_prev changes with curve[count_old] and curve[count_old-1]
						//different from linear interp you also need to bring two curve points with yourself here
						if (count >= gradientCut)
						{
							if (count_prev > gradientCut) {
								Gradient[count_prev - gradientCut - 1] += a_prev * instrument.notional *d_prev / curvepoilow_prev; //float paym sens due to d_prev
								//py::print("1",a_prev * instrument.notional *d_prev / curvepoilow_prev);
							}
							if (count_prev >= gradientCut) {
								Gradient[count_prev - gradientCut] += (1 - a_prev) * instrument.notional *d_prev / curvepoiup_prev; //float paym sens due to d_prev
								//py::print("2", (1 - a_prev) * instrument.notional *d_prev / curvepoiup_prev);
							}
							//Executive summary: below two lines for risk in floatPV, N(dprev-d), due to d; d changes with curve[count] and curve[count-1]
							if (count > gradientCut) {
								Gradient[count - gradientCut - 1] += -a * instrument.notional*d / foundfirst_cpvalue; //float paym sens due to current point discount (note: N (d_prev -d) = N( d_prev - a d_k-1 - (1-a)d_k))
								//py::print("3", -a * instrument.notional*d / foundfirst_cpvalue);
							}
							Gradient[count - gradientCut] += -(1 - a) * instrument.notional*d / k->second; //float paym sens due to current point discount
							if (count > gradientCut) {
								Gradient[count - gradientCut - 1] += -instrument.fixedRate *yFracFixed * a * instrument.notional*d / foundfirst_cpvalue; //fix paym sens due to previous point discount
								//py::print("4", -instrument.fixedRate *yFracFixed * a * instrument.notional*d / foundfirst_cpvalue);
							}
							Gradient[count - gradientCut] += -instrument.fixedRate *yFracFixed * (1 - a)* instrument.notional * d / k->second; //fix paym sens due to current point discount
						}
						//py::print("an upper lower bound located for quarter:", i);
						break;


					}
					else
					{
						if (an_unfortunate_index == points_values.size())
						{
							//discount is assumed to be this one:
							d = foundfirst_cpvalue * exp(-foundFirstZeroRate * (paydate - foundfirst_cpdate));
							Gradient[count - 1 - gradientCut] += instrument.notional *  (1 - exp(-foundFirstZeroRate * (paydate - foundfirst_cpdate))); //due to float payment
							Gradient[count - 1 - gradientCut] += -instrument.fixedRate*yFracFixed * instrument.notional *exp(-foundFirstZeroRate * (paydate - foundfirst_cpdate)); //beyond last curve point
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
					//these two linea where at the wrong place. they cost me a lot!
					curvepoilow_prev = foundfirst_cpvalue;
					curvepoiup_prev = k->second;
				}
				 
				
				if (paydate <= k->first && foundfirst_cpdate == 0)
				{
					//py::print("smaller than the smallest:");
					break;
				}
				if (paydate > k->first && foundfirst_cpdate == 0)
				{
					curvepoilow_prev = foundfirst_cpvalue;
					curvepoiup_prev = k->second;
					foundFirstZeroRate = 1000;
					foundfirst_cpdate = k->first;
					foundfirst_cpvalue = k->second;
					//py::print("Gradient size", Gradient.size());
					//py::print(Gradient[0]);
					count += 1;
				}

			}


			pvFixedLeg += instrument.fixedRate * yFracFixed * instrument.notional*d;
			pvFloatLeg += (d_prev / d - 1) * instrument.notional *d; //no time frac here, be careful in future implementations
			d_prev = d;
			count_prev = count;
			a_prev = a;
			//py::print("pvFixedLeg, pvFloatLeg ",pvFixedLeg, pvFloatLeg, instrument.fixedRate, yFracFixed, instrument.notional, d);
		}
		std::tuple<double, std::vector<double>> result(-pvFixedLeg + pvFloatLeg, Gradient);
		return result;
	}//

	std::tuple<double, std::vector<double>>  pGLinZerosSingleCurveFutures(std::map<int, double, std::less<int>>& points_values, int numOfFixings, int curveSpot) const
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
					if (count > 0) Gradient[count - 1] += instrument.notional; //float paym sens due to d_prev
					//py::print("count", count);
					//py::print("gradientNow", gradientNow);
					Gradient[count] += -instrument.notional; //float paym sens due to current point discount
					Gradient[count] += -instrument.fixedRate *yFracFixed * instrument.notional; //fix paym sens due to current point discount
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
		pvFixedLeg += instrument.fixedRate * yFracFixed * instrument.notional*d;
		//py::print("foundfirst_cpvalue", foundfirst_cpvalue);
		//py::print("d", d);
		pvFloatLeg += (foundfirst_cpvalue / d - 1) * instrument.notional *d; //no time frac here, be careful in future implementations
		//py::print("c++ pvFixedLeg: ", pvFixedLeg, "pvFloatLeg: ", pvFloatLeg, "yearF: ", yFracFixed);
		//py::print("price: ", -pvFixedLeg + pvFloatLeg, "Gradient!", Gradient);
		std::tuple<double, std::vector<double>> result(-pvFixedLeg + pvFloatLeg, Gradient);
		return result;
	}//
	
	//this one takes a discount curve and calculates another, as in USD OIS and USD LIBOR 
	std::tuple<double, std::vector<double>>  pGLinZerosCurveSequential(std::map<int, double, std::less<int>>& points_values, int numOfFixings, int curveSpot, Curve &oisCurve) const
	{
		double pvFixedLeg = 0.0;
		double pvFloatLeg = 0.0;
		 
		double d_prev = 1; //initial discount //I use this in fwd rate calc; I don't know if dprev/d sacrifices any precision in 3month rate
		double a_prev = 0.0;
		int count_prev = 0; //this will be the least upper bound for the previous payment date
		double curvepoilow_prev = 1;
		double curvepoiup_prev = 1;
		int curvepoilowd_prev;;
		int curvepoiupd_prev;
		//initialize day count calculators
		DayCountCalculator dayCCFixed = DayCountCalculator(instrument.dayCFixed);
		DayCountCalculator dayCCFloat = DayCountCalculator(instrument.dayCFloat);
		double yFracFixed = 0;
		double yFracFloat = 0;
		int paydate_prev = instrument.spotDate;
		int paydate = instrument.spotDate;

		std::vector<double> Gradient(points_values.size() - numOfFixings);

		//i is used for parity; parity is for 3M-6M setup
		for (int i = 1, e = instrument.numofPeriods; i <= e; ++i)
		{
			//shift paydate while keeping a prev for the next period
			int auxpaydate = paydate;
			if (instrument.tenorType == "3M6M" || instrument.tenorType == "6M3M" || instrument.tenorType == "3M3M")
				paydate = shiftByThreeMonths(paydate, false);
			else {
				if (instrument.tenorType == "12M6M" || instrument.tenorType == "6M12M" || instrument.tenorType == "6M6M")
					paydate = shiftBySixMonths(paydate);
				else py::print("trouble!");
			}

			int alpha = 1;
			if (instrument.tenorType == "3M6M" || instrument.tenorType == "6M3M" || instrument.tenorType == "6M12M" || instrument.tenorType == "12M6M")
				alpha = 2;

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
					TInterp phi = TInterp(foundfirst_cpdate, k->first, foundfirst_cpvalue, k->second, paydate, curveSpot);
					TInterp phi_prev = TInterp(curvepoilowd_prev, curvepoiupd_prev, curvepoilow_prev, curvepoiup_prev, paydate_prev, curveSpot, a_prev, d_prev);
					 

					if (paydate <= k->first)
					{
						 
						a = phi.a;
						d = phi.d;
						//Executive summary: below to lines for risk in N(dprev-d) due to d_prev; d_prev changes with curve[count_old] and curve[count_old-1]
						//different from linear interp you also need to bring two curve points with yourself here
						if (count >= gradientCut)
						{
							//py::print("pay period : ", i, " data: ", phi.a, phi.d, phi.sensA(), phi.sensB());
							//py::print("pay period : ", i, " data_prev: ", phi_prev.a, phi_prev.d, phi_prev.sensA(), phi_prev.sensB());
							if (count_prev > gradientCut) {
								Gradient[count_prev - gradientCut - 1] += floatSens.floatSensA(curvepoilow_prev, curvepoiup_prev, oisCurve.dRateLinZeros(paydate)) *phi_prev.sensA();
								//py::print("1", instrument.notional *phi_prev.sensA());
							}
							if (count_prev >= gradientCut) {
								Gradient[count_prev - gradientCut] += floatSens.floatSensA(curvepoilow_prev, curvepoiup_prev,oisCurve.dRateLinZeros(paydate)) *phi_prev.sensB();
								//py::print("2", instrument.notional *phi_prev.sensB());
							}
							//Executive summary: below two lines for risk in floatPV, N(dprev-d), due to d; d changes with curve[count] and curve[count-1]
							if (count > gradientCut) {
								Gradient[count - gradientCut - 1] += floatSens.floatSensB(foundfirst_cpvalue, k->second, oisCurve.dRateLinZeros(paydate)) *phi.sensA(); //float paym sens due to current point discount (note: N (d_prev -d) = N( d_prev - a d_k-1 - (1-a)d_k))
								//py::print("3", -instrument.notional *phi.sensA());
							}
							Gradient[count - gradientCut] += floatSens.floatSensB(foundfirst_cpvalue, k->second, oisCurve.dRateLinZeros(paydate)) * phi.sensB(); //float paym sens due to current point discount
							if (i % alpha == 0 && count > gradientCut) {
								Gradient[count - gradientCut - 1] += -instrument.fixedRate *yFracFixed * floatSens.fixedSensA()* phi.sensA(); //fix paym sens due to previous point discount
								//py::print("4", -instrument.fixedRate *yFracFixed * instrument.notional* phi.sensA());
							}
							if (i % alpha == 0) Gradient[count - gradientCut] += -instrument.fixedRate *yFracFixed * floatSens.fixedSensB() * phi.sensB(); //fix paym sens due to current point discount
						}
						//py::print("an upper lower bound located for quarter:", i);
						break;
					}
					else
					{
						if (an_unfortunate_index == points_values.size())
						{
							//discount is assumed to be this one:
							d = foundfirst_cpvalue * exp(-phi.zeroR * (paydate - foundfirst_cpdate));
							Gradient[count - 1 - gradientCut] += instrument.notional *  (1 - exp(-foundFirstZeroRate * (paydate - foundfirst_cpdate))); //due to float payment
							if (i % alpha == 0) Gradient[count - 1 - gradientCut] += -instrument.fixedRate*yFracFixed * instrument.notional *exp(-foundFirstZeroRate * (paydate - foundfirst_cpdate)); //beyond last curve point
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
					//these two lines were at the wrong place. they cost me a lot!
					curvepoilow_prev = foundfirst_cpvalue;
					curvepoiup_prev = k->second;
					curvepoilowd_prev = foundfirst_cpdate;
					curvepoiupd_prev = k->first;
				}

				if (paydate <= k->first && foundfirst_cpdate == 0)
				{
					//this part is not at work at the moment, a curve will always start from today
					break;
				}
				if (paydate > k->first && foundfirst_cpdate == 0)
				{
					curvepoilow_prev = foundfirst_cpvalue;
					curvepoiup_prev = k->second;
					curvepoilowd_prev = 10000; //this will be needed only for the pass thorugh first curve point; just a placeholder
					curvepoiupd_prev = 10000; //this will be needed only for the pass thorugh first curve point; just a placeholder
					foundFirstZeroRate = 1000;
					foundfirst_cpdate = k->first;
					foundfirst_cpvalue = k->second;
					count += 1;
				}

			}

			//you are done with loop over curve points, record, calculate and do to next swap paydate
			if (i % alpha == 0) pvFixedLeg += instrument.fixedRate * yFracFixed * floatSens.fixedPV(d, oisCurve.dRateLinZeros(paydate));
			pvFloatLeg += (floatSens.floatPV(d_prev, d, oisCurve.dRateLinZeros(paydate)) + instrument.spreadFloat * yFracFloat * floatSens.floatSpreadPV(d, oisCurve.dRateLinZeros(paydate))); //no time frac here, be careful in future implementations
			d_prev = d; //d_prev is assigned right after an instrument is located on the curve points grid.
						//hence dprev is the latest discount, hence the latest discount for the latest paydate for the current instrument
			count_prev = count;
			a_prev = a;
		}

		std::tuple<double, std::vector<double>> result(-pvFixedLeg + pvFloatLeg, Gradient);
		return result;
	}//


	//THIS GUY SELECTS FROM THE ONES ABOVE
	std::tuple<double, std::vector<double>> defaultPG(std::map<int, double, std::less<int>>& points_values, int numOfFixings, int curveSpot, Curve &OIS) const
	{
		if (instrument.numofPeriods == 1)
		{
			return pGLinZerosSingleCurveFutures(points_values, numOfFixings, curveSpot);
		}
		else
		{
			return pGLinZerosCurveSequential(points_values, numOfFixings, curveSpot, OIS);
		}

	}
 
};

