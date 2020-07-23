/*
d, d_prev for swap pay dates
foundfirstdisc etc. for curve points

*/
#pragma once
#include "dates.h" 

/*
OIS swaps: Fixed: 3Month Act/360, Float: 3Month, Act/360
Libor swaps: Fixed: 6Month 30/360, Float: 3Month, Act/360
*/

//LATER you might split swaps into legs
struct SwapInstrument //USD LIBOR
{
	int floatTenor; //"3M" 6M", "3M3M"
	int fixedTenor; //"3M" "6M", "3M3M"
	std::string swapType; //plain, plainFwdOnly, "fedf", 
	std::string paymType; //flofix, flofloCCS,flofloTBS 
	std::string tenorType = "6M12M"; //flofix, flofloCCS,flofloTBS

	int maturityPoint = 0; //for 30Year wap: 30*364
	int numOfPeriods; //enter num of periods
	double notional = 1'000'000;
	double fixedRate = 1;
	int spotDate;
	std::string dayCFixed;
	std::string dayCFloat;
		
	int curvePoint;
	double curvePointV = 1.1;
	
	double spreadFloat;

	bool edControl = false;
	 

	explicit SwapInstrument(int floatTenor, int fixedTenor, std::string swapType, std::string paymType, 
		double fixedRate,  int numOfPeriods, int spotDate, std::string dayCFixed, std::string dayCFloat,
		double spreadFloat, bool edControl, double curvePointV ):
		floatTenor(floatTenor), fixedTenor(fixedTenor), swapType(swapType), paymType(paymType), 
		fixedRate(fixedRate), numOfPeriods(numOfPeriods), spotDate(spotDate), dayCFixed(dayCFixed), dayCFloat(dayCFloat),  
		spreadFloat(spreadFloat), edControl(edControl), curvePointV(curvePointV)
	{
		//numof3mPeriods = maturityPoint / 91;
		//py::print("swap created: fixedRate: ", fixedRate, " maturity days: ", maturityPoint);
		int numOfYears;
		if (floatTenor == 3 || fixedTenor == 3) numOfYears = numOfPeriods / 4;
		else 
		{
			if (floatTenor == 6 || fixedTenor == 6) numOfYears = numOfPeriods / 2;
			else 
			{
				if (floatTenor == 1 || fixedTenor == 1) numOfYears = numOfPeriods / 12;
				py::print("this tenor type is not coded!");
			}
		}
		 

		maturityPoint = spotDate;
		if(edControl == false)
		{
			//py::print("not a futures device");
			for (int i = 0; i < numOfYears; ++i)
			{
				maturityPoint = shiftByOneYear(maturityPoint);
			}
		}
		else 
		{ //maturity point is important as it is the curve point for the instrument.
			//I should separate these two concepts!!!
			//py::print("a futures device");
			//maturityPoint = maturityPoint; 
			
			int yr = findTheYear(spotDate);
			int mo = findTheMonth(spotDate);
			if ((mo + 3) % 12 < mo) yr += 1;
			maturityPoint = findIMMDateGivenMonthYear((mo + 3) %12, yr);
			
		}
		curvePoint = maturityPoint;
	}
		
	//this copy constructer is trouble
	//sometimes i use this and solve problems, sometimes I turn it off and solve problems!!!
};

struct FuturesInstrument
{
	std::string tenor;
	int maturityPoint;
	int date;
	double notional = 7'000'000;
	double fixedrate = 0.0;
	int spotDate;
	//FuturesInstrument(std::string tenor, int date) : tenor(tenor), date(date) {}
};

struct FixingInstrument
{
	double curvePointV;
	int curvePoint;
	FixingInstrument(int curvePoint, double curvePointV) : curvePoint(curvePoint), curvePointV(curvePointV) {}

	FixingInstrument(const FixingInstrument& x)
	{
		curvePointV = x.curvePointV;
		curvePoint = x.curvePoint;

	}
};


struct LiborInstrument
{
	std::string tenor;
	int maturityPoint;
	int numofPeriods;
	int date;
	int spotDate;
	double notional = 1000000;
	double fixedRate = 0.025;
	LiborInstrument(std::string tenor, int date) : tenor(tenor), date(date) {}
	LiborInstrument(std::string tenor, int date, double notional) : tenor(tenor), date(date), notional(notional) {}
	LiborInstrument(int maturityPoint, double fixedRate) : maturityPoint(maturityPoint), fixedRate(fixedRate) {}
	int getmaturityPoint() const
	{
		return maturityPoint;
	}

	int getnumPeriods() const
	{
		return numofPeriods;
	}
};