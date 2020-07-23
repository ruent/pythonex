#pragma once
#include "instruments.h"
#include "curve.h"
#include "dates.h"
//this is to check that curve calibration and and curve interpolation are consistent.
//this can be verified by pricing the instruments using the calibrated curve

//this is also for CCS domestic leg PV calculation

double instrumentPrice(SwapInstrument anInst, Curve &fwdCurve, Curve &discCurve, bool finalReturn)
{
	double pvFloatLeg = 0.0;
	int paydate_prev = anInst.spotDate;
	int paydate = paydate_prev;
	if (anInst.numOfPeriods>1)
	{
		for (int i = 1, e = anInst.numOfPeriods; i <= e; ++i)
		{
			if (anInst.floatTenor == 3 || anInst.fixedTenor == 3) paydate = shiftByThreeMonths(paydate_prev, false);
			else {
				if (anInst.floatTenor == 6 || anInst.fixedTenor == 6) paydate = shiftBySixMonths(paydate);
				else py::print("this tenor type is not coded for pricer.h!");
			}		 
			pvFloatLeg += anInst.notional * (fwdCurve.dRateLogL(paydate_prev)/ fwdCurve.dRateLogL(paydate) -1 )* discCurve.dRateLogL(paydate); //no time frac here, be careful in future implementations
			paydate_prev = paydate;
		}
		double dfinal = discCurve.dRateLogL(paydate);
		pvFloatLeg += finalReturn * anInst.notional * dfinal;
		//py::print("num of years, final disc", anInst.numOfPeriods/4, dfinal);
	}
	else
	{
		if (anInst.floatTenor == 3 && anInst.fixedTenor == 3) paydate = shiftByThreeMonths(paydate_prev, true);
		else py::print("not ready for such tenors!");
		pvFloatLeg += anInst.notional * (fwdCurve.dRateLogL(paydate_prev) / fwdCurve.dRateLogL(paydate) - 1)* discCurve.dRateLogL(paydate); //no 
	}
	return pvFloatLeg;
}