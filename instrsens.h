#pragma once
#include "instruments.h" 

// assuming float payment is determined by two discount numbers: dprev and d.
//dprev is for the earlier time point
class InstrSens
{
public:
	SwapInstrument &instrument;
	InstrSens(SwapInstrument &instrument) : instrument(instrument) {	};

	double floatSensA(double dprev, double d, double dOIS)
	{
		if (instrument.swapType == "plain") return instrument.notional;
		if (instrument.swapType == "plainFwdOnly") return instrument.notional * dOIS/ d;
		if (instrument.swapType == "plainDiscOnly") return instrument.notional * (dprev / d -1);
		if (instrument.swapType == "fedf") return instrument.notional*d /dprev;
	}
	double floatSensB(double dprev, double d, double dOIS)
	{
		if (instrument.swapType == "plain") return -instrument.notional;
		if (instrument.swapType == "plainFwdOnly") return -instrument.notional * dOIS * dprev/ (d*d);
		if (instrument.swapType == "plainDiscOnly") return instrument.notional * (dprev / d - 1);
		if (instrument.swapType == "fedf") return instrument.notional*(log(dprev/d)-1);
	}

	double fixedSensA()
	{
		if (instrument.swapType == "plain") return instrument.notional;
		if (instrument.swapType == "plainFwdOnly") return 0.0;   
		if (instrument.swapType == "plainDiscOnly") return instrument.notional;
		if (instrument.swapType == "fedf") return instrument.notional;
	}
	double fixedSensB()
	{
		if (instrument.swapType == "plain") return instrument.notional;
		if (instrument.swapType == "plainFwdOnly") return 0.0;  
		if (instrument.swapType == "plainDiscOnly") return instrument.notional;
		if (instrument.swapType == "fedf") return instrument.notional;
	}

	double floatPV(double dprev, double d, double dOIS)
	{
		if (instrument.swapType == "plain") return (dprev / d - 1) * instrument.notional *d ;
		if (instrument.swapType == "plainFwdOnly") return (dprev / d - 1) * instrument.notional *dOIS;
		if (instrument.swapType == "plainDiscOnly") return (dprev / d - 1) * instrument.notional *dOIS;
		if (instrument.swapType == "fedf") return log(dprev / d) * instrument.notional *d;
	}

	double fixedPV(double d, double dOIS)
	{
		if (instrument.swapType == "plain") return  instrument.notional *d;
		if (instrument.swapType == "plainFwdOnly") return  instrument.notional *dOIS;
		if (instrument.swapType == "plainDiscOnly") return instrument.notional *dOIS;
		if (instrument.swapType == "fedf") return  instrument.notional *d;
	}

	double floatSpreadPV(double d, double dOIS)
	{
		if (instrument.swapType == "plain") return  instrument.notional *d;
		if (instrument.swapType == "plainFwdOnly") return  instrument.notional *dOIS;
		if (instrument.swapType == "plainDiscOnly") return instrument.notional *dOIS;
		if (instrument.swapType == "fedf") return  instrument.notional *d;
	}
};
