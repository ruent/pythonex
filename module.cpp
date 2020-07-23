#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
//#include <pybind11/stl_bind.h>
#include <Windows.h>
#include <cmath>
//#include <pybind11/stl_bind.h>


#include <string>
#include <map>
#include "instruments.h"
#include "instrumentpg.h"
#include "solver.h"
#include "curve.h"
#include "dates.h"
#include "interp.h"
#include "pricer.h"
 
//PYBIND11_MAKE_OPAQUE(std::map<int, double>);
namespace py = pybind11;


PYBIND11_MODULE(pythonex, m) {
	 
	m.def("isLeapYear", &isLeapYear, R"pbdoc(...)pbdoc");
	m.def("daysInAMonth", &daysInAMonth, R"pbdoc(...)pbdoc");
	m.def("findTheYear", &findTheYear, R"pbdoc(...)pbdoc");
	m.def("findTheMonth", &findTheMonth, R"pbdoc(...)pbdoc");
	m.def("findTheDayOfTheMonth", &findTheDayOfTheMonth, R"pbdoc(...)pbdoc");
	m.def("findTheMonth", &findTheMonth, R"pbdoc(...)pbdoc");
	m.def("shiftByThreeMonths", &shiftByThreeMonths, R"pbdoc(...)pbdoc");
	m.def("findTheFirstDayMonthYear", &findTheFirstDayMonthYear, R"pbdoc(...)pbdoc");
	m.def("findIMMDateGivenMonthYear", &findIMMDateGivenMonthYear, R"pbdoc(...)pbdoc");
	m.def("shiftByOneYear", &shiftByOneYear, R"pbdoc(...)pbdoc");
	m.def("shiftBySixMonths", &shiftBySixMonths, R"pbdoc(...)pbdoc");
	m.def("shiftByOneMonth", &shiftByOneMonth, R"pbdoc(...)pbdoc");
	
	m.def("instrumentPrice", &instrumentPrice, R"pbdoc(...)pbdoc");
	

	//py::class_<LogLinInterp>(m, "LogLinInterp");
	//py::class_<LinZInterp>(m, "LinZInterp");
	py::class_<GoodInterp> (m, "GoodInterp")
		.def(py::init<const std::map<int, double, std::less<int>>& >())
		//.def("finddDotDot", &GoodInterp::finddDotDot)
		.def("Interp", &GoodInterp::Interp);
	
	py::class_<Solver>(m, "SolverLogL")
		.def(py::init<>())
		.def("NewtonRaphson", &Solver::NewtonRaphson, R"pbdoc(applies Newton Raphson Method)pbdoc")
		.def("SteepestDescent", &Solver::SteepestDescent, R"pbdoc(applies steepest descent Method)pbdoc");

	//py::class_<Solver<LinZInterp>>(m, "SolverLinZ")
	//	.def(py::init<>())
	//	.def("NewtonRaphson", &Solver<LinZInterp>::NewtonRaphson, R"pbdoc(applies Newton Raphson Method)pbdoc");

	py::class_<SwapInstrument>(m, "SwapInstrument")
		.def(py::init<int, int, std::string, std::string, double, int,  int, std::string,  std::string, double, bool, double >())
		.def_readwrite("maturityPoint", &SwapInstrument::maturityPoint)
		.def_readwrite("numofPeriods", &SwapInstrument::numOfPeriods)
		.def_readwrite("floatTenor", &SwapInstrument::floatTenor)
		.def_readwrite("fixedTenor", &SwapInstrument::fixedTenor)
		.def_readwrite("paymType", &SwapInstrument::paymType)
		.def_readwrite("fixedRate", &SwapInstrument::fixedRate)
		.def_readwrite("spotDate", &SwapInstrument::spotDate)
		.def_readwrite("curvePointV", &SwapInstrument::curvePointV);

	py::class_<FixingInstrument>(m, "FixingInstrument")
		.def(py::init<int, double >())
		.def_readwrite("curvePoint", &FixingInstrument::curvePoint)
		.def_readwrite("curvePointV", &FixingInstrument::curvePointV);

	py::class_<InstrumentPG<SwapInstrument>>(m, "InstrumentPGLogL")
		.def(py::init<SwapInstrument&>())
		.def("pGLogLCurveFutures", &InstrumentPG<SwapInstrument>::pGCurveFutures)
		.def("pGLogLCurveSwap", &InstrumentPG<SwapInstrument>::pGCurveSwap)
		.def("pGLogLCurveCCS", &InstrumentPG<SwapInstrument>::pGCurveCCS)
		.def("pGLogLCurveCCSCubS", &InstrumentPG<SwapInstrument>::pGCurveCCSCubS)
		.def("pGLogLCurveTBS", &InstrumentPG<SwapInstrument>::pGCurveTBS);
	
	/*
	py::class_<InstrumentPG<SwapInstrument, LogLinInterp>>(m, "InstrumentPGLogL")
		.def(py::init<SwapInstrument&>())
		.def("pGLogLCurveFutures", &InstrumentPG<SwapInstrument, LogLinInterp>::pGCurveFutures)
		.def("pGLogLCurveSwap", &InstrumentPG<SwapInstrument, LogLinInterp>::pGCurveSwap)
		.def("pGLogLCurveCCS", &InstrumentPG<SwapInstrument, LogLinInterp>::pGCurveCCS)
		.def("pGLogLCurveTBS", &InstrumentPG<SwapInstrument, LogLinInterp>::pGCurveTBS);
		//.def("pGLogLZerosCurveSequentialFFS", &InstrumentPG<SwapInstrument, LogLinInterp>::pGLinZerosCurveSequentialFFS); //float-float swap
   */
	py::class_<Curve>(m, "Curve")
		.def(py::init< std::vector<FixingInstrument>& , std::vector<SwapInstrument>&>())
		.def("dRateLogL", &Curve::dRateLogL, R"pbdoc(enter and int for date and "3M" for tenor)pbdoc")
		.def("fRateLogL", &Curve::fRateLogL, R"pbdoc(enter and int for date and "3M" for tenor)pbdoc")
		.def("dRateCubS", &Curve::dRateCubS, R"pbdoc(enter and int for date and "3M" for tenor)pbdoc")
		.def("fRateCubS", &Curve::fRateCubS, R"pbdoc(enter and int for date and "3M" for tenor)pbdoc")
		.def_readwrite("points_values", &Curve::points_values)
		.def_readwrite("curveSpotDate", &Curve::curveSpotDate);
	
	py::class_<DayCountCalculator>(m, "DayCountCalculator")
		.def(py::init< std::string>())
		.def("yFrac", &DayCountCalculator::yFrac, R"pbdoc(calculates the year fraction)pbdoc");


#ifdef VERSION_INFO
	m.attr("__version__") = VERSION_INFO;
#else
	m.attr("__version__") = "dev";
#endif
}

 