#pragma once
#include <string>
#include <cassert>
namespace py = pybind11;


bool isLeapYear(int yr)
{
	if ((yr % 4 == 0) && (yr % 100 != 0))
	{
		return true;
	}
	else if ((yr % 100 == 0) && (yr % 400 == 0))
	{
		return true;
	}
	else if (yr % 400 == 0)
	{
		return true;
	}
	else
	{
		return false;
	}

}

//43831 is the excel date for January 1, 2020, and it is a Wednesday
//43836 is the excel date for January 6, 2020, and it is a Monday
//I'll assume Python datetime convention: 
//date.weekday()
//Return the day of the week as an integer, where Monday is 0 and Sunday is 6. For example, date(2002, 12, 4).weekday() == 2, a Wednesday.See also isoweekday().
int daysOfWeek(int date)
{
	int n = (date - 43836) % 7;
	return n;
}
//months: january =0,...., december=11
int daysInAMonth(int month, int yr)
{
	//0,2,4,6,7,9,11
	//1,3,5,8,10,
	month = month % 12;
	switch (month) {
	case 0:
		return 31;
	case 2:
		return 31;
	case 4:
		return 31;
	case 6:
		return 31;
	case 7:
		return 31;
	case 9:
		return 31;
	case 11:
		return 31;
	case 1:
		if (isLeapYear(yr))
			return 29;
		else
			return 28;
	case 3:
		return 30;
	case 5:
		return 30;
	case 8:
		return 30;
	case 10:
		return 30;
	}
}



int findTheYear(int date)  
{
	int tryThisMany = (date - 43831) / 365+1;
	int startYear = 2020;
	int startDate = 43830;
	int endDate = startDate;
	for (int i = 0; i < tryThisMany; ++i)
	{
		
		if (isLeapYear(startYear))
		{
			endDate = endDate + 366;//last day of the year
		}
		else
		{
			endDate = endDate + 365; //last day of the year
		}
		if (date <= endDate && date > startDate)
		{
			return startYear;
		}
		else
		{
			startYear += 1;
			startDate = endDate;
		}

	}
	 

}

int findTheMonth(int date) //this could go year by year and it would be faster
{
	int yr = findTheYear(date);
	int tryThisMany = (date - 43831) / 28 +1;
	int currentMonthEndDate = 43830;
	for (int i = 0; i < tryThisMany; ++i) //go month by month
	{
		int whichMonth = i % 12;
		//if (i > 0 && whichMonth == 0)
		//	yr += 1;
		int nextMonthEndDate = currentMonthEndDate + daysInAMonth(whichMonth, 2000 + i / 12);
		if (date > currentMonthEndDate && date <= nextMonthEndDate)
			return whichMonth;
		else
		{
			currentMonthEndDate = nextMonthEndDate;
		}
		 
	}
}

int findTheDayOfTheMonth(int date) //this could go year by year and it would be faster
{
	int yr = findTheYear(date);
	int tryThisMany = (date - 43831) / 28 +1;
	int currentMonthEndDate = 43830; // December 31, 2019
	for (int i = 0; i < tryThisMany; ++i) //go month by month
	{
		int whichMonth = i % 12;
		int nextMonthEndDate = currentMonthEndDate + daysInAMonth(whichMonth, 2000 + i/12);
		if (date > currentMonthEndDate && date <= nextMonthEndDate)
			return date - currentMonthEndDate;
		else
		{
			currentMonthEndDate = nextMonthEndDate;
		}
	}
}



int findTheFirstDayMonthYear(int month, int year)
{
	//found the date for the first day of a given year month
	//follow the excel date convention of this code
	//43831 is Jan 1, 2020.
	int lastdayoftheprevyear = 43830;
	int numofyears = year - 2020;
	int numofmonths = month - 0;
	int aux_year = 0;
	//find the last day of the year just before the argument year
	for (int i = 0; i < numofyears; ++i)
	{
		aux_year = 2020 + i;
		lastdayoftheprevyear += (isLeapYear(aux_year) ? 366 : 365);
	}
	//add this last day the days of the months
	//note that lastdayoftheprevyear becomes the last day of the previous month
	//don't get confused!
	for (int i = 0; i < numofmonths; ++i)
	{
		int days = daysInAMonth(i, year);
		lastdayoftheprevyear += days;
	}
	return lastdayoftheprevyear + 1;

}

int findIMMDateGivenMonthYear(int month, int year)
{
	//get the third wednesday of the found month-year
	int firstday = findTheFirstDayMonthYear(month, year);
	while (daysOfWeek(firstday) != 2)
	{
		firstday += 1;
	}
	return firstday + 14;

}

int shiftByThreeMonths(int date, bool immShift = false)
{

	int whichDay = findTheDayOfTheMonth(date);
	int whichMonth = findTheMonth(date);
	int whichYear = findTheYear(date);
	//py::print("day, month, year: ", whichDay, whichMonth, whichYear);
	int lastDayOfThePrevMonth = date - whichDay;
	if (!immShift)
	{
		for (int i = 0; i < 3; ++i)
		{
			if ((whichMonth + i) > 11)
			{
				whichYear += 1;
				whichMonth -= 12;
			}
			//py::print("lastDayOfThePrevMonth: ", daysInAMonth((whichMonth + i) % 12, whichYear));
			//py::print("lastDayOfThePrevMonth: ", lastDayOfThePrevMonth);
			lastDayOfThePrevMonth += daysInAMonth((whichMonth + i) % 12, whichYear);

		}
		return lastDayOfThePrevMonth + whichDay;
	}
	else
	{
		if ((whichMonth + 3) % 12 < whichMonth) whichYear += 1;
		//py::print("not exact Libor tenor!");
		return findIMMDateGivenMonthYear((whichMonth + 3) % 12, whichYear);
	}
	//py::print("lastDayOfThePrevMonth, whichDay", lastDayOfThePrevMonth, whichDay);

}

int shiftBySixMonths(int date)
{
	int shiftedDate = date;
	for (int i = 0; i < 2; ++i)
	{
		shiftedDate = shiftByThreeMonths(shiftedDate, false);
	}
	return shiftedDate;
}

int shiftByOneYear(int date)
{
	int shiftedDate = date;
	for (int i = 0; i < 4; ++i)
	{
		shiftedDate = shiftByThreeMonths(shiftedDate);
	}
	return shiftedDate;
}

//might need modified following here as most of the rest:
//e.g.: mo = Jan, yr = 2023, day = 30, then you shift to the march when you shift by one month!
int shiftByOneMonth(int date)
{

	int whichDay = findTheDayOfTheMonth(date);
	int whichMonth = findTheMonth(date);
	int whichYear = findTheYear(date);

	return  findTheFirstDayMonthYear(whichMonth, whichYear) + daysInAMonth(whichMonth, whichYear) - 1 + whichDay;
}

struct DayCountCalculator
{
	std::string dayCMethod;
	//int dStart;
	//int dEnd;

	//DayCountCalculator(std::string dayCMethod, int dStart, int dEnd) : dayCMethod(dayCMethod), dStart(dStart), dEnd(dEnd) {}
	DayCountCalculator(std::string dayCMethod) : dayCMethod(dayCMethod)  {}

	//all dates inclusive calculation
	double yFrac(int dStart, int dEnd)
	{
		assert(dStart > 43830);
		assert(dEnd >= dStart);
		if (dayCMethod == "30/360")
		{
			//follows ISDA 2006, Section 4.16
			int y1 = findTheYear(dStart);
			int y2 = findTheYear(dEnd+1);
			int m1 = findTheMonth(dStart);
			int m2 = findTheMonth(dEnd + 1);
			int d1 = findTheDayOfTheMonth(dStart);
			int d2 = findTheDayOfTheMonth(dEnd + 1);
			if (d1 > 30) d1 = 30;
			if (d2 == 31 && d1 > 29) d2 = 30;
			//py::print("year, month, date", y1, y2, m1, m2, d1, d2);
			return (360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1)) / 360.0;
		}

		if (dayCMethod == "Act/360")
		{
			return (dEnd - dStart) / 360.0;
		}

		if (dayCMethod == "Act/365")
		{
			return (dEnd - dStart) / 365.0;
		}

		//this is for the fixed leg for USD 6M...
		if (dayCMethod == "wrongHalf")
		{
			return 0.5;
		}
		if (dayCMethod == "wronghQuarter")
		{
			return 0.25;
		}
		else return 0;
			
	}

};


