#include <string>

using namespace std;
int tenorsToDate(const string& t)
{
	int dateVal = 0;
	switch (t)
	{
	case '1Year': dateVal = 1 * 364; break;
	case '5Year': dateVal = 5 * 364; break;
	case '30Year': dateVal = 30 * 364; break;
	default:;
	}

	return dateVal;
}