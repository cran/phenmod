#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include "pimModel.h"

int * pimModel(double *params, int *lcdoy, int *year, double *latitude, double *temperaturevec, int *modelno, int *budburstdoy){
	int doyleafcolouring, yearbudburst, modelnumber, i=0, validvalue=0, yearchanged=0, doy=0;
	int previousyeardaynumber = 365;
	double param[10];
	double temperatures[365];
	double lat, inhibitorvalue=1, promotervalue=0, daylength=0, inhibitorprevious=1, promoterprevious=0;
	
	//copy values from pointer variables
	doyleafcolouring = lcdoy[0];
	yearbudburst = year[0];
	lat = latitude[0];
	modelnumber = modelno[0];

	for (i=0; i < 10; i++){
		param[i] = params[i];
	}

	for (i=0; i < 365; i++){
		temperatures[i] = temperaturevec[i];
	}

	doy = doyleafcolouring;
	previousyeardaynumber+=isLeapYear(yearbudburst-1);
	
	for (i=0; i<365; i++){
		if (yearchanged == 0){
			daylength = getDaylength(lat, doy, yearbudburst-1);
		} else {
			daylength = getDaylength(lat, doy, yearbudburst);
		}

		inhibitorprevious = inhibitorvalue;
		promoterprevious = promotervalue;

		inhibitorvalue = inhibitorprevious + deltaInhibitor(param[0], param[1], param[4], param[5], 
			 					param[6], inhibitorprevious, 
								temperatures[i], daylength, modelnumber);
		promotervalue = promoterprevious + deltaPromoter(param[2],param[3], param[7], param[8], 
								param[9], inhibitorprevious,
								promoterprevious, temperatures[i], 
								daylength, modelnumber);	

		if (promotervalue >= 1){
			validvalue = 1;
			break;
		} else {
			doy++;
			if (doy > previousyeardaynumber){
				doy = 1;
				yearchanged = 1;
			}
		}	
	}	

	if (validvalue == 0){
		*budburstdoy = 0;
	} else {
		*budburstdoy = doy;
	}

	return budburstdoy;
}

double pimTemperature(double Tmin, double Topt, double Tmax, double T){
	double temperature=0;
	if ((Tmin <= T)&&(T <= Topt)){
		temperature = (T - Tmin) / (Topt - Tmin);
	}
	if ((Topt < T) && (T <= Tmax)){
		temperature = (Tmax - T) / (Tmax - Topt);
	}
	if ((T < Tmin)||(T > Tmax)){
		temperature = 0;
	}
	return temperature;
}	

double deltaInhibitor(double a1, double a2, double Tmini, double Topti, 
			double Tmaxi, double inhprevious, double T,
			double daylength, int modelnumber){
	double dI=0;
	double nightlength=0;

	switch(modelnumber){
		case 1:
		case 2:
		case 4:
		case 5: 
			dI = -a2 * pimTemperature(Tmini, Topti, Tmaxi, T) * inhprevious; 
			break;
		case 3:
		case 6:
			nightlength = (24.0 - daylength);
			dI = (a1 * (nightlength / 24.0)) - (a2 * pimTemperature(Tmini, Topti, Tmaxi, T) * inhprevious);
			break;
		case 7:
		case 8:
		case 10:
		case 11:
			dI = -a2 * pimTemperature(Tmini, Topti, Tmaxi, T) * inhprevious * (daylength / 24.0);
			break;
		case 9:
		case 12:
			nightlength = (24 - daylength);
			dI = (a1 * (nightlength / 24.0)) - (a2 * pimTemperature(Tmini, Topti, Tmaxi, T) * inhprevious * (daylength / 24.0));
			break;
		default:
			break;
	}
	return dI;
}

double deltaPromoter(double a3, double a4, double Tminp, double Toptp, 
			double Tmaxp, double inhprevious, double proprevious, double T,
			double daylength, int modelnumber){
	double dP=0, nightlength=0;

	switch(modelnumber){
		case 3:
		case 9:
			dP = a3 * pimTemperature(Tminp, Toptp, Tmaxp, T) * (1.0 - inhprevious);
			break;
		case 1:
		case 7:
			dP = a3 * pimTemperature(Tminp, Toptp, Tmaxp, T) * (1.0 - inhprevious) - a4 * proprevious;
			break;
		case 4:
		case 10:
			dP = a3 * pimTemperature(Tminp, Toptp, Tmaxp, T) * (1.0 - inhprevious) * (daylength / 24.0) - a4 * proprevious;
			break;
		case 2:
		case 8:
			nightlength = (24.0 - daylength);
			dP = a3 * pimTemperature(Tminp, Toptp, Tmaxp, T) * (1.0 - inhprevious) - a4 * proprevious * (nightlength / 24.0);
			break;
		case 5:
		case 11:
			nightlength = (24.0 - daylength);
			dP = a3 * pimTemperature(Tminp, Toptp, Tmaxp, T) * (1.0 - inhprevious) * (daylength / 24.0) - a4 * proprevious * (nightlength / 24.0);
			break;
		case 6:
		case 12:
			nightlength = (24.0 - daylength);
			dP = a3 * pimTemperature(Tminp, Toptp, Tmaxp, T) * (1.0 - inhprevious) * (daylength / 24.0);
			break;
		default: 
			break;
	}
	return dP;
}

//returns daylength of latitude lat at doy in year
double getDaylength(double lat, int doy, int year){
	double daylength=0;
	daylength = 2.0 * radianToDegree(acos(-tan(degreeToRadian(lat)) * tan(degreeToRadian(getDeclination(year,doy)))))/15.0;
	return daylength;
}

//returns declination of earth at doy in year
double getDeclination(int year, int doy){
	double decl=999;
	double declsol=23.43889;
	int yearlength=0;
	
	if (isLeapYear(year)==0){
		yearlength=365;	
	} else {
		yearlength=366;
	}	

	decl = sin(((doy-81)/yearlength)*2*M_PI)*declsol;

	return decl;
}

//is year YYYY a leap year ?
int isLeapYear(int year) {
	if (year % 4 == 0)
		if (!(year % 100 == 0))
			return 1;
		else
			if (year % 400 == 0)
				return 1;
	return 0;
}

double radianToDegree(double rad){
	double degree=0;
	degree = rad * (180.0 / M_PI);
	return degree;
}

double degreeToRadian(double degree){
	double radian=0;
	radian = degree * (M_PI / 180.0);
	return radian;
}

//Temperature-Sum-Model
int * tsModel(double *params, int *lcdoy, int *year, double *temperaturevec, int *budburstdoy){
	int doyleafcolouring, yearbudburst, i=0, validvalue=0, doy=0, previousyeardaynumber = 365; 
	double temperatures[365], param[2];
	double rate = 0, Tb, Fstar, t0;
	
	//copy values from pointer variables
	doyleafcolouring = lcdoy[0];
	yearbudburst = year[0];
	Tb = params[0];
	Fstar = params[1];

	for (i=0; i < 365; i++){
		temperatures[i] = temperaturevec[i];
	}

	//set starting day to first day of the year
	//and search first value in temperatur-vector
	previousyeardaynumber+=isLeapYear(yearbudburst-1);
	t0 = previousyeardaynumber - doyleafcolouring + 1;
	doy = 1;

	//iterate over temperature values
	for (i=t0-1; i<365; i++){
		if (temperatures[i] > Tb){
			rate += temperatures[i] - Tb;
		}
		if (rate >= Fstar){
			validvalue = 1;
			break;
		} else {
			doy++;
		}
	}

	//check if budburst-doy was found
	if (validvalue == 0){
		*budburstdoy = 0;
	} else {
		*budburstdoy = doy;
	}

	return budburstdoy;
}
