// Simplex.cpp: implementation of the CSimplex class.
//
//////////////////////////////////////////////////////////////////////

/*
    Copyright (C) 2008  Colin Caprani - www.colincaprani.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "pch.h"
#include "Simplex.h"
#include <iomanip>
using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CSimplex::CSimplex(int ndim, double tol, int maxIter, std::function<double(std::vector<double>)> func)
{
	m_NoDim = ndim;
	m_Tol = tol;
	m_Func = func;
	m_MaxIter = maxIter;
	m_NoEval = 0;
	m_EPS = EPS();
	m_Write = true;
}

CSimplex::~CSimplex()
{

}

double CSimplex::amoeba(std::vector<double>& vPoint)
{
	// This is the main function that performs the minimization
	// it accepts initial parameter values and returns those at the minimum
	// in the same vector. Alos, the value of the minimum is returned

	std::vector<int> vEx(3,0);
	std::vector<double> vFx(m_NoDim+1,0.0);
	std::vector<double> vMid(m_NoDim,0.0);
	std::vector<double> vLine(m_NoDim,0.0);
	
	CMatrix<double> mSimplex = make_simplex(vPoint);
	vFx = evaluate(mSimplex);

	const int WIDTH = 11;

	if(m_Write)
		cout << setw(4) << "Iter.: " << "  " << setw(WIDTH) << "Y:" << "   " << setw(WIDTH) << "X[]: "  << endl;

	int count = 0;

	while(m_NoEval < m_MaxIter)
	{
		count++;
		extremes(vFx, vEx);
		vLine = bearings(mSimplex,vEx,vMid);

		if( check_tol(vFx[ vEx[0] ], vFx[ vEx[2] ]) )
			break;

		double& fhigh = vFx[vEx[0]];
		double& fnexthigh = vFx[vEx[1]];
		double& flow = vFx[vEx[2]];

		update(mSimplex,vMid,vLine,vEx,-1.0,fhigh);

		if(fhigh < flow)
			update(mSimplex,vMid,vLine,vEx,-2.0,fhigh);
		else if(fhigh >= fnexthigh)
		{
			if( !update(mSimplex,vMid,vLine,vEx,0.5,fhigh) )
				contract(mSimplex,vEx,vFx);
		}

		if(m_Write)
		{
			cout << std::setw(4) << count << "   " << setw(WIDTH) << vFx[vEx[2]] << "   ";
			for(int i = 0; i<m_NoDim; i++)
				cout << setw(WIDTH) << mSimplex(i,i) << "   ";
			cout << endl;
		}
	}
	
	vPoint = mSimplex.slice(vEx[2]);
	double fmin = vFx[vEx[2]];

	return fmin;
}

bool CSimplex::check_tol(const double fmax, const double fmin)
{
	double delta = fabs(fmax - fmin);
	double accuracy = (fabs(fmax)+fabs(fmin))*m_Tol;

	bool converge = (delta < (accuracy + m_EPS)) ? 1 : 0;

	return converge;
}


double CSimplex::feval(const std::vector<double> vPoints)
{
	m_NoEval++;
	double ans = (m_Func)(vPoints);
	return ans;
}

CMatrix<double> CSimplex::make_simplex(const std::vector<double> vPoint)
{
	// this function creates the initial simplex, that is, the ndim+1
	// points that are to be evauluated. The starting points are perturbed
	// slightly to begint he process

	CMatrix<double> mSimplex(m_NoDim+1,m_NoDim);

	for(int i = 0; i < m_NoDim+1; i++)
	{
		for(int j = 0; j < m_NoDim; j++)
			mSimplex(i,j) = vPoint[j];
	}

	for(int i = 0; i < m_NoDim; i++)
		mSimplex(i,i) *= 1.01;		// slightly disturb the parameters

	return mSimplex;
}

std::vector<double> CSimplex::evaluate(CMatrix<double> mSimplex)
{
	// This function evaluates the function at the simplex points

	std::vector<double> vFx(m_NoDim+1,0.0);

	for(int i = 0; i < m_NoDim+1; i++)
		vFx[i] = feval(mSimplex.slice(i));

	return vFx;
}

void CSimplex::extremes(const std::vector<double> vFx, std::vector<int>& vEx)
{
	// This function finds the extremes of the function at the simplex.
	// Finds the highest (worst), next highest (sencond worst) and the 
	// lowest (best) points

	int iHigh = 0;	int iLow = 0;	int iNextHigh = 0;

	if(vFx[0] > vFx[1])
	{ iHigh = 0; iNextHigh = 1; iLow = 1; }
	else
	{ iHigh = 1; iNextHigh = 0; iLow = 0; }

	for(int i = 2; i < m_NoDim + 1; i++)
	{
		double fi = vFx[i];

		if(fi <= vFx[iLow])
		{	iLow = i; }
		else if(fi > vFx[iHigh])
		{	iNextHigh = iHigh;	iHigh = i; }
		else if(fi > vFx[iNextHigh])
		{	iNextHigh = i; }
	}

	vEx[0] = iHigh;	vEx[1] = iNextHigh;	vEx[2] = iLow;
}

std::vector<double> CSimplex::bearings(CMatrix<double> mSimplex, const std::vector<int> vEx, std::vector<double>& vMid)
{
	// This function finds the direction for optimization by drawing a line through
	// the worst point and the centroid of all the other points.

	int iWorst = vEx[0];

	std::vector<double> vLine(m_NoDim,0.0);
	vMid.assign(m_NoDim,0.0);

	for(int i = 0; i < m_NoDim + 1; i++)
	{
		if(i != iWorst)
		{
			for(int j = 0; j < m_NoDim; j++)
				vMid[j] += mSimplex(i,j);
		}
	}

	for(int j = 0; j < m_NoDim; j++)
	{
		vMid[j] /= m_NoDim;
		vLine[j] = mSimplex(iWorst,j) - vMid[j];
	}

	return vLine;
}

int CSimplex::update(CMatrix<double>& mSimplex, const std::vector<double> vMid, const std::vector<double> vLine, 
					 const std::vector<int> vEx, double scale, double& fmax)
{
	// This function updates the simplex if the bearings resulted in a new minimum
	// through reflection or expansion as sepcified by "scale"

	int update = 0;
	std::vector<double> vNext(m_NoDim,0.0);

	for(int i = 0; i < m_NoDim; i++)
		vNext[i] = vMid[i] + scale*vLine[i];

	double fx = feval(vNext);
	
	if(fx < fmax)
	{
		for(int j = 0; j < m_NoDim; j++)
			mSimplex(vEx[0],j) = vNext[j];
		fmax = fx;
		update = 1;
	}

	return update;
}

void CSimplex::contract(CMatrix<double>& mSimplex, const std::vector<int> vEx, std::vector<double>& vFx)
{
	// This function contracts the existing simplex in the hope of achieving an update

	int iLow = vEx[2];

	for(int i = 0; i < m_NoDim+1; i++)
	{
		if(i != iLow)
		{
			for(int j = 0; j < m_NoDim; j++)
				mSimplex(i,j) = (mSimplex(iLow,j)+mSimplex(i,j))*0.5;
			vFx[i] = feval(mSimplex.slice(i));
		}
	}
}

double CSimplex::EPS()
{
	// This function calculates the machine epsilon for maximal accuracy of the algorithm

	double a; double b;
	int POW = 0;
	int iEps;

	do{
		POW++;
		a = 1.0 + pow(2,-1*POW);
		b = a - 1.0;
	}while(b > 0);

	iEps= POW - 1;
	b = pow(2,-1*iEps);

	return b;
}

int CSimplex::calls()
{
	return m_NoEval;
}

void CSimplex::write()
{
	m_Write = true;
}


void CSimplex::reset()
{
	m_Write = false;
	m_NoEval = 0;
}
