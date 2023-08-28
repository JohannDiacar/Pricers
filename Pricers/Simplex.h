// Simplex.h: interface for the CSimplex class.
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

#if !defined(AFX_SIMPLEX_H__CC41AD4A_3850_49EA_93E8_AAC6139BAE1A__INCLUDED_)
#define AFX_SIMPLEX_H__CC41AD4A_3850_49EA_93E8_AAC6139BAE1A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <iostream>

#include <vector>
#include <math.h>
#include "Matrix.h"
#include <functional>
class CSimplex  
{
public:
	void reset();
	// Dans la déclaration de CSimplex
	CSimplex(int ndim, double tol, int maxIter, std::function<double(std::vector<double>)> func);
	virtual ~CSimplex();
	double amoeba(std::vector<double>& vPoint);
	int calls();
	void write();

private:
	bool m_Write;
	double m_EPS;
	double m_Tol;
	int m_NoDim;
	int m_NoEval;
	int m_MaxIter;
	std::function<double(std::vector<double>)> m_Func;
	double EPS();
	double feval(std::vector<double> vPoints);
	bool check_tol(double fmax, double fmin);
	void contract(CMatrix<double>& mSimplex, std::vector<int> vEx, std::vector<double>& vFx);
	int update(CMatrix<double>& mSimplex, const std::vector<double> vMid, 
		const std::vector<double> vLine, std::vector<int> vEx, double scale, double& fmax);
	std::vector<double> bearings(CMatrix<double> mSimplex, std::vector<int> vEx, std::vector<double>& vMid);
	void extremes(const std::vector<double> vFx, std::vector<int>& vEx);
	std::vector<double> evaluate(CMatrix<double> mSimplex);
	CMatrix<double> make_simplex(std::vector<double> vPoint);
};

#endif // !defined(AFX_SIMPLEX_H__CC41AD4A_3850_49EA_93E8_AAC6139BAE1A__INCLUDED_)
