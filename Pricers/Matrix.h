// Matrix.h: interface for the CMatrix class.
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

#if !defined(AFX_MATRIX_H__E6397425_F122_44A4_A3A8_53CD4D2523E7__INCLUDED_)
#define AFX_MATRIX_H__E6397425_F122_44A4_A3A8_53CD4D2523E7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <vector>

using namespace std;

template <class T>
class CMatrix  
{
public:
	int dim(int a);
	CMatrix(int a, int b);
	virtual ~CMatrix();
	bool checkBounds(int a, int b);
	void add(int a, int b);
	
	// Get/Set the value of one element.
	T operator()(const int, const int) const;
	// Get the value of one element.
	T& operator()(const int, const int);

	std::vector<T> slice(int row);
	
private:
	int d2;
	int d1;
	std::vector< std::vector<T> > m;
};

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

template<class T>
CMatrix<T>::CMatrix(int a, int b)
{
	d1 = a; d2 = b;
	m.resize(d1);  // Utilisez resize pour définir la taille du vecteur externe
	for (int i = 0; i < d1; i++)
		m[i].resize(d2, 0);  // Utilisez resize pour définir la taille des vecteurs internes avec une valeur par défaut de 0
}

template<class T>  
CMatrix<T>::~CMatrix()
{

}

template<class T>  
T CMatrix<T>::operator()(const int x, const int y) const
{
	if(checkBounds(x, y))
		return m[x][y];
}


template<class T>  
T& CMatrix<T>::operator()(const int x, const int y)
{
	if(checkBounds(x, y))
		return m[x][y];
	else
	{
		add(x, y);
		return m[x][y];
	}
}

template<class T>  
int CMatrix<T>::dim(int a)
{
	if(a == 1)
		return d1;
	else
		return d2;
}

template<class T>  
std::vector<T> CMatrix<T>::slice(int row)
{
	std::vector<T> vec;
	if(row <= d1)
		vec = m[row];

	return vec;		
}

template<class T>  
void CMatrix<T>::add(int a, int b)
{
	if(a > d1)
	{
		d1 = a;
		m.resize(d1);
	}
	if(b > d2)
	{
		d2 = b;
		for(int i = 0; i < d1; i++)
			m[i].resize(d2, 0);
	}
}

template<class T>  
inline bool CMatrix<T>::checkBounds(int a, int b)
{
	if(a < d1 || b < d2)
		return 1;
	else
		return 0;
}

#endif // !defined(AFX_MATRIX_H__E6397425_F122_44A4_A3A8_53CD4D2523E7__INCLUDED_)
