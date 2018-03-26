/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _vector_utils_H
#define _vector_utils_H

#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cassert>

//Vector scalar unary op
template<typename T, typename S> std::vector<T>& operator+=(std::vector<T>& a,  const S& s)
{				
	std::for_each(a.begin(), a.end(), [&s](T& item){ item += s; });
	return a;
}

template<typename T, typename S> std::vector<T>& operator-=(std::vector<T>& a, const S& s)
{
	std::for_each(a.begin(), a.end(), [&s](T& item){ item -= s; });	
	return a;
}

template<typename T, typename S> std::vector<T>& operator*=(std::vector<T>& a, const S& s)
{	
	std::for_each(a.begin(), a.end(), [&s](T& item){ item *= s;});
	return a;
}

template<typename T, typename S> std::vector<T>& operator/=(std::vector<T>& a, const S& s)
{
	std::for_each(a.begin(), a.end(), [&s](T& item){ item /= s; });	
	return a;
}

//Vector scalar binary op
template<typename T, typename S> std::vector<T> operator+(const std::vector<T>& a, const S& s)
{
	std::vector<T> b = a;
	return b += s;
}

template<typename T, typename S> std::vector<T> operator-(const std::vector<T>& a, const S& s)
{
	std::vector<T> b = a;
	return b -= s;
}

template<typename T, typename S> std::vector<T> operator*(const std::vector<T>& a, const S& s)
{
	std::vector<T> b=a;
	return b *= s;	
}

template<typename T, typename S> std::vector<T> operator/(const std::vector<T>& a, const S& s)
{
	std::vector<T> b = a;
	return b /= s;
}

template<typename T, typename S> std::vector<T> operator+(const S& s, const std::vector<T>& a)
{
	std::vector<T> b = a;
	return b += s;
}

template<typename T, typename S> std::vector<T> operator-(const S& s, const std::vector<T>& a)
{
	std::vector<T> b(a.size(),s);	
	return b -= a;
}

template<typename T, typename S> std::vector<T> operator*(const S& s, const std::vector<T>& a)
{
	std::vector<T> b = a;
	return b *= s;
}

template<typename T, typename S> std::vector<T> operator/(const S& s, const std::vector<T>& a)
{	
	std::vector<T> b(a.size());
	for (size_t i = 0; i < b.size(); i++) b[i] = s/a[i];
	return b;	
}

//Vector vector unary op
template<typename T> std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& b)
{
	for (size_t i = 0; i < a.size(); i++) a[i] += b[i];
	return a;
}

template<typename T> std::vector<T>& operator-=(std::vector<T>& a, const std::vector<T>& b)
{
	for (size_t i = 0; i < a.size(); i++) a[i] -= b[i];
	return a;
}

template<typename T> std::vector<T>& operator*=(std::vector<T>& a, const std::vector<T>& b)
{
	for (auto i = 0; i < b.size(); i++) a[i] *= b[i];
	return a;
}

template<typename T> std::vector<T>& operator/=(std::vector<T>& a, const std::vector<T>& b)
{
	for (size_t i = 0; i < a.size(); i++) a[i] /= b[i];
	return a;
}

//vector vector binary op
template<typename T> std::vector<T> operator*(const std::vector<T>& a, const std::vector<T>& b)
{
	std::vector<T> c = a;
	return c *= b;
}

template<typename T> std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
	std::vector<T> c = a;
	return c += b;	
}

template<typename T> std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
	std::vector<T> c = a;
	return c -= b;	
}

template<typename T> std::vector<T> operator/(const std::vector<T>& a, const std::vector<T>& b)
{
	std::vector<T> c = a;
	return c /= b;
}


//Functions
template<typename T> void log10_apply(std::vector<T>& v)
{
	std::for_each(v.begin(), v.end(), [](T& item){ item = std::log10(item); });
};

template<typename T> std::vector<T> log10(const std::vector<T>& v)
{
	std::vector<T> a = v;
	log10_apply(a); return a;
};

template<typename T> void logn_apply(std::vector<T>& v)
{
	std::for_each(v.begin(), v.end(), [](T& item){ item = std::log(item); });
};

template<typename T> std::vector<T> logn(const std::vector<T>& v)
{
	std::vector<T> a = v;
	logn_apply(a); return a;
};

template<typename T> void exp_apply(std::vector<T>& v)
{
	std::for_each(v.begin(), v.end(), [](T& item){ item = std::exp(item); });
};

template<typename T> std::vector<T> exp(const std::vector<T>& v)
{
	std::vector<T> a = v;
	exp_apply(a); return a;
};

template<typename T> T min(const std::vector<T>& v)
{
	return *std::min_element(v.cbegin(), v.cend());
};

template<typename T> T max(const std::vector<T>& v)
{
	return *std::max_element(v.cbegin(), v.cend());
};

template<typename T> T sum(const std::vector<T>& v)
{
	return std::accumulate(v.cbegin(), v.cend(), 0.0);
};

template<typename T> T dot(const std::vector<T>& u, const std::vector<T>& v)
{
	return std::inner_product(u.cbegin(), u.cend(), v.cbegin(), 0.0);
};

template<typename T> std::vector<T> sqrt(const std::vector<T>& v)
{
	std::vector<T> o(v.size());
	for (size_t i=0; i<v.size(); ++i) o[i]=sqrt(v[i]);
	//std::transform(v.cbegin(), v.cend(), o.cbegin(), (double(*)(double)) sqrt);
	return o;
};

template<typename T> T mean(const std::vector<T>& v)
{
	return sum(v) / v.size();
};

template<typename T> T variance(const std::vector<T>& v)
{
	T vmean = mean(v);
	T init = 0.0;
	T sum = std::accumulate(v.cbegin(), v.cend(), init, [&vmean](T& total, const T& item){ return total += std::pow(item - vmean, 2.0); });
	return sum / v.size();
};

template<typename T> T stddev(const std::vector<T>& v)
{
	return sqrt(variance(v));
};

template<typename T> T fabs(const std::vector<T>& v)
{
	std::vector<T> a = v;
	std::for_each(a.begin(), a.end(), [](T& item){ item = std::fabs(item); });
	return a;
};


template<typename T>
void append(std::vector<T>& a, const std::vector<T>& b){	
	//This does not work for Intel: a.insert(std::end(a), std::begin(b), std::end(b));	
	a.insert(a.end(), b.begin(), b.end());
}

template<typename T>
void prepend(std::vector<T>& a, const std::vector<T>& b){
	//This does not work for Intel: a.insert(std::begin(a), std::begin(b), std::end(b));
	a.insert(a.begin(), b.begin(), b.end());
}

template<typename T>
std::vector<T> concaternate(const std::vector<T>& a, const std::vector<T>& b){
	std::vector<T> c=a;
	//This does not work for Intel: c.insert(std::end(c), std::begin(b), std::end(b));
	c.insert(c.end(), b.begin(), b.end());
	return c;
}

//Resize 2d array
template<typename T> void resize(std::vector<std::vector<T>>& m, size_t nrows, size_t ncols)
{
	m.resize(nrows);
	for (size_t i = 0; i < nrows; i++){
		m[i].resize(ncols);
	}
};

template<typename TA, typename TB> void cast(const std::vector< std::vector<TA> >& a, std::vector< std::vector<TB> >& b)
{
	size_t nr = a.size();
	b.resize(nr);
	for (size_t i = 0; i<nr; i++){
		size_t nc = a[i].size();
		b[i].resize(nc);
		for (size_t j = 0; j<nc; j++){
			b[i][j] = a[i][j];
		}
	}
}

//Stats on raw pointer
template<typename T> T min(const size_t n, const T* v)
{
	T m = *std::min_element(v, v + n);
	return m;
};

template<typename T> T max(const size_t n, const T* v)
{
	T m = *std::max_element(v, v + n);
	return m;
};

template<typename T> T sum(const size_t n, const T* v)
{
	return std::accumulate(v, v + n, 0.0);
};

template<typename T> T mean(const size_t n, const T* v)
{
	return sum(n, v) / n;
};

template<typename T> T stddev(const size_t n, const T* v)
{
	T vmean = mean(n, v);
	T init = 0.0;
	T sum = std::accumulate(v, v+n, init, [&vmean](T& total, const T& item){return total += std::pow(item - vmean, 2.0);});
	return sum / n;
};

template<typename ForwardIterator, typename T>
void arange(ForwardIterator first, ForwardIterator last, T begin, T end)
{
	size_t d = std::distance(first,last);
	if (d==0) return;
	T inc = (end - begin) / (d-1);
	size_t i = 0;
	while (first != last) *first++ = begin + inc * i++;
}

template<typename T>
std::vector<T> increment(const size_t n, const T start = 0, const T inc = 1)
{
	std::vector<T> v(n);
	v[0] = start;
	for (size_t i = 1; i<n; i++) v[i] = v[i - 1] + inc;
	return v;
};

template<typename T>
std::vector<T> cumsum(std::vector<T> a)
{
	std::vector<T> cs(a.size());
	std::partial_sum(a.begin(), a.end(), cs.begin(), std::plus<T>());
	return cs;
}

template<typename T>
std::vector<std::vector<T>> matinv(std::vector<std::vector<T>>& mat){
	size_t n = mat.size(); 
	size_t i, j;
	std::vector<std::vector<T>> inv(n, std::vector<T>(n, 0));
	T determinant = 0;
	for(i = 0; i < n; i++)
		determinant = determinant + (mat[0][i] * (mat[1][(i+1)%3] * mat[2][(i+2)%3] - mat[1][(i+2)%3] * mat[2][(i+1)%3]));
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			inv[i][j] = ((mat[(j+1)%n][(i+1)%n] * mat[(j+2)%n][(i+2)%n]) - (mat[(j+1)%n][(i+2)%n] * mat[(j+2)%n][(i+1)%n]))/ determinant;
	return inv;
}
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename T>
T logdet(const std::vector<std::vector<T>>& a_const, T& sign) {
	size_t n = a_const.size();
	/*
	std::vector<std::vector<T>> l(n,std::vector<T>(n,0));
	std::vector<std::vector<T>> u(n,std::vector<T>(n,0));
	size_t i,k,j,p;
	*/
	//********** LU decomposition *****//
	/*
	for(k=0;k<n;k++)
	{
		u[k][k]=1;
		for(i=k;i<n;i++)
		{
			sum=0;
			for(p=0;p<k-1;p++)
				sum+=l[i][p]*u[p][k];
			l[i][k]=a[i][k]-sum;
		}

		for(j=k+1;j<n;j++)
		{
			sum=0;
			for(p=0;p<k-1;p++)
				sum+=l[k][p]*u[p][j];
			u[k][j]=(a[k][j]-sum)/l[k][k];
		}
	}
	*//*
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (j < i)
				l[j][i] = 0;
			else
			{
				l[j][i] = a[j][i];
				for (k = 0; k < i; k++)
				{
					l[j][i] = l[j][i] - l[j][k] * u[k][i];
				}
			}
		}
		for (j = 0; j < n; j++)
		{
			if (j < i)
				u[i][j] = 0;
			else if (j == i)
				u[i][j] = 1;
			else
			{
				u[i][j] = a[i][j] / l[i][i];
				for (k = 0; k < i; k++)
				{
					u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
				}
			}
		}
	}
	*/
	std::vector<std::vector<T>> a = a_const; // I hope this copy is deep.
	const T TINY=1.0e-20;
	int i,imax,j,k;
	T big,dum,sum,temp;

	std::vector<T> vv(n,0);
	std::vector<int> indx(n,0); // for logdet we cat remove indx and d
	T d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) return std::numeric_limits<double>::infinity();//"Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ((dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			d = -d;
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}

	// sum logs of diagonals
	T ld = 0;
	sign = 1;
	//for(k=0;k<n;k++) ld += std::log(u[k][k]) + std::log(l[k][k]);
	for(k=0;k<n;k++) {
		ld += std::log(std::fabs(a[k][k])); // log of 1 is zero
		sign *= sgn(a[k][k]);
	}
	return ld;
}

#endif

