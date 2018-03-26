/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Authors: Ross C. Brodie, Geoscience Australia.
         Laurence Davies, Geoscience Australia.
*/

#ifndef _RANDOM_H_
#define _RANDOM_H_
#include <vector>
#include <algorithm>
#include <iterator>
#include "vector_utils.h"
#include "matrix_ops.h"
using namespace std;

void seedrand();
void seedrand(unsigned int seed);
double urand();
double urand(double rmin, double rmax);
size_t irand(size_t imin, size_t imax);
int    irand(int imin, int imax);
double nrand();
std::vector<double> nrand(size_t n);
void nrand(size_t n, double* x);
void nrand(size_t n, double* x, double mean, double std);

std::vector<double> mvnrand_lowercholesky(const MatrixDouble& L);
std::vector<double> mvnrand_covariance(const MatrixDouble& C);

double gaussian_pdf(double mean, double std, double x);
double mvgaussian_pdf(const VectorDouble& m0, const MatrixDouble& C, const VectorDouble& m);

template<typename T>
std::vector<size_t> randsamplesortedsequentialindex(size_t k, std::vector<T> w){
	// return a vector of size k sampling ints from range 0::n-1 using multinomial weights w
	std::vector<size_t> ind(k);
	size_t n = w.size();
	std::vector<T> w_cumsum = cumsum(w);
	T upper = w_cumsum[w_cumsum.size()-1];
	T v = urand(0,upper/k); // FIXME use a non-typecast version of urand
	//std::cout << "Upper of weight range = " << upper << std::endl;
	for (size_t i=0;i<k;++i){
		T stepi = v + upper * (T)i/(T)k;
		ind[i] = std::distance(w_cumsum.begin(),std::lower_bound(w_cumsum.begin(),w_cumsum.end(),stepi));
		ind[i] = std::min(n-1,ind[i]); // in case the above yields n which results in a segfault.
		//std::cout << "v=" << stepi << ", ind["<<i<<"]="<<ind[i]<<"; ";
	}
	std::sort(ind.begin(),ind.end()); // Not needed.
	return ind;
}

template<typename T>
std::vector<size_t> randsamplesortedindex(size_t n, size_t k, std::vector<T> w){
	// return a vector of size k sampling ints from range 0::n-1 using multinomial weights w
	std::vector<size_t> ind(k);
	std::vector<T> w_cumsum = cumsum(w);
	T upper = w_cumsum[w_cumsum.size()-1];
	//std::cout << "Upper of weight range = " << upper << std::endl;
	//for (size_t i=0;i<w_cumsum.size();++i){
	for (size_t i=0;i<k;++i){
		T v = urand(0,upper); // FIXME use a non-typecast version of urand
		ind[i] = std::distance(w_cumsum.begin(),std::lower_bound(w_cumsum.begin(),w_cumsum.end(),v));
		//std::cout << "v=" << v << ", ind["<<i<<"]="<<ind[i]<<"; ";
	}
	std::sort(ind.begin(),ind.end());
	return ind;
}

template<typename T>
std::vector<T> randsamplesorted(size_t n, size_t k, std::vector<T> w){
	// return a vector of size k sampling ints from range 0::n-1 using multinomial weights w
	std::vector<T> samp(k);
	std::vector<T> w_cumsum = cumsum(w);
	T upper = w_cumsum[w_cumsum.size()-1];
	//std::cout << "Upper of weight range = " << upper << std::endl;
	for (size_t i=0;i<w_cumsum.size();++i){
		T v = urand(0,upper); // FIXME use a non-typecast version of urand
		auto lower = std::lower_bound(w_cumsum.begin(),w_cumsum.end(),v);
		samp[i] = *lower;
		//std::cout << "v=" << v << ", samp["<<i<<"]="<<samp[i]<<"; ";
	}
	std::sort(samp.begin(),samp.end());
	return samp;
}

#endif
