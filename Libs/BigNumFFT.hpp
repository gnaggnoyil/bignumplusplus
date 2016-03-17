#ifndef _BIG_NUM_HPP_
#error "This header must be included through BigNum.hpp"
#endif // _BIG_NUM_HPP_

#ifndef _BIG_NUM_FFT_HPP_
#define _BIG_NUM_FFT_HPP_

#include <cstddef>
#include <cmath>
#include <utility>

namespace bignum{
	
	namespace _utility{
		
		template <typename Ring, typename InFunc, typename OutFunc>
		void fft1DPower2(std::size_t sizeN, Ring omega, InFunc &&getIn, OutFunc &&getOut){
			assert((sizeN & ((~sizeN) + 1)) == sizeN);
			
			using std::pow;
			using std::swap;
			
			using Ele = typename std::remove_reference<typename std::result_of<OutFunc(std::size_t)>::type>::type;
			
			for(std::size_t i = 0;i < sizeN;++i){
				getOut(i) = getIn(i);
			}
			
			std::size_t preRev = 0, rev;
			for(std::size_t i = 1;i < sizeN;++i){
				std::size_t tmp = sizeN >> 1;
				for(rev = preRev;(rev & tmp) != 0;tmp >>= 1){
					rev -= tmp;
				}
				rev += tmp;
				preRev = rev;
				
				if(rev < i){
					swap(getOut(i), getOut(rev));
				}
			}
			
			for(std::size_t m = 2;m <= sizeN;m <<= 1){
				Ring w0 = pow(omega, sizeN / m);
				Ring w(1);
				for(std::size_t i = 0, mh = m >> 1;i < mh;++i){
					for(std::size_t j = i;j < sizeN;j += m){
						std::size_t k = j + mh;
						Ring x(w * getOut(k));
						getOut(k) = Ele(Ring(getOut(j)) - x);
						getOut(j) = Ele(Ring(getOut(j)) + x);
					}
					w *= w0;
				}
			}
			
			return ;
		}
		
	}; // namespace _utility
	
}; // namespace bignum
#endif // _BIG_NUM_FFT_HPP_