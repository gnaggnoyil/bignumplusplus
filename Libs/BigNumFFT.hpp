#ifndef _BIG_NUM_HPP_
#error "This header must be included through BigNum.hpp"
#endif // _BIG_NUM_HPP_

#ifndef _BIG_NUM_FFT_HPP_
#define _BIG_NUM_FFT_HPP_

#include <cstddef>
#include <cmath>
#include <utility>
#include <iterator>
#include <type_traits>

#include "BigNumTypeTrait.hpp"

namespace bignum{
	
	/*namespace _type{
		
		template <typename Ring>
		class RingTraits{
		public:
			template <typename RingRef1, typename RingRef2, 
				typename std::enable_if<isRLRef<Ring, RingRef1 &&>::value>::type * = nullptr, 
				typename std::enable_if<isRLRef<Ring, RingRef2 &&>::value>::type * = nullptr>
			inline static Ring add(RingRef1 &&_lhs, RingRef2 &&_rhs){
				return std::forward<RingRef1>(_lhs) + std::forward<RingRef2>(_rhs);
			}
			template <typename RingRef, 
				typename std::enable_if<isRLRef<Ring, RingRef &&>::value>::type * = nullptr>
			inline static void selfAdd(Ring &_lhs, RingRef &&_rhs){
				_lhs += std::foward<RingRef>(_rhs);
			}
			
			template <typename RingRef1, typename RingRef2, 
				typename std::enable_if<isRLRef<Ring, RingRef1 &&>::value>::type * = nullptr, 
				typename std::enable_if<isRLRef<Ring, RingRef2 &&>::value>::type * = nullptr>
			inline static Ring sub(RingRef1 &&_lhs, RingRef2 &&_rhs){
				return std::forward<RingRef1>(_lhs) - std::forward<RingRef2>(_rhs);
			}
			template <typename RingRef, 
				typename std::enable_if<isRLRef<Ring, RingRef &&>::value>::type * = nullptr>
			inline static void selfSub(Ring &_lhs, RingRef &&_rhs){
				_lhs -= std::foward<RingRef>(_rhs);
			}
			
			template <typename RingRef1, typename RingRef2, 
				typename std::enable_if<isRLRef<Ring, RingRef1 &&>::value>::type * = nullptr, 
				typename std::enable_if<isRLRef<Ring, RingRef2 &&>::value>::type * = nullptr>
			inline static Ring mult(RingRef1 &&_lhs, RingRef2 &&_rhs){
				return std::forward<RingRef1>(_lhs) * std::forward<RingRef2>(_rhs);
			}
			template <typename RingRef, 
				typename std::enable_if<isRLRef<Ring, RingRef &&>::value>::type * = nullptr>
			inline static void selfMult(Ring &_lhs, RingRef &&_rhs){
				_lhs *= std::foward<RingRef>(_rhs);
			}
			
			template <typename RingRef, 
				typename std::enable_if<isRLRef<Ring, RingRef &&>::value>::type * = nullptr>
			inline static Ring pow(RingRef &&_lhs, std::size_t _rhs){
				using std::pow;
				return pow(_lhs, _rhs);
			}
		};
		
	}; // namespace _type*/
	
	namespace _utility{
		
		/*template <typename InIter, typename RndIter, class Trait>
		inline static void fft1DPower2Impl(typename std::iterator_traits<RndIter>::size_type sizeN, 
			InIter in, RndIter out, Trait, 
			typename std::iterator_traits<RndIter>::value_type omega, 
			std::input_iterator_tag, std::random_access_iterator_tag){
			assert((sizeN & ((~sizeN) + 1)) == sizeN);
			
			using std::swap;
			
			using SizeT = typename std::iterator_traits<RndIter>::size_type;
			
			InIter inIter = in;
			RndIter outIter = out;
			for(SizeT i = 0;i < sizeN;++i){
				*outIter = *inIter;
				++outIter;
				++inIter;
			}
			
			SizeT preRev = 0; rev;
			for(SizeT i = 1;i < sizeN;++i){
				SizeT tmp = sizeN >> 1;
				for(rev = preRev;(rev & tmp) != 0;tmp >>= 1){
					rev -= tmp;
				}
				rev += tmp;
				preRev = rev;
				
				if(rev < i){
					swap(*(out + i), *(out + rev));
				}
			}
			
			using Ring = typename std::iterator_traits<RndIter>::value_type;
			
			for(SizeT m = 2;m <= sizeN;m <<= 1){
				Ring w0 = Trait::pow(omega, sizeN / m);
				Ring w = 1;
				for(SizeT i = 0, mh = m >> 1;i < mh;++i){
					for(SizeT j = i;j < sizeN;j += m){
						SizeT k = j + mh;
						Ring x = Trait::multi(w, *(out + k));
						*(out + k) = Trait::sub(*(out + j), x);
						*(out + j) = Trait::add(*(out + j), x);
					}
				}
				Trait::selfMult(w, w0);
			}
			
			return ;
		}
		
		template <typename InIter, typename RndIter>
		inline static void fft1DPower2(typename std::iterator_traits<RndIter>::size_type sizeN, 
			InIter in, RndIter out, 
			typename std::iterator_traits<RndIter>::value_type omega){
			using Ring = typename std::iterator_traits<RndIter>::value_type;
			fft1DPower2Impl(sizeN, in, out, omega, _type::RingTraits<Ring>{}, 
				typename std::iterator_traits<InIter>::iterator_category{}, 
				typename std::iterator_traits<RndIter>::iterator_category{});
		}
		
		template <typename InIter, typename RndIter, class Trait>
		inline static void fft1DPower2(typename std::iterator_traits<RndIter>::size_type sizeN, 
			InIter in, RndIter out, 
			typename std::iterator_traits<RndIter>::value_type omega, 
			Trait){
			fft1DPower2Impl(sizeN, in, out, omega, Trait{},  
				typename std::iterator_traits<InIter>::iterator_category{}, 
				typename std::iterator_traits<RndIter>::iterator_category{});
		}*/
		
		template <typename Ring, typename InFunc, typename OutFunc>
		inline static void fft1DPower2(std::size_t sizeN, Ring omega, InFunc &&getIn, OutFunc &&getOut){
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