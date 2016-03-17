#ifndef _BIG_NUM_HPP_
#error "This header must be included through BigNum.hpp"
#endif // _BIG_NUM_HPP_

#ifndef _BIG_NUM_MODULAR_RING_HPP_
#define _BIG_NUM_MODULAR_RING_HPP_

#include "BigNumTypeTrait.hpp"

namespace bignum{
	
	namespace _utility{
		
		template <typename SizeT, typename Ele, SizeT P, Ele OMEGA>
		class ModularP{
		private:
			using SqrEle = typename _type::squareType<Ele>::type;
		public:
			constexpr ModularP(Ele _num)
				:num(_num % P){}
			
			explicit constexpr operator Ele() const{
				return num;
			}
			
			constexpr ModularP &operator+=(const ModularP &_rhs){
				num = (num + _rhs.num) % P;
				return *this;
			}
			friend constexpr ModularP operator+(const ModularP &_lhs, const ModularP &_rhs){
				return ModularP((_lhs.num + _rhs.num) % P);
			}
			friend constexpr ModularP operator+(const ModularP &_lhs, Ele _rhs){
				return ModularP((_lhs.num + _rhs % P) % P);
			}
			friend constexpr ModularP operator+(Ele _lhs, const ModularP &_rhs){
				return ModularP((_lhs % P + _rhs.num) % P);
			}
			
			constexpr ModularP &operator-=(const ModularP &_rhs){
				num = (P - _rhs.num + num) % P;
				return *this;
			}
			friend constexpr ModularP operator-(const ModularP &_lhs, const ModularP &_rhs){
				return ModularP((P - _rhs.num + _lhs.num) % P);
			}
			friend constexpr ModularP operator-(const ModularP &_lhs, Ele _rhs){
				return ModularP((P - _rhs % P + _lhs.num) % P);
			}
			friend constexpr ModularP operator-(Ele _lhs, const ModularP &_rhs){
				return ModularP((P - _rhs.num + _lhs % P) % P);
			}
			
			constexpr ModularP &operator*=(const ModularP &_rhs){
				num = static_cast<Ele>((static_cast<SqrEle>(num) * static_cast<SqrEle>(_rhs.num)) % P);
				return *this;
			}
			friend constexpr ModularP operator*(const ModularP &_lhs, const ModularP &_rhs){
				return ModularP(static_cast<Ele>((static_cast<SqrEle>(_lhs.num) * static_cast<SqrEle>(_rhs.num)) % P));
			}
			friend constexpr ModularP operator*(const ModularP &_lhs, Ele _rhs){
				return ModularP(static_cast<Ele>((static_cast<SqrEle>(_lhs.num) * static_cast<SqrEle>(_rhs % P)) % P));
			}
			friend constexpr ModularP operator*(Ele _lhs, const ModularP &_rhs){
				return ModularP(static_cast<Ele>((static_cast<SqrEle>(_lhs % P) * static_cast<SqrEle>(_rhs.num)) % P));
			}
			
			friend constexpr ModularP pow(const ModularP &base, std::size_t exp){
				ModularP res(Ele(1));
				ModularP exponent = base;
				
				for(;exp > 0;exp >>= 1){
					if((exp & 1) == 1){
						res *= exponent;
					}
					exponent *= exponent;
				}
				
				return res;
			}
			
			/*constexpr ModularP inverse() const{
				//using std::swap;
				
				// exgcd
				std::intmax_t r(num), oldr(P);
				//std::intmax_t s(0), olds(1);
				std::intmax_t t(1), oldt(0);
				for(;r > 0;){
					// a*s+b*t==r
					//{
					//	std::intmax_t prov = olds - (oldr / r) * s;
					//	olds = s;
					//	s = prov;
					//}
					{
						std::intmax_t prov = oldt - (oldr / r) * t;
						oldt = t;
						t = prov;
					}
					{
						std::intmax_t prov = oldr % r;
						oldr = r;
						r = prov;
					}
				}
				// assert(oldr == 1)
				
				return ModularP((oldt < 0)? static_cast<Ele>(oldt + static_cast<std::intmax_t>(P)): static_cast<Ele>(oldt));
			}
			
			constexpr ModularP &operator/=(const ModularP &_rhs){
				num = static_cast<Ele>(num * _rhs.inverse());
				return *his;
			}
			friend constexpr ModularP operator/(const ModularP &_lhs, const ModularP &_rhs){
				return _lhs * _rhs.inverse();
			}
			friend constexpr ModularP operator/(const ModularP &_lhs, const Ele &_rhs){
				return _lhs * ModularP(_rhs).inverse();
			}
			friend constexpr ModularP operator/(const Ele &_lhs, const ModularP &_rhs){
				return _lhs * _rhs.inverse();
			}*/
			
			Ele num;
		};
		
	};// namespace _utility
	
}; // namespace bignum
#endif // _BIG_NUM_MODULAR_RING_HPP_