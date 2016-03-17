#ifndef _BIG_NUM_HPP_
#error "This header must be included through BigNum.hpp"
#endif // _BIG_NUM_HPP_

#ifndef _BIG_INT_HPP_
#define _BIG_INT_HPP_

#include <memory>
#include <utility>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <cassert>
#include <sstream>
#include <climits>
#include <ostream>
#include <iterator>
#include <locale>
#include <istream>
#include <string>

#include "Libs/BigNumTypeTrait.hpp"
#include "Libs/BigNumMemory.hpp"
#include "Libs/BigNumFFT.hpp"
#include "Libs/BigNumModularRing.hpp"
#include "Libs/BigNumGenerics.hpp"
#include "BigInt/BigIntOutput.hpp"
#include "BigInt/BigIntInput.hpp"

namespace bignum{
	
	using _type::isSigned;
	using _utility::destroyAll;
	using _type::isRLRef;
	using _utility::fft1DPower2;
	
	template <class Allocator = std::allocator<std::uint32_t>>
	class BigInt{
	private:
		using Alloc = Allocator;
		using Ptr = typename Alloc::pointer;
		
		using AllocTrait = std::allocator_traits<Alloc>;
		
		using SizeT = std::uint32_t;
		using LogSizeT = std::uint8_t;
		
		using Ele = typename Alloc::value_type;
		
		static constexpr SizeT MAX_LEN = 32768;
		static constexpr SizeT PRI_ORDER = 134217728;
		//static constexpr SizeT MAX_LEN = 134217728;	// 2 ^ 27
		static constexpr Ele P = 2013265921;	// 15 * MAX_LEN + 1
		static constexpr Ele OMEGA = 440564289;	// 31 ^ 15 mod P
		static constexpr LogSizeT ENTRY_SIZE = 8;
		//static constexpr LogSizeT ENTRY_SIZE = 15;	// floor(log(p) / 2)
		static constexpr Ele TWO_INV = 1006632961;	// 2 ^ (-1) mod P
		
		using ModularP_T = _utility::ModularP<SizeT, Ele, P, OMEGA>;
		
		// wrapper for Ptr to simpfy array operations
		struct DigitBuffer{
		public:
			DigitBuffer(Alloc *_alloc, SizeT _len)
				:alloc(_alloc), len(_len), cap(static_cast<SizeT>(std::pow(2.0, std::ceil(std::log2(len))))){
				data = alloc->allocate(static_cast<std::size_t>(cap));
			}
			DigitBuffer(Alloc *_alloc, SizeT _len, SizeT _cap)
				:alloc(_alloc), len(_len), cap(_cap){
				assert(cap == static_cast<SizeT>(std::pow(2.0, std::ceil(std::log2(len)))));
				data = alloc->allocate(static_cast<SizeT>(cap));
			}
			
			DigitBuffer(Alloc *_alloc, std::nullptr_t)
				:alloc(_alloc), data(nullptr), len(0), cap(0){}
			DigitBuffer(Alloc *_alloc, std::nullptr_t, SizeT _len, SizeT _cap)
				:alloc(_alloc), data(nullptr), len(_len), cap(_cap){}
			
			DigitBuffer(Alloc *_alloc, const DigitBuffer &_rhs)
				:alloc(_alloc), len(_rhs.len), cap(_rhs.cap){
				data = alloc->allocate(static_cast<std::size_t>(cap));
				SizeT i(0);
				try{
					for(;i < len;++i){
						alloc->construct(data + i, _rhs.data[i]);
					}
				}
				catch(...){
					destroyAll(data, data + i, *alloc);
					alloc->deallocate(data, static_cast<std::size_t>(cap));
					data = nullptr;
					throw ;
				}
			}
			DigitBuffer(Alloc *_alloc, DigitBuffer &&_rhs)
				:alloc(_alloc), data(std::move(_rhs.data)), len(std::move(_rhs.len)), cap(std::move(_rhs.cap)){}
			
			~DigitBuffer(){
				if(alloc == nullptr){
					assert(data == nullptr);
					return ;
				}
				if(data != nullptr){
					alloc->deallocate(data, static_cast<std::size_t>(cap));
					data = nullptr;
				}
				len = 0;
				cap = 0;
				alloc = nullptr;
			}
			
			DigitBuffer(const DigitBuffer &_rhs)
				:DigitBuffer(_rhs.alloc, _rhs){}
			
			DigitBuffer(DigitBuffer &&_rhs)
				:DigitBuffer(std::move(_rhs.alloc), std::move(_rhs)){}
			
			DigitBuffer &operator=(const DigitBuffer &) = delete;
			DigitBuffer &operator=(DigitBuffer &&) = delete;
			
			inline DigitBuffer &operator=(std::nullptr_t){
				destroyAll(data, data + len, *alloc);
				alloc->deallocate(data, static_cast<std::size_t>(cap));
				
				data = nullptr;
				len = 0;
				cap = 0;
				alloc = nullptr;
				return *this;
			}
			
			inline friend bool operator==(const DigitBuffer &_lhs, std::nullptr_t){
				if(_lhs.data != nullptr){
					return false;
				}
				//assert(len == 0);
				//assert(cap == 0);
				return true;
			}
			inline friend bool operator==(std::nullptr_t, const DigitBuffer &_rhs){
				return _rhs == nullptr;
			}
			inline friend bool operator!=(const DigitBuffer &_lhs, std::nullptr_t){
				return !(_lhs == nullptr);
			}
			inline friend bool operator!=(std::nullptr_t, const DigitBuffer &_rhs){
				return !(_rhs == nullptr);
			}
			
			inline void zeroLen() noexcept{
				len = 0;
				cap = 0;
			}
			
			inline DigitBuffer realloc(SizeT _len){
				SizeT _cap = static_cast<SizeT>(std::pow(2.0, std::ceil(std::log2(_len))));
				
				if(_len == len){
					return DigitBuffer(alloc, nullptr, _len, _cap);
				}
				
				if(_cap == cap){
					return DigitBuffer(alloc, nullptr, _len, _cap);
				}
				
				if(_cap > MAX_LEN){
					// maybe use some compile time techiques to caculate the string
					// literals needed
					throw std::out_of_range("BigInt::DigitBuffer::realloc");
				}
				
				return DigitBuffer(alloc, _len, _cap);
			}
			
			inline void reconstruct(DigitBuffer &tmp, SizeT _len, SizeT _cap){
				assert(static_cast<SizeT>(std::pow(2.0, std::ceil(std::log2(tmp.len)))) == tmp.cap);
				
				if(len == tmp.len){
					assert(tmp == nullptr);
					return ;
				}
				
				if(nullptr == tmp){
					// no realloc
					assert(_cap == cap);
					
					if(_len < len){
						// shrink
						destroyAll(data + _len, data + len, *alloc);
					}
					else{
						// expand
						SizeT i = len;
						try{
							for(;i < _len;++i){
								alloc->construct(data + i, Ele());
							}
						}
						catch(...){
							destroyAll(data + len, data + i, *alloc);
							throw;
						}
					}
				}
				else{
					assert(tmp.len == _len);
					assert(tmp.cap == _cap);
					assert(cap != _cap);
					
					if(_len < len){
						// shrink
						SizeT i(0);
						try{
							for(;i < _len;++i){
								tmp.alloc->construct(tmp.data + i, std::move(data[i]));
							}
						}
						catch(...){
							destroyAll(tmp.data, tmp.data + i, *(tmp.alloc));
							tmp.alloc->deallocate(tmp.data, tmp.cap);
							tmp.data = nullptr;
							throw ;
						}
					}
					else{
						// expand
						SizeT i(0);
						try{
							for(;i < len;++i){
								tmp.alloc->construct(tmp.data + i, std::move(data[i]));
							}
							for(;i < _len;++i){
								tmp.alloc->construct(tmp.data + i, Ele());
							}
						}
						catch(...){
							destroyAll(tmp.data, tmp.data + i, *(tmp.alloc));
							tmp.alloc->deallocate(tmp.data, tmp.cap);
							tmp.data = nullptr;
							throw ;
						}
					}
				}
			}
			
			inline void reassign(DigitBuffer &&tmp){
				if(nullptr != tmp){
					destroyAll(data, data + len, *alloc);
					alloc->deallocate(data, cap);
					
					data = std::move(tmp.data);
					alloc = tmp.alloc;
					
					tmp.data = nullptr;
				}
				
				len = std::move(tmp.len);
				cap = std::move(tmp.cap);
				
				tmp.zeroLen();
			}
			
			inline void resize(SizeT _len){
				DigitBuffer tmp = realloc(_len);
				reconstruct(tmp, _len, tmp.cap);
				reassign(std::move(tmp));
			}
			
			inline void setLen(SizeT _len) noexcept{
				len = _len;
				cap = static_cast<SizeT>(std::pow(2.0, std::ceil(std::log2(len))));
			}
			
			inline void setCap() noexcept{
				cap = static_cast<SizeT>(std::pow(2.0, std::ceil(std::log2(len))));
			}
			
			// 1 - thisRaw > _rhsRaw
			// 0 - thisRaw == _rhsRaw
			// -1 - thisRaw < _rhsRaw
			std::int8_t compareRaw(const DigitBuffer &_rhs) const{
				if(len > _rhs.len){
					return 1;
				}
				if(len < _rhs.len){
					return -1;
				}
				
				assert(len == _rhs.len);
				for(SizeT i(1);i <= len;++i){
					if(data[len - i] > _rhs.data[len - i]){
						return 1;
					}
					if(data[len - i] < _rhs.data[len - i]){
						return -1;
					}
				}
				
				return 0;
			}
			
			template <typename UnsignedInt, SizeT dLen>
			inline std::int8_t compareUnsignedIntBuffer(const UnsignedInt &_rhs, std::integral_constant<SizeT, dLen>) const{
				if(len > dLen){
					return 1;
				}
				
				return compareUnsignedIntRaw(_rhs, std::integral_constant<SizeT, dLen>());
			}
			// overload for small unsigned types
			template <typename UnsignedInt>
			inline std::int8_t compareUnsignedIntBuffer(const UnsignedInt &_rhs, std::integral_constant<SizeT, 1>) const{
				if(data[0] < _rhs){
					return -1;
				}
				if(data[0] > _rhs){
					return 1;
				}
				return 0;
			}
			
			Ele propagateCarryRange(SizeT st, SizeT en){
				Ele carry(0);
				for(SizeT i(st);i < en;++i){
					data[i] += carry;
					carry = data[i] >> ENTRY_SIZE;
					data[i] -= carry << ENTRY_SIZE;
				}
				return carry;
			}
			
			inline void propagateCarry(){
				Ele carry = propagateCarryRange(0, len);
				
				SizeT expand(0);
				for(Ele tmp = carry;tmp > 0;tmp >>= ENTRY_SIZE, ++expand);
				SizeT _len = len;
				// simple hack. TODO: manually resize and construct carries
				resize(len + expand);
				for(SizeT i = _len;i < len;++i){
					data[i] = carry & ((1 << ENTRY_SIZE) - 1);
					carry >>= ENTRY_SIZE;
				}
			}
			
			// thisRaw += _rhsRaw
			inline void addRaw(const DigitBuffer &_rhs){
				// simple hack
				if(len < _rhs.len){
					resize(_rhs.len);
				}
				for(SizeT i(0);i < _rhs.len;++i){
					data[i] += _rhs.data[i];
				}
				propagateCarry();
			}
			inline void addRaw(DigitBuffer &&_rhs){
				if(len < _rhs.len){
					resize(_rhs.len);
				}
				for(SizeT i = 0;i < _rhs.len;++i){
					data[i] += std::move(_rhs.data[i]);
				}
				propagateCarry();
			}
			
			inline void shrinkToFit(){
				SizeT _len = len - 1;
				for(;true;--_len){
					if(Ele(0) != data[_len]){
						break;
					}
					
					if(_len == 0){
						break;
					}
				}
				++_len;
				
				resize(_len);
			}
			
			// thisRaw -= _rhsRaw
			// assert(thisRaw >= _rhsRaw)
			void subRaw(const DigitBuffer &_rhs){
				Ele carry(0);
				assert(len >= _rhs.len);
				
				for(SizeT i(0);i < _rhs.len;++i){
					assert(data[i] + (1 << ENTRY_SIZE) >= _rhs.data[i] + carry);
					if(data[i] < (_rhs.data[i] + carry)){
						data[i] += (1 << ENTRY_SIZE) - (_rhs.data[i] + carry);
						carry = Ele(1);
					}
					else{
						data[i] -= _rhs.data[i] + carry;
						carry = Ele(0);
					}
				}
				
				for(SizeT i = _rhs.len;i < len;++i){
					assert(data[i] + (1 << ENTRY_SIZE) >= carry);
					if(data[i] < carry){
						data[i] += (1 << ENTRY_SIZE) - carry;
						carry = Ele(1);
					}
					else{
						data[i] -= carry;
						carry = Ele(0);
						break;
					}
				}
				if(Ele(0) != carry){
					throw std::underflow_error("left operand less than right operand.");
				}
				
				shrinkToFit();
			}
			void subRaw(DigitBuffer &&_rhs){
				Ele carry(0);
				assert(len >= _rhs.len);
				
				for(SizeT i(0);i < _rhs.len;++i){
					assert(data[i] + (1 << ENTRY_SIZE) >= _rhs.data[i] + carry);
					if(data[i] < (_rhs.data[i] + carry)){
						data[i] += (1 << ENTRY_SIZE) - (std::move(_rhs.data[i]) + carry);
						carry = Ele(1);
					}
					else{
						data[i] -= std::move(_rhs.data[i]) + carry;
						carry = Ele(0);
					}
				}
				
				for(SizeT i = _rhs.len;i < len;++i){
					assert(data[i] + (1 << ENTRY_SIZE) >= carry);
					if(data[i] < carry){
						data[i] += (1 << ENTRY_SIZE) - carry;
						carry = Ele(1);
					}
					else{
						data[i] -= carry;
						carry = Ele(0);
						break;
					}
				}
				if(Ele(0) != carry){
					throw std::underflow_error("left operand less than right operand.");
				}
				
				shrinkToFit();
			}
		private:
			// TODO: rewite this function using std::integer_sequence instead of raw
			// "recursion"
			template <typename UnsignedInt, SizeT dLen>
			inline std::int8_t compareUnsignedIntRaw(const UnsignedInt &_rhs, std::integral_constant<SizeT, dLen>) const{
				Ele digit = (_rhs >> ((dLen - 1) * ENTRY_SIZE)) & ((1 << ENTRY_SIZE) - 1);
				if((digit > 0)){
					if(len < dLen){
						return -1;
					}
					if(data[dLen - 1] < digit){
						return -1;
					}
					if(data[dLen - 1] > digit){
						return 1;
					}
					return compareUnsignedIntRaw(_rhs, std::integral_constant<SizeT, dLen - 1>());
				}
				else{
					if(len < dLen){
						return compareUnsignedIntRaw(_rhs, std::integral_constant<SizeT, dLen - 1>());
					}
					if(data[dLen - 1] > 0){
						return 1;
					}
					return compareUnsignedIntRaw(_rhs, std::integral_constant<SizeT, dLen - 1>());
				}
			}
			// note: this is an overload.
			template <typename UnsignedInt>
			inline std::int8_t compareUnsignedIntRaw(const UnsignedInt &_rhs, std::integral_constant<SizeT, 0>) const{
				if(data[0] < _rhs){
					return -1;
				}
				if(data[0] > _rhs){
					return 1;
				}
				return 0;
			}
		public:			
			Alloc *alloc;
			Ptr data;
			SizeT len, cap;
		};// struct DigitBuffer
	public:
		// I want a more specialized version...
		template <typename, class, typename>
		friend class _GenericRadix;
		template <typename, class, typename>
		friend class _DecimalRadix;
		template <typename, class, typename>
		friend class _SmallPower2Radix;
		template <typename, class, typename>
		friend class _LargePower2Radix;
		template <typename, class>
		friend class _ExactDigitExtract;
		
		template <typename, class>
		friend class RadixConvertEnumer;
		
		template <typename, class, typename>
		friend class _DecimalAutomatic;
		template <typename, class, typename>
		friend class _GenericAutomatic;
		template <typename, class, typename, class>
		friend class _SmallPower2RadixAutomatic;
		template <typename, class, typename, class>
		friend class _LargePower2RadixAutomatic;
		template <typename, class, typename, class>
		friend class _ExactDigitAutomatic;
		
		template <typename, class>
		friend class RadixConvertRecver;
		
		template <typename, class>
		friend class _type::DigitRecvIterator;
		
		template <class, class>
		friend class _type::LiteralParser;
		
		template <typename Digit>
		using RadixConvertEnumer = RadixConvertEnumer<Digit, BigInt>;
		
		template <typename Digit>
		using RadixConvertRecver = RadixConvertRecver<Digit, BigInt>;
		
		template <typename Digit>
		using DigitRecvIterator = _type::DigitRecvIterator<Digit, BigInt>;
		
		// default zero instead of nullptr, since BigInt should behave as normal 
		// numerial, not accepting null state
		explicit BigInt()
			:allocator(), buf(&allocator, 1), positive(true){
			allocator.construct(buf.data + 0, Ele(0));
		}
		
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		/*explicit */BigInt(Integer _rhs)
			:allocator(), buf(&allocator, nullptr){
			assignIntegral(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		
		// copy constructor
		BigInt(const BigInt &_rhs)
			:allocator(AllocTrait::select_on_container_copy_construction(_rhs.allocator)), buf(&allocator, nullptr, _rhs.buf.len, _rhs.buf.cap), positive(_rhs.positive){
			assert(buf.cap <= MAX_LEN);
			
			buf.data = allocator.allocate(static_cast<std::size_t>(buf.cap));
			SizeT i(0);
			try{
				for(;i < buf.len;++i){
					allocator.construct(buf.data + i, _rhs.buf.data[i]);
				}
			}
			catch(...){
				destroyAll(buf.data, buf.data + i, allocator);
				allocator.deallocate(buf.data, buf.cap);
				buf.data = nullptr;
				throw ;
			}
		}
		
		// move constructor
		BigInt(BigInt &&_rhs)
			:allocator(std::move(_rhs.allocator)), buf(&allocator, std::move(_rhs.buf)), positive(std::move(_rhs.positive)){
			_rhs.buf.data = nullptr;
			_rhs.buf.len = 0;
			_rhs.buf.cap = 0;
		}
		
		template <typename Digit, 
			typename std::enable_if<std::is_integral<Digit>::value>::type * = nullptr>
		explicit BigInt(const DigitRecvIterator<Digit> &_rhs)
			:BigInt(_rhs.receiver->_finish()){}
		template <typename Digit, 
			typename std::enable_if<std::is_integral<Digit>::value>::type * = nullptr>
		explicit BigInt(DigitRecvIterator<Digit> &&_rhs)
			:BigInt(std::move(*(_rhs.receiver))._finish()){}
		
		// destructor
		~BigInt(){
			if(nullptr == buf.data){
				// This object has been moved
				buf.len = 0;
				buf.cap = 0;
				return ;
			}
			
			destroyAll(buf.data, buf.data + buf.len, allocator);
			allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
			buf.data = nullptr;
			buf.len = 0;
			buf.cap = 0;
		}
		
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		BigInt &operator=(Integer _rhs){
			if(nullptr != buf.data){
				destroyAll(buf.data, buf.data + buf.len, allocator);
				allocator.deallocate(buf.data, buf.cap);
				buf.data = nullptr;
				buf.len = 0;
				buf.cap = 0;
			}
			
			assignIntegral(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
			return *this;
		}
		
		// copy assignment operator
		BigInt &operator=(const BigInt &_rhs){
			if(this != &_rhs){
				assignLv(_rhs, typename AllocTrait::propagate_on_container_copy_assignment());
				
				assert(buf.cap <= MAX_LEN);
			}
			
			return *this;
		}
		
		// move assignment operator
		BigInt &operator=(BigInt &&_rhs){
			if(this != &_rhs){
				assignRv(std::move(_rhs), typename AllocTrait::propagate_on_container_move_assignment());
				
				assert(buf.cap <= MAX_LEN);
			}
			
			return *this;
		}
		
		inline const BigInt &operator+() const{
			return *this;
		}
		
		inline BigInt operator-() const &{
			BigInt tmp = *this;
			tmp.changeSign();
			return tmp;
		}
		inline BigInt operator-() &&{
			changeSign();
			return *this;
		}
		
		// swap two BigInts
		inline void swap(BigInt &_rhs){
			if(this != &_rhs){
				swapImpl(_rhs, typename AllocTrait::propagate_on_container_swap());
			}
			
			return ;
		}
		
		// is equal
		inline friend bool operator==(const BigInt &_lhs, const BigInt &_rhs){
			return 0 == _lhs.compare(_rhs);
		}
		// less than
		inline friend bool operator<(const BigInt &_lhs, const BigInt &_rhs){
			return -1 == _lhs.compare(_rhs);
		}
		// greater than
		inline friend bool operator>(const BigInt &_lhs, const BigInt &_rhs){
			return 1 == _lhs.compare(_rhs);
		}
		inline friend bool operator!=(const BigInt &_lhs, const BigInt &_rhs){
			return 0 != _lhs.compare(_rhs);
		}
		inline friend bool operator<=(const BigInt &_lhs, const BigInt &_rhs){
			return _lhs.compare(_rhs) <= 0;
		}
		inline friend bool operator>=(const BigInt &_lhs, const BigInt &_rhs){
			return _lhs.compare(_rhs) >= 0;
		}
		
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline friend bool operator==(const BigInt &_lhs, Integer _rhs){
			return 0 == _lhs.compareInt(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline friend bool operator<(const BigInt &_lhs, Integer _rhs){
			return -1 == _lhs.compareInt(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline friend bool operator>(const BigInt &_lhs, Integer _rhs){
			return 1 == _lhs.compareInt(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline friend bool operator!=(const BigInt &_lhs, Integer _rhs){
			return 0 != _lhs.compareInt(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline friend bool operator<=(const BigInt &_lhs, Integer _rhs){
			return 0 > _lhs.compareInt(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline friend bool operator>=(const BigInt &_lhs, Integer _rhs){
			return 0 <= _lhs.compareInt(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline friend bool operator==(Integer &_lhs, const BigInt &_rhs){
			return 0 == _rhs.compareInt(_lhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline friend bool operator<(Integer &_lhs, const BigInt &_rhs){
			return 1 == _rhs.compareInt(_lhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline friend bool operator>(Integer &_lhs, const BigInt &_rhs){
			return -1 == _rhs.compareInt(_lhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline friend bool operator!=(Integer &_lhs, const BigInt &_rhs){
			return 0 != _rhs.compareInt(_lhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline friend bool operator<=(Integer &_lhs, const BigInt &_rhs){
			return 0 <= _rhs.compareInt(_lhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline friend bool operator>=(Integer &_lhs, const BigInt &_rhs){
			return 0 >= _rhs.compareInt(_lhs, std::integral_constant<bool, isSigned<Integer>::value>());
		}
		
		// shl
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		BigInt &operator<<=(const Integer &_rhs){
			shl(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
			return *this;
		}
		
		template <class BigIntRef, typename Integer, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef &&>::value && std::is_integral<Integer>::value>::type * = nullptr>
		friend BigInt operator<<(BigIntRef &&_lhs, const Integer &_rhs){
			BigInt tmp = std::forward<BigIntRef &&>(_lhs);
			tmp.shl(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
			return tmp;
		}
		
		// shr
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		BigInt &operator>>=(const Integer &_rhs){
			shl(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
			return *this;
		}
		
		template <typename BigIntRef, typename Integer, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef &&>::value && std::is_integral<Integer>::value>::type * = nullptr>
		friend BigInt operator>>(BigIntRef &&_lhs, const Integer &_rhs){
			BigInt tmp = std::forward<BigIntRef>(_lhs);
			tmp.shr(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
			return tmp;
		}
		
		// self add
		template <typename BigIntRef, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef &&>::value>::type * = nullptr>
		BigInt &operator+=(BigIntRef &&_rhs){
			if(this != &_rhs){
				add(std::forward<BigIntRef>(_rhs));
				return *this;
			}
			else{
				// this <<= 1
				shl(static_cast<unsigned int>(1), std::false_type());
				return *this;
			}
		}
		
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline BigInt &operator+=(const Integer &_rhs){
			operator+=(BigInt(_rhs));
			return *this;
		}
		
		inline friend BigInt operator+(const BigInt &_lhs, const BigInt &_rhs){
			BigInt tmp = _lhs;
			tmp += _rhs;
			return tmp;
		}
		inline friend BigInt operator+(const BigInt &_lhs, BigInt &&_rhs){
			BigInt tmp = std::move(_rhs);
			tmp += _lhs;
			return tmp;
		}
		inline friend BigInt operator+(BigInt &&_lhs, const BigInt &_rhs){
			BigInt tmp = std::move(_lhs);
			tmp += _rhs;
			return  tmp;
		}
		inline friend BigInt operator+(BigInt &&_lhs, BigInt &&_rhs){
			BigInt tmp = std::move(_lhs);
			tmp += _rhs;
			return tmp;
		}
		template <class BigIntRef, typename Integer, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef &&>::value && std::is_integral<Integer>::value>::type * = nullptr>
		inline friend BigInt operator+(BigIntRef &&_lhs, const Integer &_rhs){
			BigInt tmp = std::forward<BigIntRef>(_lhs);
			tmp += _rhs;
			return tmp;
		}
		template <class BigIntRef, typename Integer, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef &&>::value && std::is_integral<Integer>::value>::type * = nullptr>
		inline friend BigInt operator+(const Integer &_lhs, BigIntRef &&_rhs){
			BigInt tmp = std::forward<BigIntRef>(_rhs);
			tmp += _lhs;
			return tmp;
		}
		
		// self substract
		template <class BigIntRef, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef &&>::value>::type * = nullptr>
		inline BigInt &operator-=(BigIntRef &&_rhs){
			if(this != &_rhs){
				sub(std::forward<BigIntRef>(_rhs));
			}
			else{
				zerolize();
			}
			
			return *this;
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline BigInt &operator-=(const Integer &_rhs){
			sub(BigInt(_rhs));
			return *this;
		}
		
		inline friend BigInt operator-(const BigInt &_lhs, const BigInt &_rhs){
			BigInt tmp = _lhs;
			tmp -= _rhs;
			return tmp;
		}
		inline friend BigInt operator-(const BigInt &_lhs, BigInt &&_rhs){
			BigInt tmp = std::move(_rhs);
			tmp.changeSign();
			tmp += _lhs;
			return tmp;
		}
		inline friend BigInt operator-(BigInt &&_lhs, const BigInt &_rhs){
			BigInt tmp = std::move(_lhs);
			tmp -= _rhs;
			return tmp;
		}
		inline friend BigInt operator-(BigInt &&_lhs, BigInt &&_rhs){
			BigInt tmp = std::move(_lhs);
			tmp -= std::move(_rhs);
			return tmp;
		}
		template <class BigIntRef, typename Integer, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef &&>::value && std::is_integral<Integer>::value>::type * = nullptr>
		inline friend BigInt operator-(BigIntRef &&_lhs, const Integer &_rhs){
			BigInt tmp = std::forward<BigIntRef>(_lhs);
			tmp -= _rhs;
			return tmp;
		}
		template <class BigIntRef, typename Integer, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef &&>::value && std::is_integral<Integer>::value>::type * = nullptr>
		inline friend BigInt operator-(const Integer &_lhs, BigIntRef &&_rhs){
			BigInt tmp = std::forward<BigIntRef>(_rhs);
			tmp += _lhs;
			return tmp;
		}
		
		// self multiply
		template <class BigIntRef, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef &&>::value>::type * = nullptr>
		inline BigInt &operator*=(BigIntRef &&_rhs){
			multiply(std::forward<BigIntRef>(_rhs));
			
			return *this;
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		inline BigInt &operator*=(const Integer &_rhs){
			multiplyInt(_rhs, std::integral_constant<bool, isSigned<Integer>::value>());
			return *this;
		}
		
		inline friend BigInt operator*(const BigInt &_lhs, const BigInt &_rhs){
			BigInt tmp = _lhs;
			tmp *= _rhs;
			return tmp;
		}
		inline friend BigInt operator*(const BigInt &_lhs, BigInt &&_rhs){
			BigInt tmp = std::move(_rhs);
			tmp *= _lhs;
			return tmp;
		}
		inline friend BigInt operator*(BigInt &&_lhs, const BigInt &_rhs){
			BigInt tmp = std::move(_lhs);
			tmp *= _rhs;
			return tmp;
		}
		inline friend BigInt operator*(BigInt &&_lhs, BigInt &&_rhs){
			BigInt tmp = std::move(_lhs);
			tmp *= std::move(_rhs);
			return tmp;
		}
		
		// self divide
		inline BigInt &operator/=(BigInt &_rhs){
			if(this != &_rhs){
				bool _positive1 = positive;
				bool _positive2 = _rhs.positive;
				positive = true;
				_rhs.positive = true;
				*this = std::move(*this).divideBy(_rhs).first;
				positive = _positive1 == _positive2;
				_rhs.positive = _positive2;
				return *this;
			}
			else{
				destroyAll(buf.data, buf.data + buf.len, allocator);
				allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
				buf.data = allocator.allocate(1);
				allocator.construct(buf.data + 0, 1);
				buf.len = 1;
				buf.cap = 1;
				return *this;
			}
		}
		inline BigInt &operator/=(BigInt &&_rhs){
			bool _positive = positive == _rhs.positive;
			positive = true;
			_rhs.positive = true;
			*this = std::move(*this).divideBy(std::move(_rhs)).first;
			positive = _positive;
			return *this;
		}
		inline BigInt &operator/=(const BigInt &_rhs){
			if(this != &_rhs){
				if(_rhs.positive){
					bool _positive = positive;
					positive = true;
					*this = std::move(*this).divideBy(_rhs).first;
					positive = _positive;
					return *this;
				}
				else{
					bool _positive = !positive;
					positive = true;
					BigInt tmp = _rhs;
					tmp.changeSign();
					*this = std::move(*this).divideBy(tmp).first;
					positive = _positive;
					return *this;
				}
			}
			else{
				destroyAll(buf.data, buf.data + buf.len, allocator);
				allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
				buf.data = allocator.allocate(1);
				allocator.construct(buf.data + 0, 1);
				buf.len = 1;
				buf.cap = 1;
				return *this;
			}
		}
		
		// TODO: avoid unnecessary memory allocation and construction
		template <class BigIntRef1, class BigIntRef2, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef1 &&>::value && isRLRef<BigInt, BigIntRef2 &&>::value>::type * = nullptr>
		inline friend BigInt operator/(BigIntRef1 &&_lhs, BigIntRef2 &&_rhs){
			BigInt tmp = std::forward<BigIntRef1>(_lhs);
			tmp /= std::forward<BigIntRef2>(_rhs);
			return tmp;
		}
		
		// self modular
		inline BigInt &operator%=(BigInt &&_rhs){
			bool _positive = positive == _rhs.positive;
			positive = true;
			_rhs.positive = true;
			*this = std::move(*this).modularBy(std::move(_rhs)).first;
			positive = _positive;
			return *this;
		}
		inline BigInt &operator%=(BigInt &_rhs){
			if(this != &_rhs){
				bool _positive1 = positive;
				bool _positive2 = _rhs.positive;
				positive = true;
				_rhs.positive = true;
				*this = std::move(*this).modularBy(_rhs).first;
				positive = _positive1 == _positive2;
				_rhs.positive = _positive2;
				return *this;
			}
			else{
				zerolize();
				return *this;
			}
		}
		inline BigInt &operator%=(const BigInt &_rhs){
			if(this != &_rhs){
				if(_rhs.positive){
					bool _positive = positive;
					positive = true;
					*this = std::move(*this).modularBy(_rhs).first;
					positive = _positive;
					return *this;
				}
				else{
					bool _positive = !positive;
					positive = true;
					BigInt tmp = _rhs;
					tmp.changeSign();
					*this = std::move(*this).modularBy(tmp).first;
					positive = _positive;
					return *this;
				}
			}
			else{
				zerolize();
				return *this;
			}
		}
		
		template <class BigIntRef1, class BigIntRef2, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef1 &&>::value && isRLRef<BigInt, BigIntRef2 &&>::value>::type * = nullptr>
		inline friend BigInt operator%(BigIntRef1 &&_lhs, BigIntRef2 &&_rhs){
			BigInt tmp = std::forward<BigIntRef1>(_lhs);
			tmp %= std::forward<BigIntRef2>(_rhs);
			return tmp;
		}
		
		// input
		template <typename Char, class Trait>
		friend std::basic_istream<Char, Trait> &operator>>(std::basic_istream<Char, Trait> &is, BigInt &_rhs){
			using Tr = typename std::basic_istream<Char, Trait>::traits_type;
			using Int = typename std::basic_istream<Char, Trait>::int_type;
			
			using Finder = _utility::CharFindHelper<Char, Int, Trait>;
			
			_rhs.buf.resize(1);
			_rhs.buf.data[0] = 0;
			_rhs.positive = true;
			
			// this check is useless for Char = char/wchar_t
			/*if(std::has_facet<std::ctype<Char>>(is.getloc())){
				
			}*/
			std::locale inLocale = is.getloc();
			auto &ctypeFacet = std::use_facet<std::ctype<Char>>(inLocale);
			auto &npFacet = std::use_facet<std::numpunct<Char>>(inLocale);
			
			if((is.flags() & std::ios_base::skipws) != 0){
				is >> std::ws;
			}
			Int _st;
			_st = is.peek();
			if(Tr::eq_int_type(_st, Tr::eof())){
				is.setstate(std::ios_base::eofbit);
				return is;
			}
			
			if((is.flags() & std::ios_base::oct) != 0){
				bool _positive = true;
				
				if(Tr::eq_int_type(Tr::to_int_type(_BIG_NUM_GENERIC_LITERAL_(Char, '-')), _st)){
					bool _positive = true;
					is.ignore();
					_st = is.peek();
					if(Tr::eq_int_type(_st, Tr::eof())){
						is.setstate(std::ios_base::eofbit);
						is.setstate(std::ios_base::failbit);
						return is;
					}
				}
				else if(Tr::eq_int_type(Tr::to_int_type(_BIG_NUM_GENERIC_LITERAL_(Char, '+')), _st)){
					is.ignore();
					_st = is.peek();
					if(Tr::eq_int_type(_st, Tr::eof())){
						is.setstate(std::ios_base::eofbit);
						is.setstate(std::ios_base::failbit);
						return is;
					}
				}
				
				static const Char *octDigit = _BIG_NUM_GENERIC_LITERAL_(Char, "01234567");
				using Automatic = typename _type::CompareCond<LogSizeT, 3, ENTRY_SIZE, 
					_LargePower2RadixAutomatic<Ele, BigInt, SizeT, std::input_iterator_tag>, 
					_ExactDigitAutomatic<Ele, BigInt, SizeT, std::input_iterator_tag>, 
					_SmallPower2RadixAutomatic<Ele, BigInt, SizeT, std::input_iterator_tag>>::type;
				Automatic automatic(3);
				_ThousandSepParser<Char, Trait> _parser(is, octDigit, octDigit + 8, ctypeFacet, npFacet);
				const Char *digitPtr = _parser.tryNextDigit();
				while(digitPtr == octDigit){
					digitPtr = _parser.tryNextDigit();
				}
				if((digitPtr < octDigit) || (digitPtr >= octDigit + 8)){
					return is;
				}
				while((digitPtr >= octDigit) && (digitPtr < octDigit + 8)){
					automatic._readDigit(digitPtr - octDigit);
					digitPtr = _parser.tryNextDigit();
				}
				_rhs = std::move(automatic)._finish();
				if((!_positive) && (!_rhs.isZero())){
					_rhs.positive = false;
				}
				return is;
			}
			
			if((is.flags() & std::ios_base::hex) != 0){
				bool _positive = true;
				
				if(Tr::eq_int_type(Tr::to_int_type(_BIG_NUM_GENERIC_LITERAL_(Char, '-')), _st)){
					_positive = false;
					is.ignore();
					_st = is.peek();
					if(Tr::eq_int_type(_st, Tr::eof())){
						is.setstate(std::ios_base::eofbit);
						is.setstate(std::ios_base::failbit);
						return is;
					}
				}
				else if(Tr::eq_int_type(Tr::to_int_type(_BIG_NUM_GENERIC_LITERAL_(Char, '+')), _st)){
					is.ignore();
					_st = is.peek();
					if(Tr::eq_int_type(_st, Tr::eof())){
						is.setstate(std::ios_base::eofbit);
						is.setstate(std::ios_base::failbit);
						return is;
					}
				}
				
				static const Char *hexDigit = _BIG_NUM_GENERIC_LITERAL_(Char, "0123456789abcdefABCDEF");
				using Automatic = typename _type::CompareCond<LogSizeT, 4, ENTRY_SIZE, 
					_LargePower2RadixAutomatic<Ele, BigInt, SizeT, std::input_iterator_tag>, 
					_ExactDigitAutomatic<Ele, BigInt, SizeT, std::input_iterator_tag>, 
					_SmallPower2RadixAutomatic<Ele, BigInt, SizeT, std::input_iterator_tag>>::type;
				Automatic automatic(4);
				_ThousandSepParser<Char, Trait> _parser(is, hexDigit, hexDigit + 16, ctypeFacet, npFacet);
				if(Tr::eq_int_type(Tr::to_int_type(_BIG_NUM_GENERIC_LITERAL_(Char, '0')), _st)){
					is.ignore();
					_st = is.peek();
					if(Tr::eq_int_type(_st, Tr::eof())){
						is.setstate(std::ios_base::eofbit);
						return is;
					}
					if(ctypeFacet.is(std::ctype_base::space, Tr::to_char_type(_st))){
						return is;
					}
					if((Tr::eq_int_type(Tr::to_int_type(_BIG_NUM_GENERIC_LITERAL_(Char, 'x')), _st)) ||
						(Tr::eq_int_type(Tr::to_int_type(_BIG_NUM_GENERIC_LITERAL_(Char, 'X')), _st))){
						is.ignore();
					}
					else{
						if(Tr::eq_int_type(Tr::to_int_type(npFacet.thousands_sep()), _st)){
							_parser.checkInput(Tr::to_int_type(_BIG_NUM_GENERIC_LITERAL_(Char, '0')));
						}
						const Char *digitPtr = Finder::find(hexDigit, hexDigit + 22, _st);
						if((digitPtr < hexDigit) || (digitPtr >= hexDigit + 22)){
							is.setstate(std::ios_base::failbit);
							return is;
						}
					}
				}
				
				const Char *digitPtr = _parser.tryNextDigit();
				while(digitPtr == hexDigit){
					digitPtr = _parser.tryNextDigit();
				}
				if((digitPtr < hexDigit) || (digitPtr >= hexDigit + 22)){
					return is;
				}
				
				while((digitPtr >= hexDigit) && (digitPtr < hexDigit + 22)){
					if(digitPtr >= hexDigit + 16){
						digitPtr -= 6;
					}
					automatic._readDigit(digitPtr - hexDigit);
					digitPtr = _parser.tryNextDigit();
				}
				_rhs = std::move(automatic)._finish();
				if((!_positive) && (!_rhs.isZero())){
					_rhs.positive = false;
				}
				return is;
			}
			
			bool _positive = true;
			if(Tr::eq_int_type(Tr::to_int_type(_BIG_NUM_GENERIC_LITERAL_(Char, '-')), _st)){
				_positive = false;
				is.ignore();
				_st = is.peek();
				if(Tr::eq_int_type(_st, Tr::eof())){
					is.setstate(std::ios_base::eofbit);
					is.setstate(std::ios_base::failbit);
					return is;
				}
			}
			else if(Tr::eq_int_type(Tr::to_int_type(_BIG_NUM_GENERIC_LITERAL_(Char, '+')), _st)){
				is.ignore();
				_st = is.peek();
				if(Tr::eq_int_type(_st, Tr::eof())){
					is.setstate(std::ios_base::eofbit);
					is.setstate(std::ios_base::failbit);
					return is;
				}
			}
			
			static const Char *decDigit = _BIG_NUM_GENERIC_LITERAL_(Char, "0123456789");
			_DecimalAutomatic<Ele, BigInt, SizeT> automatic;
			_ThousandSepParser<Char, Trait> _parser(is, decDigit, decDigit + 10, ctypeFacet, npFacet);
			const Char *digitPtr = _parser.tryNextDigit();
			while(digitPtr == decDigit){
				digitPtr = _parser.tryNextDigit();
			}
			if((digitPtr < decDigit) || (digitPtr >= decDigit + 10)){
				return is;
			}
			while((digitPtr >= decDigit) && (digitPtr < decDigit + 10)){
				automatic._readDigit(digitPtr - decDigit);
				digitPtr = _parser.tryNextDigit();
			}
			_rhs = std::move(automatic)._finish();
			if((!_positive) && (!_rhs.isZero())){
				_rhs.positive = false;
			}
			return is;
		}
		

		
		template <typename Integer, typename Iter, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		friend BigInt getDigit(Integer radix, Iter begin, Iter end){
			RadixConvertRecver<Integer> recver(radix);
			return recver.readDigits(begin, end);
		}
		
		// output		
		// implicitly inlined
		// TODO: output thousand seperators correctly for decimal output according to the
		// current locale that os is using. This requires a precomputation of a total length
		// of the output, which will be avaliable as soon as the first(highest) digit is
		// calculated.
		template <typename Char, class Trait>
		friend std::basic_ostream<Char, Trait> &operator<<(std::basic_ostream<Char, Trait> &os, BigInt &&_rhs){
			if(!_rhs.positive){
				assert(!_rhs.isZero());
				os.put(_BIG_NUM_GENERIC_LITERAL_(Char, '-'));
				_rhs.positive = true;
			}
			else{
				if((os.flags() & std::ios_base::showpos) != 0){
					os.put(_BIG_NUM_GENERIC_LITERAL_(Char, '+'));
				}
			}
			
			if((os.flags() & std::ios_base::oct) != 0){
				if((os.flags() & std::ios_base::showbase) != 0){
					os.put(_BIG_NUM_GENERIC_LITERAL_(Char, '0'));
				}
				
				if(_rhs.isZero()){
					os.put(_BIG_NUM_GENERIC_LITERAL_(Char, '0'));
					return os;
				}
				
				static const Char *octDigit = _BIG_NUM_GENERIC_LITERAL_(Char, "01234567");
				using Producer = typename _type::CompareCond<LogSizeT, 3, ENTRY_SIZE, 
					_LargePower2Radix<SizeT, BigInt>, 
					_ExactDigitExtract<SizeT, BigInt>, 
					_SmallPower2Radix<SizeT, BigInt>>::type;
				Producer producer(std::move(_rhs), 3);
				for(SizeT digit = producer._start();true;digit = producer._next()){
					assert(digit < 8);
					os.put(octDigit[digit]);
					if(!producer._hasNext()){
						break;
					}
				}
				
				return os;
			}
			
			if((os.flags() & std::ios_base::hex) != 0){
				if((os.flags() & std::ios_base::showbase) != 0){
					if((os.flags() & std::ios_base::uppercase) != 0){
						os.write(_BIG_NUM_GENERIC_LITERAL_(Char, "0X"), 2);
					}
					else{
						os.write(_BIG_NUM_GENERIC_LITERAL_(Char, "0x"), 2);
					}
				}
				
				if(_rhs.isZero()){
					os.put(_BIG_NUM_GENERIC_LITERAL_(Char, '0'));
					return os;
				}
				
				static const Char *hexDigit = nullptr;
				if((os.flags() & std::ios_base::uppercase) != 0){
					hexDigit = _BIG_NUM_GENERIC_LITERAL_(Char, "0123456789ABCDEF");
				}
				else{
					hexDigit = _BIG_NUM_GENERIC_LITERAL_(Char, "0123456789abcdef");
				}
				
				using Producer = typename _type::CompareCond<LogSizeT, 4, ENTRY_SIZE, 
					_LargePower2Radix<SizeT, BigInt>, 
					_ExactDigitExtract<SizeT, BigInt>, 
					_SmallPower2Radix<SizeT, BigInt>>::type;
				Producer producer(std::move(_rhs), 4);
				for(SizeT digit = producer._start();true;digit = producer._next()){
					assert(digit < 16);
					os.put(hexDigit[digit]);
					if(!producer._hasNext()){
						break;
					}
				}
				
				return os;
			}
			
			if(_rhs.isZero()){
				os.put(_BIG_NUM_GENERIC_LITERAL_(Char, '0'));
				return os;
			}
			
			static const Char *decDigit = _BIG_NUM_GENERIC_LITERAL_(Char, "0123456789");
			_DecimalRadix<SizeT, BigInt> producer(std::move(_rhs));
			for(SizeT digit = producer._start();true;digit = producer._next()){
				assert(digit < 10);
				os.put(decDigit[digit]);
				if(!producer._hasNext()){
					break;
				}
			}
			
			if((os.flags() & std::ios_base::unitbuf) != 0){
				os.flush();
			}
			
			return os;
		}
		template <typename Char, class Trait>
		friend std::basic_ostream<Char, Trait> &operator<<(std::basic_ostream<Char, Trait> &os, const BigInt &_rhs){
			return os << BigInt(_rhs);
		}
		
		// coroutine simulator
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		RadixConvertEnumer<Integer> getDigitEnumer(Integer radix = 10) &&{
			return RadixConvertEnumer<Integer>(std::move(*this), radix);
		}
		template <typename Integer, 
			typename std::enable_if<std::is_integral<Integer>::value>::type * = nullptr>
		RadixConvertEnumer<Integer> getDigitEnumer(Integer radix = 10) const &{
			return RadixConvertEnumer<Integer>(*this, radix);
		}
		
#ifdef _BIG_NUM_DEBUG_
		template <class Ostream>
		void output(Ostream &out) const{
			out << "len=" << (int)buf.len << "\tcap=" << (int)buf.cap << "\tpositive=" << positive;
			out << "\tdata=\"";
			SizeT i(0);
			for(SizeT i = buf.len - 1;i > 0;--i){
				out << std::hex << (int)buf.data[i] << " ";
			}
			out << std::hex << buf.data[0];
			out << std::dec;
			out << "\"";
			return ;
		}
#endif // _BIG_NUM_DEBUG_
	private:
		explicit BigInt(std::nullptr_t)
			:buf(&allocator, nullptr), positive(false){}
		
		// case for unsigned integer
		template <typename Integer>
		void assignIntegral(Integer _rhs, std::false_type){
			if(Integer(0) == _rhs){
				buf.len = 1;
				buf.cap = 1;
				buf.data = allocator.allocate(1);
				allocator.construct(buf.data + 0, Ele(0));
				positive = true;
				return ;
			}
			
			using Common = typename std::common_type<Ele, Integer>::type;
			
			{
				Common tmp = _rhs;
				for(buf.len = 0;tmp > Common(0);++buf.len, tmp >>= ENTRY_SIZE);
				buf.setCap();
				buf.data = allocator.allocate(static_cast<std::size_t>(buf.cap));
			}
			{
				Common tmp = _rhs;
				SizeT i(0);
				try{
					for(;tmp > Common(0);++i, tmp >>= ENTRY_SIZE){
						allocator.construct(buf.data + i, static_cast<Ele>(tmp & ((1 << ENTRY_SIZE) - 1)));
					}
				}
				catch(...){
					destroyAll(buf.data, buf.data + i, allocator);
					allocator.deallocate(buf.data, buf.cap);
					buf.data = nullptr;
					buf.zeroLen();
					throw ;
				}
			}
			
			positive = true;
		}
		// case for signed integer
		template <typename Integer>
		void assignIntegral(Integer _rhs, std::true_type){
			if(Integer(0) == _rhs){
				buf.len = 1;
				buf.cap = 1;
				buf.data = allocator.allocate(1);
				allocator.construct(buf.data + 0, Ele(0));
				positive = true;
				return ;
			}
			
			if(Integer(1 << (sizeof(Integer) * CHAR_BIT - 1)) == _rhs){
				buf.setLen(static_cast<SizeT>((sizeof(Integer) * CHAR_BIT - 1 + ENTRY_SIZE - 1) / ENTRY_SIZE));
				buf.data = allocator.allocate(static_cast<std::size_t>(buf.cap));
				SizeT i(0);
				try{
					for(;i + 1 < buf.len;++i){
						allocator.construct(buf.data + i, Ele(0));
					}
					i = buf.len - 1;
					allocator.construct(buf.data + buf.len - 1, 1 << ((sizeof(Integer) + ENTRY_SIZE - 1) % ENTRY_SIZE));
				}
				catch(...){
					destroyAll(buf.data, buf.data + i, allocator);
					allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
					buf.data = nullptr;
					buf.zeroLen();
					throw ;
				}
				positive = false;
				return ;
			}
			
			using Common = typename std::common_type<Ele, Integer>::type;
			
			positive = _rhs > Integer(0);
			{
				Common tmp = positive? _rhs: -_rhs;
				for(buf.len = 0;tmp > Common(0);++buf.len, tmp >>= ENTRY_SIZE);
				buf.setCap();
				buf.data = allocator.allocate(static_cast<std::size_t>(buf.cap));
			}
			{
				Common tmp = positive? _rhs: -_rhs;
				SizeT i(0);
				try{
					for(;tmp > Common(0);++i, tmp >>= ENTRY_SIZE){
						allocator.construct(buf.data + i, static_cast<Ele>(tmp & ((1 << ENTRY_SIZE) - 1)));
					}
				}
				catch(...){
					destroyAll(buf.data, buf.data + i, allocator);
					allocator.deallocate(buf.data, buf.cap);
					buf.data = nullptr;
					buf.zeroLen();
					throw ;
				}
			}
		}
		
		void assignLv(const BigInt &_rhs, std::false_type){
			if(buf.data == nullptr){
				buf.data = allocator.allocate(static_cast<std::size_t>(_rhs.buf.cap));
			}
			else{
				destroyAll(buf.data, buf.data + buf.len, allocator);
				if(buf.cap != _rhs.buf.cap){
					allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
					
					buf.data = allocator.allocate(static_cast<std::size_t>(_rhs.buf.cap));
				}
			}
			
			SizeT i(0);
			try{
				for(;i < _rhs.buf.len;++i){
					allocator.construct(buf.data + i, _rhs.buf.data[i]);
				}
			}
			catch(...){
				destroyAll(buf.data, buf.data + i, allocator);
				allocator.deallocate(buf.data, static_cast<std::size_t>(_rhs.buf.cap));
				buf.data = nullptr;
				buf.zeroLen();
				throw ;
			}
			
			buf.len = _rhs.buf.len;
			buf.cap = _rhs.buf.cap;
			positive = _rhs.positive;
		}
		void assignLv(const BigInt &_rhs, std::true_type){
			if(buf.data == nullptr){
				allocator = _rhs.allocator;
				buf.data = allocator.allocate(static_cast<std::size_t>(_rhs.buf.cap));
			}
			else{
				destroyAll(buf.data, buf.data + buf.len, allocator);
				
				if((buf.cap != _rhs.buf.cap) || (allocator != _rhs.allocator)){
					allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
					
					allocator = _rhs.allocator;
					buf.data = allocator.allocate(static_cast<std::size_t>(buf.cap));
				}
				else{
					allocator = _rhs.allocator;
				}
			}
			
			SizeT i(0);
			try{
				for(;i < _rhs.buf.len;++i){
					allocator.construct(buf.data + i, _rhs.buf.data[i]);
				}
			}
			catch(...){
				destroyAll(buf.data, buf.data + i, allocator);
				allocator.deallocate(buf.data, static_cast<std::size_t>(_rhs.buf.cap));
				buf.data = nullptr;
				buf.zeroLen();
				throw ;
			}
			
			buf.len = _rhs.buf.len;
			buf.cap = _rhs.buf.cap;
			positive = _rhs.positive;
		}
		
		void assignRv(BigInt &&_rhs, std::true_type){
			if(nullptr != buf.data){
				destroyAll(buf.data, buf.data + buf.len, allocator);
				allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
			}
			
			allocator = std::move(_rhs.allocator);
			buf.len = std::move(_rhs.buf.len);
			buf.cap = std::move(_rhs.buf.cap);
			buf.data = std::move(_rhs.buf.data);
			positive = std::move(_rhs.positive);
			
			_rhs.buf.data = nullptr;
			_rhs.buf.zeroLen();
		}
		void assignRv(BigInt &&_rhs, std::false_type){
			if(allocator == _rhs.allocator){
				if(nullptr != buf.data){
					destroyAll(buf.data, buf.data + buf.len, allocator);
					allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
				}
				
				buf.len = std::move(_rhs.buf.len);
				buf.cap = std::move(_rhs.buf.cap);
				buf.data = std::move(_rhs.buf.data);
				positive = std::move(_rhs.positive);
				
				_rhs.buf.data = nullptr;
				_rhs.buf.zeroLen();
			}
			else{
				if(buf.data == nullptr){
					buf.data = allocator.allocate(static_cast<std::size_t>(_rhs.buf.cap));
				}
				else{
					destroyAll(buf.data, buf.data + buf.len, allocator);
					if(buf.cap != _rhs.buf.cap){
						allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
						
						buf.data = allocator.allocate(static_cast<std::size_t>(_rhs.buf.cap));
					}
				}
				
				SizeT i(0);
				try{
					for(;i < _rhs.buf.len;++i){
						allocator.construct(buf.data + i, std::move(_rhs.buf.data[i]));
					}
				}
				catch(...){
					destroyAll(buf.data, buf.data + i, allocator);
					allocator.deallocate(buf.data, static_cast<std::size_t>(_rhs.buf.cap));
					buf.data = nullptr;
					buf.zeroLen();
					throw ;
				}
				
				buf.len = std::move(_rhs.buf.len);
				buf.cap = std::move(_rhs.buf.cap);
				positive = std::move(_rhs.positive);
			}
		}
		
		inline void swapImpl(BigInt &_rhs, std::true_type){
			assert(buf.data != nullptr);
			assert(_rhs.buf.data != nullptr);
			
			using std::swap;
			
			swap(allocator, _rhs.allocator);
			swap(buf.len, _rhs.buf.len);
			swap(buf.cap, _rhs.buf.cap);
			swap(buf.data, _rhs.buf.data);
			swap(positive, _rhs.positive);
		}
		inline void swapImpl(BigInt &_rhs, std::false_type){
			assert(buf.data != nullptr);
			assert(_rhs.buf.data != nullptr);
			
			using std::swap;
			
			if(allocator == _rhs.allocator){
				swap(buf.len, _rhs.buf.len);
				swap(buf.cap, _rhs.buf.cap);
				swap(buf.data, _rhs.buf.data);
				swap(positive, _rhs.positive);
				
				return ;
			}
			else{
				// undefined behaviour, but we try to do something meaningful
				if(buf.len == _rhs.buf.len){
					assert(buf.cap == _rhs.buf.cap);
					
					for(SizeT i(0);i < buf.len;++i){
						swap(buf.data[i], _rhs.buf.data[i]);
					}
					
					return ;
				}
				
				DigitBuffer *l = nullptr;
				DigitBuffer *s = nullptr;
				if(buf.len < _rhs.buf.len){
					l = &_rhs.buf;
					s = &buf;
				}
				else{
					l = &buf;
					s = &_rhs.buf;
				}
				
				DigitBuffer stmp = s->realloc(l->len);
				if(nullptr == stmp){
					assert(l->cap == s->cap);
					// both buffers need not be reconstructed
					for(SizeT i(0);i < s->len;++i){
						swap(l->data[i], s->data[i]);
					}
					
					SizeT i(s->len);
					try{
						for(;i < l->len;++i){
							s->allocator->construct(s->data + i, std::move(l->data[i]));
						}
					}
					catch(...){
						for(SizeT j = s->len;j < i;++j){
							l->data[j] = std::move(s->data[j]);
						}
						destroyAll(s->data + s->len, s->data + i, *(s->allocator));
						throw ;
					}
					
					destroyAll(l->data + s->len, l->data + l->len, *(l->allocator));
					
					swap(buf.len, _rhs.buf.len);
					swap(positive, _rhs.positive);
				}
				else{
					// both need to be reconstructed
					SizeT i(0);
					try{
						for(;i < l->len;++i){
							stmp.allocator->construct(stmp.data + i, std::move(l->data[i]));
						}
					}
					catch(...){
						for(SizeT j(0);j < i;++j){
							l->data[i] = std::move(stmp.data[i]);
						}
						
						destroyAll(stmp.data, stmp.data + i, *(stmp.allocator));
						stmp.allocator->deallocate(stmp.data, static_cast<std::size_t>(stmp.cap));
						stmp.data = nullptr;
						stmp.zeroLen();
						
						throw ;
					}
					
					DigitBuffer ltmp = l->realloc(s->len);
					try{
						for(i = 0;i < s->len;++i){
							ltmp.allocator.construct(ltmp.data + i, std::move(s->data[i]));
						}
					}
					catch(...){
						for(SizeT j(0);j < s->len;j++){
							s->data[i] = std::move(ltmp.data[i]);
						}
						
						destroyAll(ltmp.data, ltmp + i, *(ltmp.allocator));
						ltmp.allocator->deallocate(ltmp.data, static_cast<std::size_t>(ltmp.cap));
						ltmp.data = nullptr;
						ltmp.zeroLen();
						
						throw ;
					}
					
					destroyAll(s->data, s->data + s->len, *(s->allocator));
					s->allocator->deallocate(s->data, static_cast<std::size_t>(s->cap));
					s->data = stmp.data;
					s->len = stmp.len;
					s->cap = stmp.cap;
					stmp.data = nullptr;
					stmp.zeroLen();
					
					destroyAll(l->data, l->data + l->len, *(l->allocator));
					l->allocator->deallocate(l->data, static_cast<std::size_t>(l->cap));
					l->data = ltmp.data;
					l->len = ltmp.len;
					l->cap = ltmp.cap;
					ltmp.data = nullptr;
					ltmp.zeroLen();
				}
				
				swap(positive, _rhs.positive);
				return ;
			}
		}
		
		inline std::int8_t compare(const BigInt &_rhs) const{
			return buf.compareRaw(_rhs.buf);
		}
		
		// compare with unsigned integer
		template <typename Integer>
		inline std::int8_t compareInt(const Integer &_rhs, std::false_type) const{
			if(!positive){
				return -1;
			}
			
			return buf.compareUnsignedIntBuffer(_rhs, std::integral_constant<SizeT, (sizeof(Integer) * CHAR_BIT + ENTRY_SIZE - 1) / ENTRY_SIZE>());
		}
		// compare with signed integer
		template <typename Integer>
		inline std::int8_t compareInt(const Integer &_rhs, std::true_type) const{
			if(_rhs < Integer(0)){
				if(positive){
					return 1;
				}
				// a simple hack. may contain bugs
				return -buf.compareUnsignedIntBuffer(static_cast<typename std::make_unsigned<Integer>::type>((~_rhs) + 1), std::integral_constant<SizeT, (sizeof(Integer) * CHAR_BIT + ENTRY_SIZE - 1) / ENTRY_SIZE>());
			}
			
			return buf.compareUnsignedIntBuffer(static_cast<typename std::make_unsigned<Integer>::type>(_rhs), std::integral_constant<SizeT, (sizeof(Integer) * CHAR_BIT + ENTRY_SIZE - 1) / ENTRY_SIZE>());
		}
		
		inline static BigInt truncateFrom(BigInt &&from, SizeT k){
			assert(k >= 0);
			if(0 == k){
				return std::move(from);
			}
			if(k >= from.buf.len){
				return BigInt();
			}
			
			BigInt tmp;
			destroyAll(tmp.buf.data, tmp.buf.data + tmp.buf.len, tmp.allocator);
			tmp.allocator.deallocate(tmp.buf.data, static_cast<std::size_t>(tmp.buf.cap));
			tmp.positive = from.positive;
			tmp.buf.setLen(from.buf.len - k);
			tmp.buf.data = tmp.allocator.allocate(static_cast<std::size_t>(tmp.buf.cap));
			SizeT i(0);
			try{
				for(;i < tmp.buf.len;++i){
					tmp.allocator.construct(tmp.buf.data + i, std::move(from.buf.data[i + k]));
				}
			}
			catch(...){
				for(SizeT j(0);j < i;++j){
					from.buf.data[i + k] = std::move(tmp.buf.data[i]);
				}
				
				destroyAll(tmp.buf.data, tmp.buf.data + i, tmp.allocator);
				tmp.allocator.deallocate(tmp.buf.data, static_cast<std::size_t>(tmp.buf.cap));
				tmp.buf.data = nullptr;
				tmp.buf.zeroLen();
				throw ;
			}
			
			return tmp;
		}
		inline static BigInt truncateFrom(const BigInt &from, SizeT k){
			assert(k >= 0);
			if(k == 0){
				return from;
			}
			if(k >= from.buf.len){
				return BigInt();
			}
			
			BigInt tmp;
			destroyAll(tmp.buf.data, tmp.buf.data + tmp.buf.len, tmp.allocator);
			tmp.allocator.deallocate(tmp.buf.data, static_cast<std::size_t>(tmp.buf.cap));
			tmp.positive = from.positive;
			tmp.buf.setLen(from.buf.len - k);
			tmp.buf.data = tmp.allocator.allocate(static_cast<std::size_t>(tmp.buf.cap));
			SizeT i(0);
			try{
				for(;i < tmp.buf.len;++i){
					tmp.allocator.construct(tmp.buf.data + i, from.buf.data[i + k]);
				}
			}
			catch(...){
				destroyAll(tmp.buf.data, tmp.buf.data + i, tmp.allocator);
				tmp.allocator.deallocate(tmp.buf.data, static_cast<std::size_t>(tmp.buf.cap));
				tmp.buf.data = nullptr;
				tmp.buf.zeroLen();
				throw ;
			}
			
			return tmp;
		}
		
		// propagate_on_XXXX_assignment?
		inline BigInt subStr(SizeT st, SizeT en) const &{
			if((0 == en) || (en > buf.len)){
				en = buf.len;
			}
			if(st >= en){
				return BigInt();
			}
			if((st >= buf.len) || (en > MAX_LEN)){
				throw std::out_of_range("BigInt::subStr");
			}
			
			BigInt tmp;
			// maybe we should add a private constructor construting null BigInt?
			destroyAll(tmp.buf.data, tmp.buf.data + tmp.buf.len, tmp.allocator);
			tmp.allocator.deallocate(tmp.buf.data, static_cast<std::size_t>(tmp.buf.cap));
			
			tmp.buf.setLen(en - st);
			tmp.buf.data = tmp.allocator.allocate(static_cast<std::size_t>(tmp.buf.cap));
			SizeT i(0);
			try{
				for(;i < tmp.buf.len;++i){
					tmp.allocator.construct(tmp.buf.data + i, buf.data[i + st]);
				}
			}
			catch(...){
				destroyAll(tmp.buf.data, tmp.buf.data + i, tmp.allocator);
				tmp.allocator.deallocate(tmp.buf.data, static_cast<std::size_t>(tmp.buf.cap));
				tmp.buf.data = nullptr;
				tmp.buf.zeroLen();
				throw ;
			}
			
			// a monkey hack
			tmp.positive = true;
			return tmp;
		}
		// TODO: consider the case when Ele is not trivally destructible
		inline BigInt subStr(SizeT st, SizeT en) &&{
			if((0 == en) || (en > buf.len)){
				en = buf.len;
			}
			if(st >= en){
				return BigInt();
			}
			if((st >= buf.len) || (en > MAX_LEN)){
				throw std::out_of_range("BigInt::subStr");
			}
			
			DigitBuffer tmp = buf.realloc(en - st);
			if(nullptr == tmp){
				// in-place
				for(SizeT i(0);(i + st) < en;++i){
					buf.data[i] = std::move(buf.data[i + st]);
				}
				
				destroyAll(buf.data + (en - st), buf.data + buf.len, allocator);
				allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
				buf.setLen(en - st);
				
				return std::move(*this);
			}
			else{
				// reconstruct
				SizeT i(0);
				try{
					for(;(i + st) < en;++i){
						allocator.construct(tmp.data + i, std::move(buf.data[i + st]));
					}
				}
				catch(...){
					destroyAll(tmp.data, tmp.data + i, allocator);
					allocator.deallocate(tmp.data, static_cast<std::size_t>(tmp.cap));
					tmp.data = nullptr;
					tmp.zeroLen();
					throw ;
				}
				
				destroyAll(buf.data, buf.data + buf.len, allocator);
				allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
				buf.data = std::move(tmp.data);
				buf.len = tmp.len;
				buf.cap = tmp.cap;
				
				return std::move(*this);
			}
		}
		
		// shl throws the exceeding bits away
		template <typename Integer>
		void shl(const Integer &_rhs, std::false_type){
			if(Integer(_rhs / ENTRY_SIZE) >= MAX_LEN){
				std::ostringstream out;
				out << "Cannot left shift a BigInt for more than " << ENTRY_SIZE * MAX_LEN << " bits";
				throw std::domain_error(out.str());
				// errno = EDOM;
			}
			
			SizeT padLen = static_cast<SizeT>(_rhs / ENTRY_SIZE);
			LogSizeT shLen = static_cast<LogSizeT>(_rhs % ENTRY_SIZE);
			Ele overflow = buf.data[buf.len - 1] >> (ENTRY_SIZE - shLen);
			
			SizeT _len = buf.len + padLen +((overflow > 0)? 1: 0);
			bool exceed = _len > MAX_LEN;
			if(exceed){
				_len = MAX_LEN;
			}
			
			DigitBuffer tmp = buf.realloc(_len);
			if(nullptr == tmp){
				// case for in-place reconstruction
				
				// exceeding
				assert(!exceed || (_len == MAX_LEN));
				
				if((overflow > 0) && !exceed){
					// overflow is considered eXpiring
					allocator.construct(buf.data + buf.len + padLen, std::move(overflow));
				}
				
				// note that since cap(len) == cap(_len), it holds that padLen == _len - len
				// <= cap(_len) - len == cap(len) - len < 2* len - len == len, even if the
				// result exceeds
				SizeT i = exceed?(_len - 1): (buf.len + padLen - 1);
				try{
					for(;i >= buf.len;--i){
						allocator.construct(buf.data + i, ((buf.data[i - padLen] << shLen) | (buf.data[i - padLen - 1] >> (ENTRY_SIZE - shLen))) & ((1 << ENTRY_SIZE) - 1));
					}
				}
				catch(...){
					destroyAll(buf.data + i - 1, buf.data + _len, allocator);
					throw ;
				}
				for(i = buf.len - 1;i > padLen;--i){
					buf.data[i] = ((buf.data[i - padLen] << shLen) | (buf.data[i - padLen - 1] >> (ENTRY_SIZE - shLen))) & ((1 << ENTRY_SIZE) - 1);
				}
				buf.data[padLen] = (buf.data[0] << shLen) & ((1 << ENTRY_SIZE) - 1);
				for(i = 0; i < padLen;++i){
					buf.data[i] = Ele(0);
				}
				
				buf.setLen(_len);
				return ;
			}
			else{
				// case for a realloc
				assert(buf.cap < tmp.cap);
				// exceeding
				assert(!exceed || (tmp.len == MAX_LEN));
				assert(!exceed || (tmp.cap == MAX_LEN));
				
				if((overflow > 0) && !exceed){
					// overflow is considered eXpiring
					allocator.construct(tmp.data + buf.len + padLen, std::move(overflow));
				}
				SizeT i = exceed? (tmp.len - 1): (buf.len + padLen - 1);
				try{
					for(;i > padLen;--i){
						allocator.construct(tmp.data + i, ((buf.data[i - padLen] << shLen) | (buf.data[i - padLen - 1] >> (ENTRY_SIZE - shLen))) & ((1 << ENTRY_SIZE) - 1));
					}
					i = padLen;
					allocator.construct(tmp.data + padLen, (buf.data[0] << shLen) & ((1 << ENTRY_SIZE) - 1));
					if(i > 0){
						for(--i;true;--i){
							allocator.construct(tmp.data + i, Ele(0));
							
							if(i == 0){
								break;
							}
						}
					}
				}
				catch(...){
					destroyAll(tmp.data + i + 1, tmp.data + tmp.len, allocator);
					allocator.deallocate(tmp.data, tmp.cap);
					tmp.data = nullptr;
					throw ;
				}
				
				destroyAll(buf.data, buf.data + buf.len, allocator);
				allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
				buf.data = std::move(tmp.data);
				buf.setLen(tmp.len);
				tmp.data = nullptr;
				tmp.zeroLen();
				
				return ;
			}
		}
		// orignally designed for integer literals that are deduced to signed integer types
		template <typename Integer>
		void shl(const Integer &_rhs, std::true_type){
			if(_rhs < Integer(0)){
				throw std::domain_error("Cannot left shift a BigInt for less than 0 bits");
				// errno = EDOM;
			}
			else{
				// quick hack
				shl(static_cast<typename std::make_unsigned<Integer>::type>(_rhs), std::false_type{});
				
				return ;
			}
		}
		
		// unsigned
		template <typename Integer>
		void shr(const Integer &_rhs, std::false_type){
			if(Integer(_rhs / ENTRY_SIZE) >= MAX_LEN){
				std::ostringstream out;
				out << "Cannot right shift a BigInt for more than " << ENTRY_SIZE * MAX_LEN << " bits";
				throw std::domain_error(out.str());
				// errno = EDOM;
			}
			
			SizeT padLen = _rhs / ENTRY_SIZE;
			LogSizeT shLen = _rhs % ENTRY_SIZE;
			Ele leftover = buf.data[buf.len - 1] >> shLen;
			if((padLen + ((leftover > 0)? 0: 1)) >= buf.len){
				zerolize();
				return ;
			}
			
			SizeT _len = buf.len - padLen - ((leftover > 0)? 0: 1);
			assert(buf.len > 0);
			DigitBuffer tmp = buf.realloc(_len);
			if(nullptr == tmp){
				if((_len < 2) && (leftover > 0)){
					// in-place reconstruction would only happen when original data
					// has a length of 1
					assert(1 == buf.cap);
					assert(1 == buf.len);
					
					buf.data[0] = std::move(leftover);
					return ;
				}
				else{
					assert(buf.len >= (padLen + 1));
					for(SizeT i(0);i < buf.len - padLen - 1;++i){
						buf.data[i] = ((buf.data[i + padLen] >> shLen) | (buf.data[i + padLen + 1] << (ENTRY_SIZE - shLen))) & ((1 << ENTRY_SIZE) - 1);
					}
					
					if(leftover > 0){
						buf.data[buf.len - padLen - 1] = std::move(leftover);
					}
				}
				
				destroyAll(buf.data + _len, buf.data + buf.len, allocator);
				buf.len = _len;
				return ;
			}
			else{
				assert(_len == tmp.len);
				
				if((tmp.len < 2) && (leftover > 0)){
					assert(1 == tmp.cap);
					
					allocator.construct(tmp.data + 0, std::move(leftover));
				}
				else{
					assert(buf.len >= (padLen + 1));
					SizeT i(0);
					try{
						for(;i < buf.len - padLen - 1;++i){
							allocator.construct(tmp.data + i, ((buf.data[i + padLen] >> shLen) | (buf.data[i + padLen + 1] << (ENTRY_SIZE - shLen))) & ((1 << ENTRY_SIZE) - 1));
						}
						
						assert(i == buf.len - padLen - 1);
						if(leftover > 0){
							allocator.construct(tmp.data + i, std::move(leftover));
						}
					}
					catch(...){
						destroyAll(tmp.data, tmp.data + i, allocator);
						allocator.deallocate(tmp.data, static_cast<std::size_t>(tmp.cap));
						tmp.data = nullptr;
						tmp.zeroLen();
						throw ;
					}
				}
				
				destroyAll(buf.data, buf.data + buf.len, allocator);
				allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
				
				buf.data = std::move(tmp.data);
				buf.len = tmp.len;
				buf.cap = tmp.cap;
				tmp.data = nullptr;
				tmp.zeroLen();
				return ;
			}
		}
		// signed
		template <typename Integer>
		void shr(const Integer &_rhs, std::true_type){
			if(_rhs < Integer(0)){
				throw std::domain_error("Cannot right shift a BigInt for less than 0 bits");
				// errno = EDOM;
			}
			else{
				// quick hack
				shr(static_cast<typename std::make_unsigned<Integer>::type>(_rhs), std::false_type{});
				
				return ;
			}
		}
		
		inline SizeT lenOfBinary() const{
			for(LogSizeT i = ENTRY_SIZE;i > 0;--i){
				if((buf.data[buf.len - 1] & (1 << (i - 1))) != 0){
					return (buf.len - 1) * ENTRY_SIZE + i;
				}
			}
			assert(false);
		}
		
		inline void add(const BigInt &_rhs){
			try{
				if(positive == _rhs.positive){
					buf.addRaw(_rhs.buf);
				}
				else{
					if(buf.compareRaw(_rhs.buf) >= 0){
						buf.subRaw(_rhs.buf);
					}
					else{
						// hack
						BigInt tmp = _rhs;
						tmp.buf.subRaw(std::move(buf));
						*this = std::move(tmp);
					}
				}
			}
			catch(const std::underflow_error &e){
				assert(false);
			}
		}
		inline void add(BigInt &&_rhs){
			try{
				if(positive == _rhs.positive){
					buf.addRaw(std::move(_rhs.buf));
				}
				else{
					if(buf.compareRaw(_rhs.buf) >= 0){
						buf.subRaw(std::move(_rhs.buf));
					}
					else{
						// hack
						BigInt tmp = std::move(_rhs);
						tmp.buf.subRaw(std::move(buf));
						*this = std::move(tmp);
					}
				}
			}
			catch(const std::underflow_error &e){
				assert(false);
			}
		}
		
		inline void sub(const BigInt &_rhs){
			try{
				if(positive != _rhs.positive){
					buf.addRaw(_rhs.buf);
				}
				else{
					if(buf.compareRaw(_rhs.buf) >= 0){
						buf.subRaw(_rhs.buf);
					}
					else{
						BigInt tmp = _rhs;
						tmp.buf.subRaw(std::move(buf));
						*this = std::move(tmp);
					}
				}
			}
			catch(const std::underflow_error &e){
				assert(false);
			}
		}
		inline void sub(BigInt &&_rhs){
			try{
				if(positive != _rhs.positive){
					buf.addRaw(std::move(_rhs.buf));
				}
				else{
					if(buf.compareRaw(_rhs.buf) >= 0){
						buf.subRaw(std::move(_rhs.buf));
					}
					else{
						BigInt tmp = std::move(_rhs);
						tmp.buf.subRaw(std::move(buf));
						*this = std::move(tmp);
					}
				}
			}
			catch(const std::underflow_error &e){
				assert(false);
			}
		}
		
		inline void zerolize(){
			if(nullptr != buf.data){
				destroyAll(buf.data, buf.data + buf.len, allocator);
				allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
			}
			
			buf.data = allocator.allocate(1);
			allocator.construct(buf.data + 0, Ele(0));
			buf.len = 1;
			buf.cap = 1;
			positive = true;
		}
		
		inline bool isZero() const{
			if(1 != buf.len){
				return false;
			}
			assert(1 == buf.cap);
			assert(nullptr != buf.data);
			if(Ele(0) != buf.data[0]){
				return false;
			}
#ifdef _BIG_NUM_DEBUG_
			if(!positive){
				int a = 1;
			}
#endif // _BIG_NUM_DEBUG_
			assert(positive);
			return true;
		}
		
		inline void changeSign(){
			if(isZero()){
				assert(positive);
				return ;
			}
			positive = !positive;
		}
		
		void selfMultiply(){
			positive = true;
			buf.resize(buf.cap << 1);
			assert(buf.len == buf.cap);
			
			ModularP_T root = pow(ModularP_T(OMEGA), PRI_ORDER / buf.len);
			fft1DPower2(static_cast<std::size_t>(buf.len), root, [this](std::size_t i){
				return buf.data[static_cast<SizeT>(i)];
			}, [this](std::size_t i) -> Ele &{
				return buf.data[static_cast<SizeT>(i)];
			});
			
			for(SizeT i = 0;i < buf.len;++i){
				buf.data[i] = Ele(ModularP_T(buf.data[i]) * buf.data[i]);
			}
			
			ModularP_T(invRoot) = pow(ModularP_T(OMEGA), PRI_ORDER / buf.len * (buf.len - 1));
			fft1DPower2(static_cast<std::size_t>(buf.len), invRoot, [this](std::size_t i){
				return buf.data[static_cast<SizeT>(i)];
			}, [this](std::size_t i) -> Ele &{
				return buf.data[static_cast<SizeT>(i)];
			});
			
			ModularP_T inverseN = pow(ModularP_T(TWO_INV), std::round(std::log2(buf.len)));
			for(SizeT i = 0;i < buf.len;++i){
				buf.data[i] = static_cast<Ele>(buf.data[i] * inverseN);
			}
			
			buf.propagateCarry();
			buf.shrinkToFit();
		}
		
		void multiplyMedium(const BigInt &_rhs){
			positive = (positive == _rhs.positive);
			
			DigitBuffer tmp(&allocator, nullptr);
			tmp.cap = static_cast<SizeT>(std::pow(2.0, std::ceil(std::log2(buf.len + _rhs.buf.len))));
			tmp.len = tmp.cap;
			tmp.data = allocator.allocate(static_cast<std::size_t>(tmp.cap));
			
#ifdef _BIG_NUM_DEBUG_
			SizeT _len1 = buf.len;
			SizeT _len2 = _rhs.buf.len;
#endif // _BIG_NUM_DEBUG_
			
			SizeT i(0);
			try{
				for(;i < _rhs.buf.len;++i){
					allocator.construct(tmp.data + i, _rhs.buf.data[i]);
				}
				for(i = _rhs.buf.len;i < tmp.len;++i){
					allocator.construct(tmp.data + i, Ele(0));
				}
			}
			catch(...){
				destroyAll(tmp.data, tmp.data + i, allocator);
				allocator.deallocate(tmp.data, static_cast<std::size_t>(tmp.cap));
				tmp.data = nullptr;
				tmp.zeroLen();
				throw ;
			}
			
			buf.resize(tmp.len);
			assert(buf.len == tmp.len);
			assert(tmp.len == tmp.cap);
			
			ModularP_T root = pow(ModularP_T(OMEGA), PRI_ORDER / buf.len);
			fft1DPower2(static_cast<std::size_t>(buf.len), root, [this](std::size_t i){
				return buf.data[static_cast<SizeT>(i)];
			}, [this](std::size_t i) -> Ele &{
				return buf.data[static_cast<SizeT>(i)];
			});
			fft1DPower2(static_cast<std::size_t>(tmp.len), root, [&tmp](std::size_t i){
				return tmp.data[static_cast<SizeT>(i)];
			}, [&tmp](std::size_t i) -> Ele &{
				return tmp.data[static_cast<SizeT>(i)];
			});
			
			for(SizeT i(0);i < buf.len;++i){
				buf.data[i] = Ele(ModularP_T(buf.data[i]) * ModularP_T(tmp.data[i]));
			}
			
			ModularP_T invRoot = pow(ModularP_T(OMEGA), PRI_ORDER / buf.len * (buf.len - 1));
			fft1DPower2(static_cast<std::size_t>(buf.len), invRoot, [this](std::size_t i){
				return buf.data[static_cast<SizeT>(i)];
			}, [this](std::size_t i) -> Ele &{
				return buf.data[static_cast<SizeT>(i)];
			});
			
			ModularP_T inverseN = pow(ModularP_T(TWO_INV), std::round(std::log2(buf.len)));
			//ModularP_T inverseN = ModularP_T(ModularP_T(static_cast<Ele>(buf.len))).inverse();
			for(SizeT i(0);i < buf.len;++i){
				buf.data[i] = Ele(buf.data[i] * inverseN);
			}
			
#ifdef _BIG_NUM_DEBUG_
			assert(((_len1 + _len2) == buf.cap) || (buf.data[_len1 + _len2]== Ele(0)));
#endif // _BIG_NUM_DEBUG_
			
			destroyAll(tmp.data, tmp.data + tmp.len, allocator);
			allocator.deallocate(tmp.data, static_cast<std::size_t>(tmp.cap));
			tmp.data = nullptr;
			tmp.zeroLen();
			
			buf.propagateCarry();
			buf.shrinkToFit();
		}
		void multiplyMedium(BigInt &&_rhs){
			positive = (positive == _rhs.positive);
			
			_rhs.buf.resize(static_cast<SizeT>(std::pow(2.0, std::ceil(std::log2(buf.len + _rhs.buf.len)))));
			assert(_rhs.buf.len == _rhs.buf.cap);
			buf.resize(_rhs.buf.len);
			assert(buf.len == _rhs.buf.len);
			
			ModularP_T root = pow(ModularP_T(OMEGA), PRI_ORDER / buf.len);
			fft1DPower2(static_cast<std::size_t>(buf.len), root, [this](std::size_t i){
				return buf.data[static_cast<SizeT>(i)];
			}, [this](std::size_t i) -> Ele &{
				return buf.data[static_cast<SizeT>(i)];
			});
			fft1DPower2(static_cast<std::size_t>(_rhs.buf.len), root, [&_rhs](std::size_t i){
				return _rhs.buf.data[static_cast<SizeT>(i)];
			}, [&_rhs](std::size_t i) -> Ele &{
				return _rhs.buf.data[static_cast<SizeT>(i)];
			});
			
			for(SizeT i(0);i < buf.len;++i){
				buf.data[i] = Ele(ModularP_T(buf.data[i]) * ModularP_T(_rhs.buf.data[i]));
			}
			
			ModularP_T invRoot = pow(ModularP_T(OMEGA), PRI_ORDER / buf.len * (buf.len - 1));
			fft1DPower2(static_cast<std::size_t>(buf.len), invRoot, [this](std::size_t i){
				return buf.data[static_cast<SizeT>(i)];
			}, [this](std::size_t i) -> Ele &{
				return buf.data[static_cast<SizeT>(i)];
			});
			
			ModularP_T inverseN = pow(ModularP_T(TWO_INV), std::round(std::log2(buf.len)));
			for(SizeT i(0);i < buf.len;++i){
				buf.data[i] = Ele(buf.data[i] * inverseN);
			}
			
			buf.propagateCarry();
			buf.shrinkToFit();
		}
		
		// originally designed for small value _rhs converted from integers, 
		// however, it might be used in other situations
		void multiplySmall(BigInt &&_rhs){
			positive = (positive == _rhs.positive);
			
			SizeT _len = buf.len + _rhs.buf.len - 1;
			if(_len > MAX_LEN){
				throw std::out_of_range("BigInt::mutiplySmall");
				// errno = ERANGE;
			}
			
			// caculate the best spiling size LogSizeT
			long double _tmin = std::numeric_limits<long double>::max();
			SizeT N;
			SizeT _N = (_rhs.buf.cap == _rhs.buf.len)? (_rhs.buf.cap << 1): _rhs.buf.cap;
			// smallest power of 2 greater than _rhs.len: _rhs.cap;
			for(;_N <= _len;_N <<= 1){
				// ceil(buf.len / (_N - _rhs.buf.len + 1))
				long double cost = (buf.len + _N - _rhs.buf.len) / (_N + 1 - _rhs.buf.len) * std::log2(_N + 1) * _N;
				if(_tmin > cost){
					_tmin = cost;
					N = _N;
				}
			}
			SizeT L = N - _rhs.buf.len + 1;
			//assert((buf.len + L - 1) / L > 1);
			
			SizeT _repe = buf.len / L;
			SizeT _left = buf.len % L;
			
			ModularP_T root = pow(ModularP_T(OMEGA), PRI_ORDER / N);
			ModularP_T invRoot = pow(ModularP_T(OMEGA), PRI_ORDER / N * (N - 1));
			ModularP_T inverseN = pow(ModularP_T(TWO_INV), std::round(std::log2(N)));
			
			_rhs.buf.resize(N);
			fft1DPower2(static_cast<std::size_t>(N), root, [&_rhs](std::size_t i){
				return _rhs.buf.data[static_cast<SizeT>(i)];
			}, [&_rhs](std::size_t i) -> Ele &{
				return _rhs.buf.data[static_cast<SizeT>(i)];
			});
			
			DigitBuffer tmp(&allocator, nullptr);
			tmp.len = N;
			tmp.cap = N;
			tmp.data = allocator.allocate(static_cast<std::size_t>(tmp.cap));
			SizeT i = 0;
			try{
				for(;i < tmp.len;++i){
					allocator.construct(tmp.data + i, 0);
				}
			}
			catch(...){
				destroyAll(tmp.data, tmp.data + i, allocator);
				allocator.deallocate(tmp.data, static_cast<std::size_t>(tmp.cap));
				tmp.data = nullptr;
				tmp.zeroLen();
				
				throw;
			}
			_BIG_NUM_ADD_FINALIZER_CODE_(
				destroyAll(tmp.data, tmp.data + tmp.len, allocator);
				allocator.deallocate(tmp.data, static_cast<std::size_t>(tmp.cap));
				tmp.data = nullptr;
				tmp.zeroLen();
			);
			
			// since _repe > 0
			// nL ... nL + L - 1
			fft1DPower2(static_cast<std::size_t>(N), root, [this, &L](std::size_t i){
				if(static_cast<SizeT>(i) < L){
					return buf.data[static_cast<SizeT>(i)];
				}
				else{
					return Ele(0);
				}
			}, [&tmp, this](std::size_t i) -> Ele&{
				return tmp.data[static_cast<SizeT>(i)];
			});
			
			for(SizeT i(0);i < N;++i){
				tmp.data[i] = Ele(ModularP_T(tmp.data[i]) * ModularP_T(_rhs.buf.data[i]));
			}
			fft1DPower2(static_cast<std::size_t>(tmp.len), invRoot, [&tmp](std::size_t i){
				return tmp.data[static_cast<SizeT>(i)];
			}, [&tmp](std::size_t i) -> Ele &{
				return tmp.data[static_cast<SizeT>(i)];
			});
			
			DigitBuffer _buf(&allocator, _len);
			i = 0;
			try{
				for(;i < N;++i){
					allocator.construct(_buf.data + i, Ele(tmp.data[i] * inverseN));
				}
			}
			catch(...){
				destroyAll(_buf.data, _buf.data + i, allocator);
				allocator.deallocate(_buf.data, static_cast<std::size_t>(_buf.cap));
				_buf.data = nullptr;
				_buf.zeroLen();
				
				throw ;
			}
			
			for(SizeT k(1);k < _repe;++k){
				fft1DPower2(static_cast<std::size_t>(N), root, [this, &L, &k](std::size_t i){
					if(static_cast<SizeT>(i) < L){
						return buf.data[k * L + static_cast<SizeT>(i)];
					}
					else{
						return Ele(0);
					}
				}, [&tmp](std::size_t i) -> Ele&{
					return tmp.data[static_cast<SizeT>(i)];
				});
				for(SizeT i(0);i < N;++i){
					tmp.data[i] = Ele(ModularP_T(tmp.data[i]) * ModularP_T(_rhs.buf.data[i]));
				}
				fft1DPower2(static_cast<std::size_t>(tmp.len), invRoot, [&tmp](std::size_t i){
					return tmp.data[static_cast<SizeT>(i)];
				}, [&tmp](std::size_t i) -> Ele &{
					return tmp.data[static_cast<SizeT>(i)];
				});
				
				for(SizeT i(0);i < N - L;++i){
					_buf.data[k * L + i] += Ele(tmp.data[i] * inverseN);
				}
				SizeT i = N - L;
				try{
					for(;i < N;++i){
						allocator.construct(_buf.data + k * L + i, Ele(tmp.data[i] * inverseN));
					}
				}
				catch(...){
					destroyAll(_buf.data, _buf.data + k * L + i, allocator);
					allocator.deallocate(_buf.data, static_cast<std::size_t>(_buf.cap));
					_buf.data = nullptr;
					_buf.zeroLen();
					
					throw ;
				}
			}
			
			if(_left > 0){
				fft1DPower2(static_cast<std::size_t>(N), root, [this, &_left, &_repe, &L](std::size_t i){
					if(static_cast<std::size_t>(i) < _left){
						return buf.data[_repe * L + static_cast<SizeT>(i)];
					}
					else{
						return Ele(0);
					}
				}, [&tmp](std::size_t i) -> Ele&{
					return tmp.data[static_cast<SizeT>(i)];
				});
				for(SizeT i(0);i < N;++i){
					tmp.data[i] = Ele(ModularP_T(tmp.data[i]) * ModularP_T(_rhs.buf.data[i]));
				}
				fft1DPower2(static_cast<std::size_t>(tmp.len), invRoot, [&tmp](std::size_t i){
					return tmp.data[static_cast<SizeT>(i)];
				}, [&tmp](std::size_t i) -> Ele &{
					return tmp.data[static_cast<SizeT>(i)];
				});
				for(SizeT i(0);i < N - L;++i){
					_buf.data[_repe * L + i] += Ele(tmp.data[i] * inverseN);
				}
				i = N - L;
				try{
					assert(_repe * L + N >= _len);
					// the rest are discarded since they are sure to be zero.
					for(;i < _len - _repe * L;++i){
						allocator.construct(_buf.data + _repe * L + i, Ele(tmp.data[i] * inverseN));
					}
				}
				catch(...){
					destroyAll(_buf.data, _buf.data + _repe * L + i, allocator);
					allocator.deallocate(_buf.data, static_cast<std::size_t>(_buf.cap));
					_buf.data = nullptr;
					_buf.zeroLen();
					
					throw ;
				}
			}
			
			_buf.propagateCarry();
			_buf.shrinkToFit();
			
			destroyAll(buf.data, buf.data + buf.len, allocator);
			allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
			buf.data = _buf.data;
			buf.len = _buf.len;
			buf.cap = _buf.cap;
			_buf.data = nullptr;
			_buf.zeroLen();
			
			return ;
		}
		void multiplySmall(const BigInt &_rhs){
			positive = (positive == _rhs.positive);
			
			SizeT _len = buf.len + _rhs.buf.len - 1;
			if(_len > MAX_LEN){
				throw std::out_of_range("BigInt::mutiplySmall");
				// errno = ERANGE;
			}
			
			long double _tmin = std::numeric_limits<long double>::max();
			SizeT N;
			SizeT _N = (_rhs.buf.cap == _rhs.buf.len)? (_rhs.buf.cap << 1): _rhs.buf.cap;
			for(;_N <= _len;_N <<= 1){
				long double cost = (buf.len + _N - _rhs.buf.len) / (_N + 1 - _rhs.buf.len) * std::log2(_N + 1) * _N;
				if(_tmin > cost){
					_tmin = cost;
					N = _N;
				}
			}
			SizeT L = N - _rhs.buf.len + 1;
			//assert((buf.len + L - 1) / L > 1);
			
			SizeT _repe = buf.len / L;
			SizeT _left = buf.len % L;
			
			ModularP_T root = pow(ModularP_T(OMEGA), PRI_ORDER / N);
			ModularP_T invRoot = pow(ModularP_T(OMEGA), PRI_ORDER / N * (N - 1));
			ModularP_T inverseN = pow(ModularP_T(TWO_INV), std::round(std::log2(N)));
			
			DigitBuffer rBuf(&allocator, nullptr);
			rBuf.cap = N;
			rBuf.len = N;
			rBuf.data = allocator.allocate(static_cast<std::size_t>(rBuf.cap));
			// TODO: change functors to iterators
			SizeT i = 0;
			try{
				for(;i < rBuf.len;++i){
					allocator.construct(rBuf.data + i, 0);
				}
			}
			catch(...){
				destroyAll(rBuf.data, rBuf.data + i, allocator);
				allocator.deallocate(rBuf.data, static_cast<std::size_t>(rBuf.cap));
				rBuf.data = nullptr;
				rBuf.zeroLen();
				
				throw;
			}
			_BIG_NUM_ADD_FINALIZER_CODE_(
				destroyAll(rBuf.data, rBuf.data + rBuf.len, allocator);
				allocator.deallocate(rBuf.data, static_cast<std::size_t>(rBuf.cap));
				rBuf.data = nullptr;
				rBuf.zeroLen();
			);
			
			fft1DPower2(static_cast<std::size_t>(N), root, [&_rhs](std::size_t i){
				if(i < _rhs.buf.len){
					return _rhs.buf.data[static_cast<SizeT>(i)];
				}else{
					return Ele(0);
				}
			}, [&rBuf, this](std::size_t i) -> Ele &{
				return rBuf.data[static_cast<SizeT>(i)];
			});
			
			DigitBuffer tmp(&allocator, nullptr);
			tmp.len = N;
			tmp.cap = N;
			tmp.data = allocator.allocate(static_cast<std::size_t>(tmp.cap));
			i = 0;
			try{
				for(;i < tmp.len;++i){
					allocator.construct(tmp.data + i, 0);
				}
			}
			catch(...){
				destroyAll(tmp.data, tmp.data + i, allocator);
				allocator.deallocate(tmp.data, static_cast<std::size_t>(tmp.cap));
				tmp.data = nullptr;
				tmp.zeroLen();
				
				throw;
			}
			_BIG_NUM_ADD_FINALIZER_CODE_(
				destroyAll(tmp.data, tmp.data + tmp.len, allocator);
				allocator.deallocate(tmp.data, static_cast<std::size_t>(tmp.cap));
				tmp.data = nullptr;
				tmp.zeroLen();
			);
			
			// since _repe > 0
			// nL ... nL + L - 1
			fft1DPower2(static_cast<std::size_t>(N), root, [this, &L](std::size_t i){
				if(static_cast<SizeT>(i) < L){
					return buf.data[static_cast<SizeT>(i)];
				}
				else{
					return Ele(0);
				}
			}, [&tmp, this](std::size_t i) -> Ele &{
				return tmp.data[static_cast<SizeT>(i)];
			});
			
			for(SizeT i(0);i < N;++i){
				tmp.data[i] = Ele(ModularP_T(tmp.data[i]) * ModularP_T(rBuf.data[i]));
			}
			fft1DPower2(static_cast<std::size_t>(tmp.len), invRoot, [&tmp](std::size_t i){
				return tmp.data[static_cast<SizeT>(i)];
			}, [&tmp](std::size_t i) -> Ele &{
				return tmp.data[static_cast<SizeT>(i)];
			});
			
			DigitBuffer _buf(&allocator, _len);
			i = 0;
			try{
				for(;i < N;++i){
					allocator.construct(_buf.data + i, Ele(tmp.data[i] * inverseN));
				}
			}
			catch(...){
				destroyAll(_buf.data, _buf.data + i, allocator);
				allocator.deallocate(_buf.data, static_cast<std::size_t>(_buf.cap));
				_buf.data = nullptr;
				_buf.zeroLen();
				
				throw ;
			}
			
			for(SizeT k(1);k < _repe;++k){
				fft1DPower2(static_cast<std::size_t>(N), root, [this, &L, &k](std::size_t i){
					if(static_cast<SizeT>(i) < L){
						return buf.data[k * L + static_cast<SizeT>(i)];
					}
					else{
						return Ele(0);
					}
				}, [&tmp](std::size_t i) -> Ele &{
					return tmp.data[static_cast<SizeT>(i)];
				});
				for(SizeT i(0);i < N;++i){
					tmp.data[i] = Ele(ModularP_T(tmp.data[i]) * ModularP_T(rBuf.data[i]));
				}
				fft1DPower2(static_cast<std::size_t>(tmp.len), invRoot, [&tmp](std::size_t i){
					return tmp.data[static_cast<SizeT>(i)];
				}, [&tmp](std::size_t i) -> Ele &{
					return tmp.data[static_cast<SizeT>(i)];
				});
				
				for(SizeT i(0);i < N - L;++i){
					_buf.data[k * L + i] += Ele(tmp.data[i] * inverseN);
				}
				SizeT i = N - L;
				try{
					for(;i < N;++i){
						allocator.construct(_buf.data + k * L + i, Ele(tmp.data[i] * inverseN));
					}
				}
				catch(...){
					destroyAll(_buf.data, _buf.data + k * L + i, allocator);
					allocator.deallocate(_buf.data, static_cast<std::size_t>(_buf.cap));
					_buf.data = nullptr;
					_buf.zeroLen();
					
					throw ;
				}
			}
			
			if(_left > 0){
				fft1DPower2(static_cast<std::size_t>(N), root, [this, &_left, &_repe, &L](std::size_t i){
					if(static_cast<std::size_t>(i) < _left){
						return buf.data[_repe * L + static_cast<SizeT>(i)];
					}
					else{
						return Ele(0);
					}
				}, [&tmp](std::size_t i) -> Ele &{
					return tmp.data[static_cast<SizeT>(i)];
				});
				for(SizeT i(0);i < N;++i){
					tmp.data[i] = Ele(ModularP_T(tmp.data[i]) * ModularP_T(rBuf.data[i]));
				}
				fft1DPower2(static_cast<std::size_t>(tmp.len), invRoot, [&tmp](std::size_t i){
					return tmp.data[static_cast<SizeT>(i)];
				}, [&tmp](std::size_t i) -> Ele &{
					return tmp.data[static_cast<SizeT>(i)];
				});
				for(SizeT i(0);i < N - L;++i){
					_buf.data[_repe * L + i] += Ele(tmp.data[i] * inverseN);
				}
				i = N - L;
				try{
					assert(_repe * L + N >= _len);
					for(;i < _len - _repe * L;++i){
						allocator.construct(_buf.data + _repe * L + i, Ele(tmp.data[i] * inverseN));
					}
				}
				catch(...){
					destroyAll(_buf.data, _buf.data + _repe * L + i, allocator);
					allocator.deallocate(_buf.data, static_cast<std::size_t>(_buf.cap));
					_buf.data = nullptr;
					_buf.zeroLen();
					
					throw ;
				}
			}
			_buf.propagateCarry();
			_buf.shrinkToFit();
			
			destroyAll(buf.data, buf.data + buf.len, allocator);
			allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
			buf.data = _buf.data;
			buf.len = _buf.len;
			buf.cap = _buf.cap;
			_buf.data = nullptr;
			_buf.zeroLen();
			
			return ;
		}
#ifdef _BIG_NUM_DEBUG_
		
		void trivalMultiply(const BigInt &_rhs){
			positive = (positive == _rhs.positive);
			SizeT _len = buf.len + _rhs.buf.len;
			if(_len > MAX_LEN){
				throw std::out_of_range("BigInt::trivalMultiply");
			}
			DigitBuffer tmp(&allocator, nullptr);
			tmp.setLen(_len);
			tmp.data = allocator.allocate(static_cast<std::size_t>(tmp.cap));
			SizeT i = 0, j = 0;
			try{
				for(;i < tmp.len;++i){
					allocator.construct(tmp.data + i, 0);
				}
			}
			catch(...){
				destroyAll(tmp.data, tmp.data + i, allocator);
				allocator.deallocate(tmp.data, static_cast<std::size_t>(tmp.cap));
				tmp.data = nullptr;
				tmp.zeroLen();
				throw ;
			}
			for(i = 0;i < buf.len;++i){
				for(j= 0;j < _rhs.buf.len;++j){
					tmp.data[i + j] += buf.data[i] * _rhs.buf.data[j];
				}
				tmp.propagateCarry();
			}
			
			tmp.propagateCarry();
			tmp.shrinkToFit();
			
			destroyAll(buf.data, buf.data + buf.len, allocator);
			allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
			buf.data = tmp.data;
			buf.len = tmp.len;
			buf.cap = tmp.cap;
			tmp.data = nullptr;
			tmp.zeroLen();
			return ;
		}
#endif // _BIG_NUM_DEBUG_
		
		void multiply(BigInt &&_rhs){
			constexpr auto _SMALL_MEDIUM_THRESHOLD_ = 2;
			
			if(isZero()){
				assert(positive);
				return ;
			}
			if(_rhs.isZero()){
				zerolize();
				return ;
			}
			
			if(buf.len >= _rhs.buf.len * _SMALL_MEDIUM_THRESHOLD_){
				multiplySmall(std::move(_rhs));
				return ;
			}
			if(_rhs.buf.len < buf.len * _SMALL_MEDIUM_THRESHOLD_){
				multiplyMedium(std::move(_rhs));
				return ;
			}
			assert(_rhs.buf.len >= buf.len * _SMALL_MEDIUM_THRESHOLD_);
			std::move(_rhs).multiplySmall(std::move(*this));
			*this = std::move(_rhs);
			return ;
			//multiplyMedium(std::move(_rhs));
			//trivalMultiply(_rhs);
		}
		void multiply(const BigInt &_rhs){
			constexpr auto _SMALL_MEDIUM_THRESHOLD_ = 2;
			
			if(isZero()){
				assert(positive);
				return ;
			}
			if(_rhs.isZero()){
				zerolize();
				return ;
			}
			
			if(buf.len >= _rhs.buf.len * _SMALL_MEDIUM_THRESHOLD_){
				multiplySmall(_rhs);
				return ;
			}
			if(_rhs.buf.len < buf.len * _SMALL_MEDIUM_THRESHOLD_){
				multiplyMedium(_rhs);
				return ;
			}
			assert(_rhs.buf.len >= buf.len * _SMALL_MEDIUM_THRESHOLD_);
			BigInt tmp = _rhs;
			tmp.multiplySmall(std::move(*this));
			*this = std::move(tmp);
			return ;
			//multiplyMedium(_rhs);
			//trivalMultiply(_rhs);
		}
		
		// case for unsigned int type small enough to be hold in a uint{ENTRY_SIZE}_t
		// so that the iterate integer multiply won't overflow since Ele is 
		// uleast{2*ENTRY_SIZE}_t
		template <typename UnsignedInt>
		inline void multiplyUnsignedInt(const UnsignedInt &_rhs, std::true_type){
			for(SizeT i(0);i < buf.len;++i){
				buf.data[i] *= _rhs;
			}
			
			buf.propagateCarry();
		}
		// case for large integer types
		template <typename UnsignedInt>
		inline void multiplyUnsignedInt(const UnsignedInt &_rhs, std::false_type){
			if(_rhs < static_cast<UnsignedInt>(1 << ENTRY_SIZE)){
				multiplyUnsignedInt<Ele>(static_cast<Ele>(_rhs), std::true_type{});
				return ;
			}
			else{
				multiplySmall(BigInt(_rhs));
				return ;
			}
		}
		
		// unsigned integer
		template <typename Integer>
		inline void multiplyInt(const Integer &_rhs, std::false_type){
			multiplyUnsignedInt(_rhs, std::integral_constant<bool, sizeof(Integer) <= ENTRY_SIZE>());
		}
		// signed integer
		template <typename Integer>
		inline void multiplyInt(const Integer &_rhs, std::true_type){
			if(_rhs < Integer(0)){
				positive = !positive;
				
				if(Integer(1 << (sizeof(Integer) * CHAR_BIT - 1)) == _rhs){
					shl(sizeof(Integer) * CHAR_BIT - 1, std::false_type());
					return ;
				}
				
				multiplyUnsignedInt(static_cast<typename std::make_unsigned<Integer>::type>(-_rhs), std::false_type());
				return ;
			}
			else{
				multiplyUnsignedInt(static_cast<typename std::make_unsigned<Integer>::type>(_rhs), std::false_type());
				return ;
			}
		}
		
		// this = (this * _rhs) / beta^{k}
		// TODO: optimization since the lowest k digits is sure to be truncated
		inline void multiplyShr(BigInt &&_rhs, SizeT k){
			bool _positive = positive == _rhs.positive;
			
			if(buf.len + _rhs.buf.len <= MAX_LEN){
				multiply(std::move(_rhs));
				shr(k, std::false_type{});
				positive = _positive;
				
				return ;
			}
			
			if(buf.len <= MAX_LEN / 2){
				assert(_rhs.buf.len > MAX_LEN / 2);
				BigInt r1 = _rhs.subStr(MAX_LEN >> 1, 0);
				BigInt r2 = std::move(_rhs).subStr(0, MAX_LEN >> 1);
				
				BigInt m1 = (*this) * std::move(r1);
				if(k < MAX_LEN / 2 * ENTRY_SIZE){
					m1.shl(MAX_LEN / 2 * ENTRY_SIZE - k, std::false_type{});
				}
				else{
					m1.shr(k - MAX_LEN / 2 * ENTRY_SIZE, std::false_type{});
				}
				
				BigInt m2 = std::move(*this) * std::move(r2);
				m2.shr(k, std::false_type{});
				
				*this = m1 + m2;
				positive = _positive;
				
				return ;
			}
			if(_rhs.buf.len <= MAX_LEN / 2){
				assert(buf.len > MAX_LEN / 2);
				BigInt l1 = subStr(MAX_LEN >> 1, 0);
				BigInt l2 = std::move(*this).subStr(0, MAX_LEN >> 1);
				
				BigInt m1 = std::move(l1) * _rhs;
				if(k < MAX_LEN / 2 * ENTRY_SIZE){
					m1.shl(MAX_LEN / 2 * ENTRY_SIZE - k, std::false_type{});
				}
				else{
					m1.shr(k - MAX_LEN / 2 * ENTRY_SIZE, std::false_type{});
				}
				
				BigInt m2 = std::move(l2) * std::move(_rhs);
				m2.shr(k, std::false_type{});
				
				*this = m1 + m2;
				positive = _positive;
				
				return ;
			}
			
			assert(buf.len > MAX_LEN / 2);
			assert(_rhs.buf.len > MAX_LEN / 2);
			
			BigInt l1 = subStr(MAX_LEN / 2, 0);
			BigInt r1 = _rhs.subStr(MAX_LEN / 2, 0);
			
			BigInt m1 = l1 * r1;
			if(k < MAX_LEN * ENTRY_SIZE){
				m1.shl(MAX_LEN * ENTRY_SIZE - k, std::false_type{});
			}
			else{
				m1.shr(k - MAX_LEN * ENTRY_SIZE, std::false_type{});
			}
			
			BigInt r2 = std::move(_rhs).subStr(0, MAX_LEN / 2);
			BigInt m2 = std::move(l1) * r2;
			BigInt l2 = std::move(*this).subStr(0, MAX_LEN / 2);
			BigInt m3 = l2 * std::move(r1);
			if(k < MAX_LEN / 2 * ENTRY_SIZE){
				m2.shl(MAX_LEN / 2 * ENTRY_SIZE - k, std::false_type{});
				m3.shl(MAX_LEN / 2 * ENTRY_SIZE - k, std::false_type{});
			}
			else{
				m2.shr(k - MAX_LEN / 2 * ENTRY_SIZE, std::false_type{});
				m3.shr(k - MAX_LEN / 2 * ENTRY_SIZE, std::false_type{});
			}
			
			BigInt m4 = std::move(l2) * std::move(r2);
			m4.shr(k, std::false_type{});
			
			*this = m1 + m2 + m3 + m4;
			positive = _positive;
			
			return ;
		}
		inline void multiplyShr(const BigInt &_rhs, SizeT k){
			bool _positive = positive == _rhs.positive;
			
			if(buf.len + _rhs.buf.len <= MAX_LEN){
				multiply(_rhs);
				shr(k, std::false_type{});
				positive = _positive;
				return ;
			}
			
			if(buf.len <= MAX_LEN / 2){
				assert(_rhs.buf.len > MAX_LEN / 2);
				BigInt r1 = _rhs.subStr(MAX_LEN >> 1, 0);
				BigInt r2 = _rhs.subStr(0, MAX_LEN >> 1);
				
				BigInt m1 = (*this) * std::move(r1);
				if(k < MAX_LEN / 2 * ENTRY_SIZE){
					m1.shl(MAX_LEN / 2 * ENTRY_SIZE - k, std::false_type{});
				}
				else{
					m1.shr(k - MAX_LEN / 2 * ENTRY_SIZE, std::false_type{});
				}
				
				BigInt m2 = std::move(*this) * std::move(r2);
				m2.shr(k, std::false_type{});
				
				*this = m1 + m2;
				positive = _positive;
				
				return ;
			}
			if(_rhs.buf.len <= MAX_LEN / 2){
				assert(buf.len > MAX_LEN / 2);
				BigInt l1 = subStr(MAX_LEN >> 1, 0);
				BigInt l2 = std::move(*this).subStr(0, MAX_LEN >> 1);
				
				BigInt m1 = std::move(l1) * _rhs;
				if(k < MAX_LEN / 2 * ENTRY_SIZE){
					m1.shl(MAX_LEN / 2 * ENTRY_SIZE - k, std::false_type{});
				}
				else{
					m1.shr(k - MAX_LEN / 2 * ENTRY_SIZE, std::false_type{});
				}
				
				BigInt m2 = std::move(l2) * _rhs;
				m2.shr(k, std::false_type{});
				
				*this = m1 + m2;
				positive = _positive;
				
				return ;
			}
			
			assert(buf.len > MAX_LEN / 2);
			assert(_rhs.buf.len > MAX_LEN / 2);
			
			BigInt l1 = subStr(MAX_LEN / 2, 0);
			BigInt r1 = _rhs.subStr(MAX_LEN / 2, 0);
			
			BigInt m1 = l1 * r1;
			if(k < MAX_LEN * ENTRY_SIZE){
				m1.shl(MAX_LEN * ENTRY_SIZE - k, std::false_type{});
			}
			else{
				m1.shr(k - MAX_LEN * ENTRY_SIZE, std::false_type{});
			}
			
			BigInt r2 = _rhs.subStr(0, MAX_LEN / 2);
			BigInt m2 = std::move(l1) * r2;
			BigInt l2 = std::move(*this).subStr(0, MAX_LEN / 2);
			BigInt m3 = l2 * std::move(r1);
			if(k < MAX_LEN / 2 * ENTRY_SIZE){
				m2.shl(MAX_LEN / 2 * ENTRY_SIZE - k, std::false_type{});
				m3.shl(MAX_LEN / 2 * ENTRY_SIZE - k, std::false_type{});
			}
			else{
				m2.shr(k - MAX_LEN / 2 * ENTRY_SIZE, std::false_type{});
				m3.shr(k - MAX_LEN / 2 * ENTRY_SIZE, std::false_type{});
			}
			
			BigInt m4 = std::move(l2) * std::move(r2);
			m4.shr(k, std::false_type{});
			
			*this = m1 + m2 + m3 + m4;
			positive = _positive;
			
			return ;
		}
		
		// TODO: full implemetation
		inline void multiplyTruncate(BigInt &&_rhs, SizeT k){
			multiplyShr(std::move(_rhs), k * ENTRY_SIZE);
		}
		inline void multiplyTruncate(const BigInt &_rhs, SizeT k){
			multiplyShr(_rhs, k * ENTRY_SIZE);
		}
		
		// floor(beta^{k} / this)
		// here we specify beta = 2 to converage more quickly
		inline BigInt newtonInverse(SizeT k) const{
			SizeT lenBin = lenOfBinary();
			if(k < lenBin){
				return BigInt();
			}
			if(k == lenBin){
				return BigInt(static_cast<Ele>(1));
			}
			
			// once x[n] = q while q = floor(beta^{k}) and beta^{k} = q * this + r,
			// x[n + 1] will be q + floor(q * r / beta^{k}) which equals to q if
			// use the iteration algorithm below
			SizeT i(0);
			BigInt x[2];
			// initial guess
			// x_{init} = 1.5d for 0.5 < d < 1
			x[1] = static_cast<Ele>(3);
			if((lenBin + 1) < k){
				x[1].shl(k - lenBin - 1, std::false_type{});
			}
			
			BigInt pow2 = 1;
			pow2.shl(k + 1, std::false_type{});
			do{
				x[i & 1] = x[1 - (i & 1)];
				x[i & 1].multiplyShr(pow2 - x[1 - (i & 1)] * (*this), k);
				if(x[i & 1] == x[1 - (i & 1)]){
					break;
				}
#ifdef _BIG_NUM_DEBUG_
				/*std::cout << "Newton Inverse " << i + 1 << " steps:";
				x[i & 1].output(std::cout);
				std::cout << std::endl;*/
#endif // _BIG_NUM_DEBUG_
				++i;
			}while(true);
			// x[] is not a name of an automatic object, thus does not meet the requirement
			// of copy-elision. Move construction for return value object using lvalue shall
			// meet the requirement of copy-elision first, otherwise copy construction is 
			// performed. n4140 12.5.32
			return std::move(x[i & 1]);
		}
		
		// assumes this and _rhs is positive
		// returns std::pair(Q, R)
		// for two non-negative integer a, b such that miu - a <
		// beta^{k} / _rhs < miu + b, Barret reduction iterates at most
		// max(a + 2, b + 1) times to return the corret result
		template <class BigIntRef2, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef2 &&>::value>::type * = nullptr>
		inline std::pair<BigInt, BigInt> barretReduction(const BigInt &_rhs, BigIntRef2 &&miu) &&{
			assert(_rhs.buf.len > 0);
			assert(buf.len >= _rhs.buf.len);
			assert(_rhs.buf.len * 2 >= buf.len);
#ifdef _BIG_NUM_DEBUG_
				/*std::cout << "*this:\t" << std::endl;
				output(std::cout);
				std::cout << std::endl << "_rhs:\t" << std::endl;
				_rhs.output(std::cout);
				std::cout << std::endl << "miu:\t" << std::endl;
				miu.output(std::cout);
				std::cout << std::endl << "_rhs * miu:\t" << std::endl;
				(_rhs * miu).output(std::cout);
				std::cout << std::endl << "_rhs * (miu + 1):\t" << std::endl;
				(_rhs * (miu + 1)).output(std::cout);
				BigInt tmp1 = *this;
				tmp1.shl(static_cast<SizeT>(buf.len * ENTRY_SIZE), std::false_type{});
				BigInt tmp2 = (*this) * (_rhs * miu);
				std::cout << std::endl << "this * _rhs * miu:\t" << std::endl;
				tmp2.output(std::cout);
				std::cout << std::endl << "this * beta ^ m:\t" << std::endl;
				tmp1.output(std::cout);
				std::cout << std::endl << "this * beta ^ m - this * _rhs * miu:\t" << std::endl;
				(tmp1 - tmp2).output(std::cout);
				BigInt tmp3 = static_cast<Ele>(1);
				tmp3.shl(static_cast<SizeT>(buf.len * ENTRY_SIZE), std::false_type{});
				std::cout << std::endl << "beta ^ m - _rhs * miu:\t" << std::endl;
				(tmp3 - _rhs * miu).output(std::cout);
				std::cout << std::endl << "this * (beta ^ m - _rhs * miu):\t" << std::endl;
				((*this) * (tmp3 - _rhs * miu)).output(std::cout);
				std::cout.flush();*/
#endif // _BIG_NUM_DEBUG_
			
			std::pair<BigInt, BigInt> res;
			res.first = truncateFrom(*this, _rhs.buf.len - 1);
			res.first.multiplyTruncate(std::forward<BigIntRef2>(miu), buf.len - _rhs.buf.len + 1);
#ifdef _BIG_NUM_DEBUG_
				/*std::cout << std::endl << "q * _rhs:\t" << std::endl;
				(_rhs * res.first).output(std::cout);
				std::cout.flush();*/
#endif // _BIG_NUM_DEBUG_
			res.second = std::move(*this) - _rhs * res.first;
#ifdef _BIG_NUM_DEBUG_
				/*std::cout << std::endl << "r:\t" << std::endl;
				res.second.output(std::cout);
				std::cout.flush();*/
#endif // _BIG_NUM_DEBUG_
			do{
#ifdef _BIG_NUM_DEBUG_
				/*res.first.output(std::cout);
				std::cout << std::endl;
				res.second.output(std::cout);
				std::cout << std::endl;
				std::cout << "-----------------------------" << std::endl;*/
#endif // _BIG_NUM_DEBUG_
				if(!res.second.positive){
					res.second += _rhs;
					res.first -= 1;
				}
				else if(!(res.second < _rhs)){
					res.second -= _rhs;
					res.first += 1;
				}
				else{
#ifdef _BIG_NUM_DEBUG_
					//std::cout << "end loop" << std::endl;
#endif // _BIG_NUM_DEBUG_
					break;
				}
			}while(true);
			
			return res;
		}
		template <class BigIntRef2, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef2 &&>::value>::type * = nullptr>
		inline std::pair<BigInt, BigInt> barretReduction(const BigInt &_rhs, BigIntRef2 &&miu) const &{
			assert(_rhs.buf.len > 0);
			assert(buf.len >= _rhs.buf.len);
			assert(_rhs.buf.len * 2 >= buf.len);
			
			std::pair<BigInt, BigInt> res;
			res.first = truncateFrom(*this, _rhs.buf.len - 1);
			res.first.multiplyTruncate(std::forward<BigIntRef2>(miu), buf.len - _rhs.buf.len + 1);
			res.second = *this - _rhs * res.first;
			do{
				if(!res.second.positive){
					res.second += _rhs;
					res.first -= 1;
				}
				else if(!(res.second < _rhs)){
					res.second -= _rhs;
					res.first += 1;
				}
				else{
					break;
				}
			}while(true);
			
			return res;
		}
		template <class BigIntRef2, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef2 &&>::value>::type * = nullptr>
		inline BigInt barretResident(const BigInt &_rhs, BigIntRef2 &&miu) &&{
			assert(_rhs.buf.len > 0);
			assert(buf.len >= _rhs.buf.len);
			assert(_rhs.buf.len * 2 >= buf.len);
			
			BigInt Q, res;
			Q = truncateFrom(*this, _rhs.buf.len - 1);
			Q.multiplyTruncate(std::forward<BigIntRef2>(miu), buf.len - _rhs.buf.len + 1);
			res = std::move(*this) - _rhs * Q;
			do{
				if(!res.positive){
					res += _rhs;
				}
				else if(!(res < _rhs)){
					res -= _rhs;
				}
				else{
					break;
				}
			}while(true);
			
			return res;
		}
		template <class BigIntRef2, 
			typename std::enable_if<isRLRef<BigInt, BigIntRef2 &&>::value>::type * = nullptr>
		inline BigInt barretResident(const BigInt &_rhs, BigIntRef2 &&miu) const &{
			assert(_rhs.buf.len > 0);
			assert(buf.len >= _rhs.buf.len);
			assert(_rhs.buf.len * 2 >= buf.len);
			
			BigInt Q, res;
			Q = truncateFrom(*this, _rhs.buf.len - 1);
			Q.multiplyTruncate(std::forward<BigIntRef2>(miu), buf.len - _rhs.buf.len + 1);
			res = *this - _rhs * Q;
			do{
				if(!res.positive){
					res += _rhs;
				}
				else if(!(res < _rhs)){
					res -= _rhs;
				}
				else{
					break;
				}
			}while(true);
			
			return res;
		}
		
		// assume this and _rhs are non-negative.
		inline std::pair<BigInt, BigInt> divideByMedium(const BigInt &_rhs) &&{
			assert(positive);
			assert(_rhs.positive);
			
			if(_rhs.isZero()){
				throw std::domain_error("divide by zero");
				// errno = ERANGE;
			}
			
			if(isZero()){
				return std::make_pair<BigInt, BigInt>(BigInt(0), BigInt(0));
			}
			
			assert(buf.len >= _rhs.buf.len);
			assert(buf.len <= _rhs.buf.len * 2);
			
			//SizeT lenBin = lenOfBinary();
			//return std::move(*this).barretReduction(_rhs, _rhs.newtonInverse(lenBin));
			return std::move(*this).barretReduction(_rhs, _rhs.newtonInverse(buf.len * ENTRY_SIZE));
		}
		inline std::pair<BigInt, BigInt> divideByMedium(const BigInt &_rhs) const &{
			assert(positive);
			assert(_rhs.positive);
			
			if(_rhs.isZero()){
				throw std::domain_error("divide by zero");
				// errno = ERANGE;
			}
			
			if(isZero()){
				return std::make_pair<BigInt, BigInt>(BigInt(0), BigInt(0));
			}
			
			assert(buf.len >= _rhs.buf.len);
			assert(buf.len <= _rhs.buf.len * 2);
			
			//SizeT lenBin = lenOfBinary();
			//return barretReduction(_rhs, _rhs.newtonInverse(lenBin));
			return barretReduction(_rhs, _rhs.newtonInverse(buf.len * ENTRY_SIZE));
		}
		
		inline std::pair<BigInt, BigInt> divideBy(const BigInt &_rhs) &&{
			assert(positive);
			assert(_rhs.positive);
			if(_rhs.isZero()){
				throw std::domain_error("divide by zero");
				// errno = ERANGE;
			}
			
			if(isZero()){
				return std::make_pair<BigInt, BigInt>(BigInt(0), BigInt(0));
			}
			if(buf.len < _rhs.buf.len){
				return std::make_pair<BigInt, BigInt>(BigInt(0), std::move(*this));
			}
			
			if(buf.len <= _rhs.buf.len * 2){
				try{
					return std::move(*this).divideByMedium(_rhs);
				}
				catch(std::domain_error &){
					assert(false);
				}
			}
			
			std::pair<BigInt, BigInt> res;
			// TODO: avoid unnecessary memory allocation
			BigInt miu = _rhs.newtonInverse(2 * _rhs.buf.len * ENTRY_SIZE);
			SizeT L = buf.len;
			SizeT finish = 0;
			
			{
				// I am looking forward to P1044R1 merged into offical C++ standard
				std::pair<BigInt, BigInt> qr = subStr(buf.len - 2 * _rhs.buf.len, buf.len).barretReduction(_rhs, miu);
				
				destroyAll(res.first.buf.data, res.first.buf.data + res.first.buf.len, res.first.allocator);
				res.first.allocator.deallocate(res.first.buf.data, static_cast<std::size_t>(res.first.buf.cap));
				res.first.buf.setLen(buf.len - 2 * _rhs.buf.len + qr.first.buf.len);
				res.first.buf.data = res.first.allocator.allocate(static_cast<std::size_t>(res.first.buf.cap));
				
				SizeT i = buf.len - 2 * _rhs.buf.len;
				try{
					for(;i < buf.len - 2 * _rhs.buf.len + qr.first.buf.len;++i){
						res.first.allocator.construct(res.first.buf.data + i, std::move(qr.first.buf.data[i - (buf.len - 2 * _rhs.buf.len)]));
					}
				}
				catch(...){
					destroyAll(res.first.buf.data + (buf.len - 2 * _rhs.buf.len), res.first.buf.data + i, res.first.allocator);
					res.first.allocator.deallocate(res.first.buf.data, static_cast<std::size_t>(res.first.buf.cap));
					res.first.buf.data = nullptr;
					res.first.buf.zeroLen();
					
					throw ;
				}
				finish = buf.len - 2 * _rhs.buf.len;
				
				if(qr.second.isZero()){
					/*for(i = buf.len - 2 * _rhs.buf.len;i < buf.len;++i){
						buf.data[i] = 0;
					}*/
					destroyAll(buf.data + (buf.len - 2 * _rhs.buf.len), buf.data + buf.len, allocator);
					L -= 2 * _rhs.buf.len;
				}
				else{
					for(i = 0;i < qr.second.buf.len;++i){
						buf.data[buf.len - 2 * _rhs.buf.len + i] = std::move(qr.second.buf.data[i]);
					}
					destroyAll(buf.data + (buf.len - 2 * _rhs.buf.len + qr.second.buf.len), buf.data + buf.len, allocator);
					/*for(i = buf.len - 2 * _rhs.buf.len + qr.second.buf.len;i < buf.len;++i){
						buf.data[i] = 0;
					}*/
					L -= 2 * _rhs.buf.len - qr.second.buf.len;
				}
			}
			
			_BIG_NUM_ADD_FINALIZER_CODE_(
				destroyAll(buf.data, buf.data + L, allocator);
				allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
				buf.data = nullptr;
				buf.zeroLen();
			);
			
			while(L >= 2 * _rhs.buf.len){
				std::pair<BigInt, BigInt> qr = subStr(L - 2 * _rhs.buf.len, L).barretReduction(_rhs, miu);
				
				SizeT i = L - 2 * _rhs.buf.len;
				assert(i + qr.first.buf.len <= finish);
				try{
					for(;i < L - 2 * _rhs.buf.len + qr.first.buf.len;++i){
						res.first.allocator.construct(res.first.buf.data + i, std::move(qr.first.buf.data[i - (L - 2 * _rhs.buf.len)]));
					}
					for(i = L - 2 * _rhs.buf.len + qr.first.buf.len;i < finish;++i){
						res.first.allocator.construct(res.first.buf.data + i, Ele(0));
					}
				}
				catch(...){
					destroyAll(res.first.buf.data + L - 2 * _rhs.buf.len, res.first.buf.data + i, res.first.allocator);
					destroyAll(res.first.buf.data + finish, res.first.buf.data + res.first.buf.len, res.first.allocator);
					res.first.allocator.deallocate(res.first.buf.data, static_cast<std::size_t>(res.first.buf.cap));
					res.first.buf.data = nullptr;
					res.first.buf.zeroLen();
					
					throw ;
				}
				finish = L - 2 * _rhs.buf.len;
				
				if(qr.second.isZero()){
					/*for(i = L - 2 * _rhs.buf.len;i < L;++i){
						buf.data[i] = 0;
					}*/
					destroyAll(buf.data + (L - 2 * _rhs.buf.len), buf.data + L, allocator);
					L -= 2 * _rhs.buf.len;
				}
				else{
					for(i = 0;i < qr.second.buf.len;++i){
						buf.data[L - 2 * _rhs.buf.len + i] = std::move(qr.second.buf.data[i]);
					}
					destroyAll(buf.data + (L - 2 * _rhs.buf.len + qr.second.buf.len), buf.data + L, allocator);
					/*for(i = L - 2 * _rhs.buf.len + qr.second.buf.len;i < L;++i){
						buf.data[i] = 0;
					}*/
					L -= 2 * _rhs.buf.len - qr.second.buf.len;
				}
			}
			
			assert(L < 2 * _rhs.buf.len);
			if(L >= _rhs.buf.len){
				BigInt q;
				BigInt tmp = subStr(0, L);
				SizeT tmpLenBin = tmp.lenOfBinary();
				std::tie(q, res.second) = std::move(tmp).barretReduction(_rhs, _rhs.newtonInverse(tmpLenBin));
				assert((q.buf.len <= finish) || ((0 == finish) && (q.isZero())));
				SizeT i = 0;
				try{
					for(;i < q.buf.len;++i){
						res.first.allocator.construct(res.first.buf.data + i, std::move(q.buf.data[i]));
					}
					for(i = q.buf.len;i < finish;++i){
						res.first.allocator.construct(res.first.buf.data + i, 0);
					}
				}
				catch(...){
					destroyAll(res.first.buf.data, res.first.buf.data + i, res.first.allocator);
					destroyAll(res.first.buf.data + finish, res.first.buf.data + res.first.buf.len, res.first.allocator);
					res.first.allocator.deallocate(res.first.buf.data, static_cast<std::size_t>(res.first.buf.cap));
					res.first.buf.data = nullptr;
					res.first.buf.zeroLen();
					
					throw ;
				}
			}
			else{
				destroyAll(res.second.buf.data, res.second.buf.data + res.second.buf.len, res.second.allocator);
				res.second.allocator.deallocate(res.second.buf.data, static_cast<std::size_t>(res.second.buf.cap));
				res.second.buf.setLen(L);
				res.second.buf.data = res.second.allocator.allocate(static_cast<std::size_t>(res.second.buf.cap));
				SizeT i(0);
				try{
					for(;i < L;++i){
						res.second.allocator.construct(res.second.buf.data + i, std::move(buf.data[i]));
					}
				}
				catch(...){
					destroyAll(res.second.buf.data, res.second.buf.data + i, res.second.allocator);
					res.second.allocator.deallocate(res.second.buf.data, static_cast<std::size_t>(res.second.buf.cap));
					res.second.buf.data = nullptr;
					res.second.buf.zeroLen();
					
					throw ;
				}
				
				try{
					for(i = 0;i < finish;++i){
						res.first.allocator.construct(res.first.buf.data + i, 0);
					}
				}
				catch(...){
					destroyAll(res.first.buf.data, res.first.buf.data + i, res.first.allocator);
					destroyAll(res.first.buf.data + finish, res.first.buf.data + res.first.buf.len, res.first.allocator);
					res.first.buf.data = nullptr;
					res.first.buf.zeroLen();
					
					throw ;
				}
			}
			
			return res;
		}
		inline std::pair<BigInt, BigInt> divideBy(const BigInt &_rhs) const &{
			assert(positive);
			assert(_rhs.positive);
			
			return BigInt(*this).divideBy(_rhs);
		}
		
		inline BigInt modularByMedium(const BigInt &_rhs) &&{
			assert(positive);
			assert(_rhs.positive);
			if(_rhs.isZero()){
				throw std::domain_error("divide by zero");
				// errno = ERANGE;
			}
			if(isZero()){
				return BigInt();
			}
			assert(buf.len >= _rhs.buf.len);
			assert(buf.len <= _rhs.buf.len * 2);
			//SizeT lenBin = lenOfBinary();
			//return std::move(*this).barretResident(_rhs, _rhs.newtonInverse(lenBin));
			return std::move(*this).barretResident(_rhs, _rhs.newtonInverse(buf.len * ENTRY_SIZE));
		}
		inline BigInt modularByMedium(const BigInt &_rhs) const &{
			assert(positive);
			assert(_rhs.positive);
			if(_rhs.isZero()){
				throw std::domain_error("divide by zero");
				// errno = ERANGE;
			}
			if(isZero()){
				return BigInt();
			}
			assert(buf.len >= _rhs.buf.len);
			assert(buf.len <= _rhs.buf.len * 2);
			//SizeT lenBin = lenOfBinary();
			//return barretResident(_rhs, _rhs.newtonInverse(lenBin));
			return barretResident(_rhs, _rhs.newtonInverse(buf.len * ENTRY_SIZE));
		}
		
		inline BigInt modularBy(const BigInt &_rhs) &&{
			assert(positive);
			assert(_rhs.positive);
			if(_rhs.isZero()){
				throw std::domain_error("divide by zero");
				// errno = ERANGE;
			}
			
			if(isZero()){
				return BigInt();
			}
			if(buf.len < _rhs.buf.len){
				return std::move(*this);
			}
			
			if(buf.len <= _rhs.buf.len * 2){
				try{
					return std::move(*this).modularByMedium(_rhs);
				}
				catch(std::domain_error &){
					assert(false);
				}
			}
			
			// TODO: avoid unnecessary memory allocation
			BigInt miu = _rhs.newtonInverse(2 * _rhs.buf.len * ENTRY_SIZE);
			SizeT L = buf.len;
			while(L >= 2 * _rhs.buf.len){
				BigInt R = subStr(L - 2 * _rhs.buf.len, L).barretResident(_rhs, miu);
				
				if(R.isZero()){
					destroyAll(buf.data + (L - 2 * _rhs.buf.len), buf.data + L, allocator);
					L -= 2 * _rhs.buf.len;
				}
				else{
					for(SizeT i(0);i < R.buf.len;++i){
						buf.data[L - 2 * _rhs.buf.len + i] = std::move(R.buf.data[i]);
					}
					destroyAll(buf.data + (L - 2 * _rhs.buf.len + R.buf.len), buf.data + L, allocator);
					L -= 2 * _rhs.buf.len - R.buf.len;
				}
			}
			
			_BIG_NUM_ADD_FINALIZER_CODE_(
				destroyAll(buf.data, buf.data + L, allocator);
				allocator.deallocate(buf.data, static_cast<std::size_t>(buf.cap));
				buf.data = nullptr;
				buf.zeroLen();
			);
			
			assert(L < 2 * _rhs.buf.len);
			if(L >= _rhs.buf.len){
				BigInt tmp = subStr(0, L);
				SizeT tmpLenBin = tmp.lenOfBinary();
				BigInt resident = std::move(tmp).barretResident(_rhs, _rhs.newtonInverse(tmpLenBin));
				
				return resident;
			}
			else{
				BigInt resident;
				
				destroyAll(resident.buf.data, resident.buf.data + resident.buf.len, resident.allocator);
				resident.allocator.deallocate(resident.buf.data, static_cast<std::size_t>(resident.buf.cap));
				resident.buf.setLen(L);
				resident.buf.data = resident.allocator.allocate(static_cast<std::size_t>(resident.buf.cap));
				SizeT i(0);
				try{
					for(;i < L;++i){
						resident.allocator.construct(resident.buf.data + i, std::move(buf.data[i]));
					}
				}
				catch(...){
					destroyAll(resident.buf.data, resident.buf.data + i, resident.allocator);
					resident.allocator.deallocate(resident.buf.data, static_cast<std::size_t>(resident.buf.cap));
					resident.buf.data = nullptr;
					resident.zeroLen();
					
					throw ;
				}
				
				return resident;
			}
		}
		inline BigInt modularBy(const BigInt &_rhs) const &{
			assert(positive);
			assert(_rhs.positive);
			
			return BigInt(*this).modularBy(_rhs);
		}
		
		template <typename UnsignedInt>
		inline UnsignedInt convertSingleDigit() const{
			UnsignedInt digit(0);
			for(SizeT i = buf.len - 1;true;--i){
				assert((digit < (digit << ENTRY_SIZE) + static_cast<UnsignedInt>(buf.data[i]))
					|| (isZero()));
				digit = (digit << ENTRY_SIZE) + static_cast<UnsignedInt>(buf.data[i]);
				if(i == 0){
					break;
				}
			}
			return digit;
		}
		
		inline static std::vector<BigInt> decimalBaseInit(){
			std::vector<BigInt> res;
			res.emplace_back(static_cast<Ele>(10));
			return res;
		}
		
		inline static std::vector<BigInt> &getDecimalBase(){
			static std::vector<BigInt> decimalBase = decimalBaseInit();
			return decimalBase;
		}
		
		Alloc allocator;
		DigitBuffer buf;
		bool positive; // true - 0 or positive, false - negative
	}; // class BigInt
	
	// static members passed as an argument for a reference parameter is odr-used, thus
	// have to be defined outside class scope. For constexpr static members declared in
	// a template class, defination can be like those below:
	//template <class Allocator>
	//constexpr typename BigInt<Allocator>::Ele BigInt<Allocator>::OMEGA;
	//template <class Allocator>
	//constexpr typename BigInt<Allocator>::Ele BigInt<Allocator>::TWO_INV;
	
	template <class Alloc>
	inline void swap(BigInt<Alloc> &_lhs, BigInt<Alloc> &_rhs){
		_lhs.swap(_rhs);
		return ;
	}
	
	using bigint_t = BigInt<>;
	
	template <char ...Args>
	bigint_t operator ""_bigint(){
		using List = _type::StaticList<std::integral_constant<char, Args>...>;
		using Parser = _type::LiteralParser<bigint_t, List>;
		return Parser::getBI();
	}
	
#define _BIG_INT_MAKE_LITERAL_OPERATOR_(biginttype, suffix) \
namespace bignum{ \
	\
	template <char ...Args> \
	biginttype operator "" suffix(){ \
		using List = _type::StaticList<std::integral_constant<char, Args>...>; \
		using Parser = _type::LiteralParser<biginttype, List>; \
		return Parser::getBI(); \
	} \
	\
};
	
}; // namespace bignum
#endif // _BIG_INT_HPP_