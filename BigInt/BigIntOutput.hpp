#ifndef _BIG_NUM_HPP_
#error "This header must be included through BigNum.hpp"
#endif // _BIG_NUM_HPP_

#ifndef _BIG_INT_OUTPUT_HPP_
#define _BIG_INT_OUTPUT_HPP_

#include <stack>
#include <vector>
#include <tuple>
#include <type_traits>
#include <utility>
#include <iterator>
#include <stdexcept>
#include <cmath>
#include <memory>

#include "../Libs/BigNumTypeTrait.hpp"

namespace bignum{
	
	namespace _type{
		
		template <typename Digit, class BI>
		class DigitProducer{
		public:
			//virtual DigitProducer *clone() const = 0;
			virtual Digit _start() = 0;
			virtual Digit _next() = 0;
			virtual bool _hasNext() const = 0;
			// destructors must have a function body to reference even though
			// they are pure virtual
			virtual ~DigitProducer(){}
		};
		
		template <typename Digit, class BI, class Derived>
		class DigitProducerCRTP:public DigitProducer<Digit, BI>{
		/*public:
			virtual DigitProducer<Digit, BI> *clone() const{
				return new Derived(static_cast<const Derived &>(*this));
			}*/
		};
		
		template <typename Digit, class BI>
		class DigitEnumIterator
			:public std::iterator<
				std::input_iterator_tag, 	// iterator_category
				Digit, 						// value_type
				std::ptrdiff_t, 			// difference_type
				const Digit *, 				// pointer
				const Digit &>{				// reference
		private:
			using Producer = DigitProducer<Digit, BI>;
		public:
			using iterator_category = std::input_iterator_tag;
			using value_type = Digit;
			using difference_type = std::ptrdiff_t;
			using pointer = const Digit *;
			using reference = const Digit &;
		private:
			class ValueProxy{
			public:
				explicit ValueProxy(Digit &_value): value(_value){}
				
				ValueProxy(const ValueProxy &) = delete;
				ValueProxy(ValueProxy &&) = default;
				
				ValueProxy &operator=(const ValueProxy &) = delete;
				ValueProxy &operator=(ValueProxy &&) = delete;
				
				~ValueProxy() = delete;
				
				value_type operator*() const{
					return value;
				}
			private:
				Digit value;
			};
		public:
			DigitEnumIterator(const DigitEnumIterator &_rhs) = default;
			DigitEnumIterator(DigitEnumIterator &&_rhs) noexcept
				:producer(std::move(_rhs.producer)), digit(_rhs.digit){
				_rhs.producer = nullptr;
			}
			
			// we must ensure all the _producer arguments in all instances contain 
			// different addresses or we are f**ked.
			template <class Derived>
			explicit DigitEnumIterator(Derived *_producer)
				:producer(_producer), digit(producer->_start()){}
			
			explicit DigitEnumIterator(std::nullptr_t)
				:producer(nullptr), digit(){}
			
			DigitEnumIterator &operator=(const DigitEnumIterator &_rhs) = default;
			DigitEnumIterator &operator=(DigitEnumIterator &&_rhs) = default;
			
			~DigitEnumIterator() = default;
			
			void swap(DigitEnumIterator &_rhs){
				using std::swap;
				swap(producer, _rhs.producer);
				swap(digit, _rhs.digit);
			}
			
			// only a proxy, since it takes some RTTI support to truly determine equality
			// for subtyping
			friend bool operator==(const DigitEnumIterator &i, const DigitEnumIterator &j){
				return (i.producer == j.producer) && (i.digit == j.digit);
			}
			friend bool operator!=(const DigitEnumIterator &i, const DigitEnumIterator &j){
				return (i.producer != j.producer) || (i.digit != j.digit);
			}
			
			reference operator*() const{
				return digit;
			}
			
			// pre-increment
			DigitEnumIterator &operator++(){
				if(!producer->_hasNext()){
					digit = Digit();
					producer = nullptr;
				}
				else{
					digit = producer->_next();
				}
				return *this;
			}
			// post-increment
			ValueProxy operator++(int){
				ValueProxy tmp(digit);
				if(!producer->_hasNext()){
					digit = Digit();
					producer = nullptr;
				}
				else{
					digit = producer->_next();
				}
				return tmp;
			}
		private:
			std::shared_ptr<Producer> producer;
			Digit digit;
		};
		
		template <typename Digit, class BI>
		inline void swap(DigitEnumIterator<Digit, BI> &_lhs, DigitEnumIterator<Digit, BI> &_rhs){
			_lhs.swap(_rhs);
		}
		
	}; // namespace _type
	
	// assume coroutine parameters are non-negative
	template <typename, class, typename = void>
	class _GenericRadix;
	template <typename, class, typename = void>
	class _DecimalRadix;
	template <typename, class, typename = void>
	class _SmallPower2Radix;
	template <typename, class, typename = void>
	class _LargePower2Radix;
	template <typename, class>
	class _ExactDigitExtract;
	
	using _type::DigitProducerCRTP;
	using _type::isRLRef;
	using _type::isSigned;
	using _type::DigitProducer;
	
	// unsigned Digit
	template <typename Digit, class BI>
	class _GenericRadix<Digit, BI, 
		typename std::enable_if<!isSigned<Digit>::value>::type>
		:public DigitProducerCRTP<Digit, BI, _GenericRadix<Digit, BI>>{
	private:
		using VSize = typename std::vector<BI>::size_type;
		using SizeT = typename BI::SizeT;
	public:
		_GenericRadix(const _GenericRadix &) = default;
		_GenericRadix(_GenericRadix &&) = delete;
		
		_GenericRadix &operator=(const _GenericRadix &) = delete;
		_GenericRadix &operator=(_GenericRadix &&) = delete;
		
		~_GenericRadix() = default;
		
		// coroutine arguments
		template <class BIRef, 
			typename std::enable_if<isRLRef<BI, BIRef &&>::value>::type * = nullptr>
		explicit _GenericRadix(BIRef &&_lhs, Digit _radix)
			:fd(0), md(_radix), rd(0), ed(_radix), rf(0), rr(0), radix(_radix){
			assert(_lhs.positive);
			
			if(_lhs.isZero()){
				fd = 0;
				rf = 1;
				md = radix;
				rd = 0;
				rr = 0;
				ed = radix;
				return ;
			}
			if(_lhs.compareInt(Digit(radix), std::false_type()) == -1){
				fd = _lhs.template convertSingleDigit<Digit>();
				rf = 1;
				md = radix;
				rd = 0;
				rr = 0;
				ed = radix;
				return ;
			}
			
			base.emplace_back(radix);
			do{
				BI next = base.back();
				next.selfMultiply();
				std::int8_t comp = next.compare(_lhs);
				if(-1 == comp){
					base.push_back(std::move(next));
					continue;
				}
				if(0 == comp){
					rf = 0;
					fd = 0;
					md = 1;
					rr = SizeT(1 << base.size());
					rd = 0;
					ed = radix;
					return ;
				}
				if(1 == comp){
					break;
				}
			}while(true);
			
			dStack.emplace(std::forward<BIRef &&>(_lhs), base.size(), 0);
		}
		
		virtual Digit _start(){
			return _next();
		}
		
		virtual Digit _next(){
			assert(_hasNext());
			
			{
				Digit res = trivalDigits();
				if(res != radix){
					return res;
				}
			}
			
			while(true){
				assert(!dStack.empty());
				// isZero() and <radix will only happen when top() is pushed
				// as a reminder of some division rather than quatiot, since
				// quatiot will not be pushed into stack as shown in the end
				// of the loop content
				if(std::get<0>(dStack.top()).isZero()){
					// this means q * radix ^ {t} + 0, that is, print t zeroes when 
					// facing this reminder, where t = 2 ^ {k}
					assert(std::get<2>(dStack.top()) == static_cast<SizeT>(1 << std::get<1>(dStack.top())));
					ed = radix;
					rd = 0;
					rr = 0;
					md = radix;
					fd = 0;
					rf = std::get<2>(dStack.top());
					dStack.pop();
					return trivalDigits();
				}
				if(std::get<0>(dStack.top()).compareInt(Digit(radix), std::false_type()) == -1){
					Digit tmp = std::move(std::get<0>(dStack.top())).template convertSingleDigit<Digit>();
					rr = 0;
					rd = 0;
					ed = radix;
					if(std::get<2>(dStack.top()) == 0){
						// highest digit
						md = radix;
						rf = 1;
						fd = tmp;
					}
					else{
						md = tmp;
						rf = std::get<2>(dStack.top()) - 1;
						fd = 0;
					}
					dStack.pop();
					return trivalDigits();
				}
				
				// linear search fits well with the divide tree depth
				// binary search does not
				VSize k = std::get<1>(dStack.top());
				std::int8_t comp;
				do{
					--k;
					comp = base[k].compare(std::get<0>(dStack.top()));
					if(1 == comp){
						continue;
					}
					else{
						break;
					}
					assert(k != 0);
				}while(true);
				
				if(0 == comp){
					ed = radix;
					md = 1;
					rd = 0;
					rr = static_cast<SizeT>(1 << k);
					fd = 0;
					rf = std::get<2>(dStack.top());
					if(rf != 0){
						rf -= 1 + rr;
					}
					dStack.pop();
					return trivalDigits();
				}
				
				assert(-1 == comp);
				std::pair<BI, BI> qr = std::move(std::get<0>(dStack.top())).divideByMedium(base[k]);
				SizeT rLen = SizeT(1 << k);
				SizeT qLen = std::get<2>(dStack.top());
				if(qLen != 0){
					assert(qLen > rLen);
					qLen -= rLen;
				}
				dStack.pop();
				
				// same as the procession at the beginning of the loop content. we add this just to
				// avoid some unnecessary move assignments
				if(qr.first.compareInt(Digit(radix), std::false_type{}) == -1){
					Digit tmp = std::move(qr.first).template convertSingleDigit<Digit>();
					assert(tmp != 0);
					if(qLen == 0){
						rf = 0;
					}
					else{
						rf = qLen - 1;
					}
					
					if(qr.second.isZero()){
						ed = radix;
						rd = 0;
						rr = rLen;
						md = tmp;
						fd = 0;
						return trivalDigits();
					}
					
					if(qr.second.compareInt(Digit(radix), std::false_type()) == -1){
						ed = std::move(qr.second).template convertSingleDigit<Digit>();
						rd = 0;
						rr = rLen - 1;
						md = tmp;
						fd = 0;
						return trivalDigits();
					}
					else{
						dStack.emplace(std::move(qr.second), k, rLen);
						ed = radix;
						rd = 0;
						rr = 0;
						md = tmp;
						fd = 0;
						return trivalDigits();
					}
				}
				else{
					dStack.emplace(std::move(qr.second), k, rLen);
					dStack.emplace(std::move(qr.first), k, qLen);
					continue;
				}
			}// while(true)
			// should never be reached
			assert(false);
			return 0;
		}
		
		virtual bool _hasNext() const{
			if(rf > 0){
				assert((fd >= 0) && (fd < radix));
				return true;
			}
			if((md >= 0) && (md < radix)){
				return true;
			}
			if(rr > 0){
				assert((rd >= 0) && (rd < radix));
				return true;
			}
			if((ed >= 0) && (ed < radix)){
				return true;
			}
			return !dStack.empty();
		}
	private:
		inline Digit trivalDigits(){
			if(rf > 0){
				--rf;
				return fd;
			}
			if((md >= 0) && (md < radix)){
				Digit res = md;
				md = radix;
				return res;
			}
			if(rr > 0){
				--rr;
				return rd;
			}
			if((ed >= 0) && (ed < radix)){
				Digit res = ed;
				ed = radix;
				return res;
			}
			return radix;
		}
		
		// coroutine states
		std::vector<BI> base;
		// BI: current BigInt to producer
		// VSize: BI > base[VSize] > sqrt(BI)
		// SizeT: the output string shall be SizeT long, or 0 if no leading zero needed
		std::stack<std::tuple<BI, VSize, SizeT>> dStack;
		// {fd}^(rf)+{md}+{rd}^(rf)+{ed}
		Digit fd, md, rd, ed;
		SizeT rf, rr;
		
		// coroutine arguments
		Digit radix;
	}; // class _GenericRadix<Digit, BI, void>
	// signed Digit
	template <typename Digit, class BI>
	class _GenericRadix<Digit, BI, 
		typename std::enable_if<isSigned<Digit>::value>::type>
		:public DigitProducerCRTP<Digit, BI, _GenericRadix<Digit, BI>>{
	private:
		using Unsigned = typename std::make_unsigned<Digit>::type;
		using base = _GenericRadix<Unsigned, BI>;
	public:
		_GenericRadix(const _GenericRadix &) = default;
		_GenericRadix(_GenericRadix &&) = delete;
		
		_GenericRadix &operator=(const _GenericRadix &) = delete;
		_GenericRadix &operator=(_GenericRadix &&) = delete;
		
		~_GenericRadix() = default;
		
		template <class BIRef, 
			typename std::enable_if<isRLRef<BI, BIRef &&>::value>::type * = nullptr>
		explicit _GenericRadix(BIRef &&_lhs, Digit radix)
			:producer(std::forward<BIRef>(_lhs), Unsigned(radix)){}
		
		virtual Digit _start(){
			Unsigned tmp = producer._start();
			return tmp;
		}
		
		virtual Digit _next(){
			Unsigned tmp = producer._next();
			return tmp;
		}
		
		virtual bool _hasNext() const{
			return producer._hasNext();
		}
	private:
		base producer;
	};
	
	// unsigned digit
	template <typename Digit, class BI>
	class _DecimalRadix<Digit, BI, 
		typename std::enable_if<!isSigned<Digit>::value>::type>
		:public DigitProducerCRTP<Digit, BI, _DecimalRadix<Digit, BI>>{
	private:
		using VSize = typename std::vector<BI>::size_type;
		using SizeT = typename BI::SizeT;
		
		constexpr static Digit radix = 10;
	public:
		_DecimalRadix(const _DecimalRadix &) = default;
		_DecimalRadix(_DecimalRadix &&) = delete;
		
		_DecimalRadix &operator=(const _DecimalRadix &) = delete;
		_DecimalRadix &operator=(_DecimalRadix &) = delete;
		
		~_DecimalRadix() = default;
		
		template <class BIRef, 
			typename std::enable_if<isRLRef<BI, BIRef &&>::value>::type * = nullptr>
		explicit _DecimalRadix(BIRef &&_lhs)
			:decimalBase(BI::getDecimalBase()), fd(0), md(radix), rd(0), ed(radix), rf(0), rr(0){
			assert(_lhs.positive);
			
			if(_lhs.isZero()){
				fd = 0;
				rf = 1;
				md = radix;
				rd = 0;
				rr = 0;
				ed = radix;
				return ;
			}
			
			if(_lhs.buf.len == 1){
				if(_lhs.buf.data[0] < radix){
					fd = Digit(_lhs.buf.data[0]);
					rf = 1;
					md = radix;
					rd = 0;
					rr = 0;
					ed = radix;
					return ;
				}
			}
			
			for(VSize k(0);k < decimalBase.size();++k){
				std::int8_t comp = decimalBase[k].compare(_lhs);
				if(0 == comp){
					fd = 0;
					rf = 0;
					md = 1;
					rd = 0;
					rr = SizeT(1 << k);
					ed = radix;
					return ;
				}
				if(1 == comp){
					dStack.emplace(std::forward<BIRef>(_lhs), k, 0);
					return ;
				}
			}
			
			if(dStack.empty()){
				do{
					BI next = decimalBase.back();
					next.selfMultiply();
					std::int8_t comp = next.compare(_lhs);
					if(-1 == comp){
						decimalBase.push_back(std::move(next));
						continue;
					}
					if(0 == comp){
						fd = 0;
						rf = 0;
						md = 1;
						rf = 0;
						rr = SizeT(1 << decimalBase.size());
						ed = radix;
						return ;
					}
					if(1 == comp){
						dStack.emplace(std::forward<BIRef>(_lhs), decimalBase.size(), 0);
						break;
					}
				}while(true);
			}
			return ;
		}
		
		virtual Digit _start(){
			return _next();
		}
		
		virtual Digit _next(){
			assert(_hasNext());
			
			{
				Digit res = trivalDigits();
				if(res != radix){
					return res;
				}
			}
			
			while(true){
				assert(!dStack.empty());
				if(std::get<0>(dStack.top()).isZero()){
					assert(std::get<2>(dStack.top()) == static_cast<SizeT>(1 << std::get<1>(dStack.top())));
					ed = radix;
					rd = 0;
					rr = 0;
					md = radix;
					fd = 0;
					rf = std::get<2>(dStack.top());
					dStack.pop();
					return trivalDigits();
				}
				if((std::get<0>(dStack.top()).buf.len == 1) && (std::get<0>(dStack.top()).buf.data[0] < radix)){
					Digit tmp = std::get<0>(dStack.top()).buf.data[0];
					rr = 0;
					rd = 0;
					ed = radix;
					if(std::get<2>(dStack.top()) == 0){
						md = radix;
						rf = 1;
						fd = tmp;
					}
					else{
						md = tmp;
						rf = std::get<2>(dStack.top()) - 1;
						fd = 0;
					}
					dStack.pop();
					return trivalDigits();
				}
				
				VSize k = std::get<1>(dStack.top());
				std::int8_t comp;
				do{
					--k;
					comp = decimalBase[k].compare(std::get<0>(dStack.top()));
					if(1 == comp){
						continue;
					}
					else{
						break;
					}
					assert(k != 0);
				}while(true);
				
				if(0 == comp){
					ed = radix;
					md = 1;
					rd = 0;
					rr = static_cast<SizeT>(1 << k);
					fd = 0;
					rf = std::get<2>(dStack.top());
					if(rf != 0){
						rf -= 1 + rr;
					}
					dStack.pop();
					return trivalDigits();
				}
				
				assert(-1 == comp);
				std::pair<BI, BI> qr = std::move(std::get<0>(dStack.top())).divideByMedium(decimalBase[k]);
				SizeT rLen = SizeT(1 << k);
				SizeT qLen = std::get<2>(dStack.top());
				if(qLen != 0){
					assert(qLen > rLen);
					qLen -= rLen;
				}
				dStack.pop();
				
				if((qr.first.buf.len == 1) && (qr.first.buf.data[0] < radix)){
					assert(!qr.first.isZero());
					Digit tmp = qr.first.buf.data[0];
					if(qLen == 0){
						rf = 0;
					}
					else{
						rf = qLen - 1;
					}
					
					if(qr.second.isZero()){
						ed = radix;
						rd = 0;
						rr = rLen;
						md = tmp;
						fd = 0;
						return trivalDigits();
					}
					if((qr.second.buf.len == 1) && (qr.second.buf.data[0] < radix)){
						ed = qr.second.buf.data[0];
						rd = 0;
						rr = rLen - 1;
						md = tmp;
						fd = 0;
						return trivalDigits();
					}
					else{
						dStack.emplace(std::move(qr.second), k, rLen);
						ed = radix;
						rd = 0;
						rr = 0;
						md = tmp;
						fd = 0;
						return trivalDigits();
					}
				}
				else{
					dStack.emplace(std::move(qr.second), k, rLen);
					dStack.emplace(std::move(qr.first), k, qLen);					
					continue;
				}
			}// while(true)
			// should never be reached
			assert(false);
			return 0;
		}
		
		virtual bool _hasNext() const{
			if(rf > 0){
				assert((fd >= 0) && (fd < radix));
				return true;
			}
			if((md >= 0) && (md < radix)){
				return true;
			}
			if(rr > 0){
				assert((rd >= 0) && (rd < radix));
				return true;
			}
			if((ed >= 0) && (ed < radix)){
				return true;
			}
			return !dStack.empty();
		}
	private:
		inline Digit trivalDigits(){
			if(rf > 0){
				--rf;
				return fd;
			}
			if((md >= 0) && (md < radix)){
				Digit res = md;
				md = radix;
				return res;
			}
			if(rr > 0){
				--rr;
				return rd;
			}
			if((ed >= 0) && (ed < radix)){
				Digit res = ed;
				ed = radix;
				return res;
			}
			return radix;
		}
		
		// coroutine state
		// TODO: make write operations on decimalBase thread-safe
		std::vector<BI> &decimalBase;
		std::stack<std::tuple<BI, VSize, SizeT>> dStack;
		Digit fd, md, rd, ed;
		SizeT rf, rr;
	};// class _DecimalRadix<Digit, BI, void>
	// signed digit
	template <typename Digit, class BI>
	class _DecimalRadix<Digit, BI, 
		typename std::enable_if<isSigned<Digit>::value>::type>
		:public DigitProducerCRTP<Digit, BI, _DecimalRadix<Digit, BI>>{
	private:
		using Unsigned = typename std::make_unsigned<Digit>::type;
		using base = _DecimalRadix<Unsigned, BI>;
	public:
		_DecimalRadix(const _DecimalRadix &) = default;
		_DecimalRadix(_DecimalRadix &&) = delete;
		
		_DecimalRadix &operator=(const _DecimalRadix &) = delete;
		_DecimalRadix &operator=(_DecimalRadix &&) = delete;
		
		~_DecimalRadix() = default;
		
		template <class BIRef, 
			typename std::enable_if<isRLRef<BI, BIRef &&>::value>::type * = nullptr>
		explicit _DecimalRadix(BIRef &&_lhs)
			:producer(std::forward<BIRef>(_lhs)){}
		
		virtual Digit _start(){
			Unsigned tmp = producer._start();
			return tmp;
		}
		
		virtual Digit _next(){
			Unsigned tmp = producer._next();
			return tmp;
		}
		
		virtual bool _hasNext() const{
			return producer._hasNext();
		}
	private:
		base producer;
	};
	
	// unsigned digit
	template <typename Digit, class BI>
	class _SmallPower2Radix<Digit, BI, 
		typename std::enable_if<!isSigned<Digit>::value>::type>
		:public DigitProducerCRTP<Digit, BI, _SmallPower2Radix<Digit, BI>>{
	private:
		using SizeT = typename BI::SizeT;
		using Ele = typename BI::Ele;
		using LogSizeT = typename BI::LogSizeT;
		
		constexpr static Ele ENTRY_SIZE = BI::ENTRY_SIZE;
	public:
		_SmallPower2Radix(const _SmallPower2Radix &) = default;
		_SmallPower2Radix(_SmallPower2Radix &&) = delete;
		
		_SmallPower2Radix &operator=(const _SmallPower2Radix &) = delete;
		_SmallPower2Radix &operator=(_SmallPower2Radix &&) = delete;
		
		~_SmallPower2Radix() = default;
		
		explicit _SmallPower2Radix(const BI &_lhs, SizeT _exp)
			:lhs(_lhs), exp(_exp), curI(_lhs.buf.len), curBit(0){
			assert(exp < ENTRY_SIZE);
		}
		
		virtual Digit _start(){
			if(lhs.isZero()){
				curI = 0;
				curBit = 0;
				return static_cast<Digit>(0);
			}
			
			LogSizeT highest = (lhs.buf.len * ENTRY_SIZE) % exp;
			Ele highestLeft = lhs.buf.data[lhs.buf.len - 1] >> (ENTRY_SIZE - highest);
			curI = lhs.buf.len - 1;
			curBit = ENTRY_SIZE - highest;
			if(highestLeft > 0){
				return static_cast<Digit>(highestLeft);
			}
			do{
				Ele digit = 0;
				if(curBit >= exp){
					digit = (lhs.buf.data[curI] >> (curBit - exp)) & ((1 << exp) - 1);
					curBit -= exp;
				}
				else{
					assert(curI > 0);
					digit = ((lhs.buf.data[curI] & ((1 << curBit) - 1)) << (exp - curBit)) | (lhs.buf.data[curI - 1] >> (curBit + ENTRY_SIZE - exp));
					curBit += ENTRY_SIZE - exp;
					--curI;
					assert(digit > 0);
				}
				if(digit > 0){
					return static_cast<Digit>(digit);
				}
			}while(true);
		}
		
		virtual Digit _next(){
			assert(_hasNext());
			
			Ele digit = 0;
			if(curBit >= exp){
				digit = (lhs.buf.data[curI] >> (curBit - exp)) & ((1 << exp) - 1);
				curBit -= exp;
			}
			else{
				assert(curI > 0);
				digit = ((lhs.buf.data[curI] & ((1 << curBit) - 1)) << (exp - curBit)) | (lhs.buf.data[curI - 1] >> (curBit + ENTRY_SIZE - exp));
				curBit += ENTRY_SIZE - exp;
				--curI;
			}
			return static_cast<Digit>(digit);
		}
		
		virtual bool _hasNext() const{
			return (curBit > 0) || (curI > 0);
		}
	private:
		// coroutine arguments
		// lhs dangling after its parent RadixConvertEnumer objects invalidates
		// pretty much like what pointer behaves after its referencing memory
		// invalidates
		const BI &lhs;
		SizeT exp;
		
		// coroutine states
		LogSizeT curBit;
		SizeT curI;
	};// class _SmallPower2Radix<Digit, BI, void>
	// signed digit
	template <typename Digit, class BI>
	class _SmallPower2Radix<Digit, BI, 
		typename std::enable_if<isSigned<Digit>::value>::type>
		:public DigitProducerCRTP<Digit, BI, _SmallPower2Radix<Digit, BI>>{
	private:
		using Unsigned = typename std::make_unsigned<Digit>::type;
		using base = _SmallPower2Radix<Unsigned, BI>;
		
		using SizeT = typename BI::SizeT;
	public:
		_SmallPower2Radix(const _SmallPower2Radix &) = default;
		_SmallPower2Radix(_SmallPower2Radix &&) = delete;
		
		_SmallPower2Radix &operator=(const _SmallPower2Radix &) = delete;
		_SmallPower2Radix &operator=(_SmallPower2Radix &&) = delete;
		
		~_SmallPower2Radix() = default;
		
		explicit _SmallPower2Radix(const BI &_lhs, SizeT _exp)
			:producer(_lhs, _exp){}
		
		virtual Digit _start(){
			Unsigned tmp = producer._start();
			return tmp;
		}
		
		virtual Digit _next(){
			Unsigned tmp = producer._next();
			return tmp;
		}
		
		virtual bool _hasNext() const{
			return producer._hasNext();
		}
	private:
		base producer;
	};
	
	// unsigned digit
	template <typename Digit, class BI>
	class _LargePower2Radix<Digit, BI, 
		typename std::enable_if<!isSigned<Digit>::value>::type>
		:public DigitProducerCRTP<Digit, BI, _LargePower2Radix<Digit, BI>>{
	private:
		using SizeT = typename BI::SizeT;
		using Ele = typename BI::Ele;
		using LogSizeT = typename BI::LogSizeT;
		
		constexpr static Ele ENTRY_SIZE = BI::ENTRY_SIZE;
	public:
		_LargePower2Radix(const _LargePower2Radix &) = default;
		_LargePower2Radix(_LargePower2Radix &&) = delete;
		
		_LargePower2Radix &operator=(const _LargePower2Radix &) = delete;
		_LargePower2Radix &operator=(_LargePower2Radix &&) = delete;
		
		~_LargePower2Radix() = default;
		
		explicit _LargePower2Radix(const BI &_lhs, SizeT _exp)
			:lhs(_lhs), exp(_exp), curI(_lhs.buf.len), curBit(0){
			assert(exp > ENTRY_SIZE);
		}
		
		virtual Digit _start(){
			curI = lhs.buf.len - 1;
			curBit = ENTRY_SIZE;
			
			SizeT leftBit = (lhs.buf.len * ENTRY_SIZE) % exp;
			SizeT leftI = leftBit / ENTRY_SIZE;
			LogSizeT leftLeft = leftBit % ENTRY_SIZE;
			
			Digit highest = 0;
			for(SizeT i = 0;i < leftI;++i){
				highest = (highest << ENTRY_SIZE) | static_cast<Digit>(lhs.buf.data[curI - i]);
			}
			highest = (highest << leftLeft) | static_cast<Digit>(lhs.buf.data[curI - leftI] >> (ENTRY_SIZE - leftLeft));
			assert((leftI == 0) || (highest > 0));
			
			assert((leftLeft == 0) || (curI >= leftI));
			assert((leftLeft > 0) || (curI >= leftI - 1));
			if(leftLeft > 0){
				curI -= leftI;
			}
			else{
				curI -= leftI - 1;
			}
			curBit -= leftLeft;
			
			if(highest > 0){
				return highest;
			}
			else{
				return _next();
			}
		}
		
		virtual Digit _next(){
			assert(_hasNext());
			
			SizeT stepI = (exp - curBit) / ENTRY_SIZE;
			LogSizeT leftStep = (exp - curBit) % ENTRY_SIZE;
			Digit digit = lhs.buf.data[curI] & ((Digit(1) << curBit) - 1);
			for(SizeT i(1);i <= stepI;++i){
				digit = (digit << ENTRY_SIZE) | static_cast<Digit>(lhs.buf.data[curI - i]);
			}
			if(leftStep > 0){
				assert(curI > stepI);
				digit = (digit << leftStep) | (static_cast<Digit>(lhs.buf.data[curI - stepI - 1]) >> (ENTRY_SIZE - leftStep));
				curI -= stepI + 1;
				curBit = ENTRY_SIZE - leftStep;
			}
			else{
				assert(curI >= stepI);
				curI -= stepI;
				curBit = 0;
			}
			
			return digit;
		}
		
		virtual bool _hasNext() const{
			return (curBit > 0) || (curI > 0);
		}
	private:
		// coroutine arguments
		const BI &lhs;
		SizeT exp;
		
		// coroutine states
		SizeT curI;
		LogSizeT curBit;
	};// class _LargePower2Radix<Digit, BI, void>
	// signed digit
	template <typename Digit, class BI>
	class _LargePower2Radix<Digit, BI, 
		typename std::enable_if<isSigned<Digit>::value>::type>
		:public DigitProducerCRTP<Digit, BI, _LargePower2Radix<Digit, BI>>{
	private:
		using Unsigned = typename std::make_unsigned<Digit>::type;
		using base = _LargePower2Radix<Unsigned, BI>;
		
		using SizeT = typename BI::SizeT;
	public:
		_LargePower2Radix(const _LargePower2Radix &) = default;
		_LargePower2Radix(_LargePower2Radix &&) = delete;
		
		_LargePower2Radix &operator=(const _LargePower2Radix &) = delete;
		_LargePower2Radix &operator=(_LargePower2Radix &&) = delete;
		
		~_LargePower2Radix() = default;
		
		explicit _LargePower2Radix(const BI &_lhs, SizeT _exp)
			:producer(_lhs, _exp){}
		
		virtual Digit _start(){
			Unsigned tmp = producer._start();
			return tmp;
		}
		
		virtual Digit _next(){
			Unsigned tmp = producer._next();
			return tmp;
		}
		
		virtual bool _hasNext() const{
			return producer._hasNext();
		}
	private:
		base producer;
	};
	
	// TODO: do some meta-programming trick to generate code directly operating on
	// raw pointers
	template <typename Digit, class BI>
	class _ExactDigitExtract:public DigitProducerCRTP<Digit, BI, _ExactDigitExtract<Digit, BI>>{
	private:
		using Ptr = typename BI::Ptr;
		using SizeT = typename BI::SizeT;
	public:
		_ExactDigitExtract(const _ExactDigitExtract &) = default;
		_ExactDigitExtract(_ExactDigitExtract &&) = delete;
		
		_ExactDigitExtract &operator=(const _ExactDigitExtract &) = delete;
		_ExactDigitExtract &operator=(_ExactDigitExtract &&) = delete;
		
		~_ExactDigitExtract() = default;
		
		explicit _ExactDigitExtract(Ptr _begin, SizeT _len)
			:begin(_begin), end(begin + _len), current(end){}
		
		virtual Digit _start(){
			return _next();
		}
		
		virtual Digit _next(){
			assert(_hasNext());
			--current;
			return static_cast<Digit>(*current);
		}
		
		virtual bool _hasNext() const{
			return current != begin;
		}
	private:
		Ptr begin, end, current;
	};// class _ExactDigitExtract
	
	template <typename Digit, class BI>
	class RadixConvertEnumer{
	private:
		using SizeT = typename BI::SizeT;
		
		using Ele = typename BI::Ele;
		constexpr static Ele ENTRY_SIZE = BI::ENTRY_SIZE;
	public:
		using iterator = _type::DigitEnumIterator<Digit, BI>;
		
		RadixConvertEnumer(const RadixConvertEnumer &) = default;
		RadixConvertEnumer(RadixConvertEnumer &&) = default;
		
		RadixConvertEnumer &operator=(const RadixConvertEnumer &) = default;
		RadixConvertEnumer &operator=(RadixConvertEnumer &&) = default;
		
		~RadixConvertEnumer() = default;
		
		// TODO: avoid unnecessary copy/moves
		template <class BIRef, 
			typename std::enable_if<isRLRef<BI, BIRef &&>::value>::type * = nullptr>
		explicit RadixConvertEnumer(BIRef _num, Digit _radix)
			:num(std::forward<BIRef>(_num)), radix(_radix){
			if(!num.positive){
				num.positive = true;
				assert(!num.isZero());
			}
		}
		
		iterator begin() const{
			return beginImpl(std::integral_constant<bool, isSigned<Digit>::value>{});
		}
		
		iterator end() const{
			return iterator(nullptr);
		}
	private:
		// unsigned radix
		iterator beginImpl(std::false_type) const{
			if(radix < 2){
				throw std::domain_error("a radix less than 2 is not accepted.");
				// errno = EDOM;
			}
			
			DigitProducer<Digit, BI> *producer = nullptr;
			if(radix == 10){
				producer = new _DecimalRadix<Digit, BI>(std::move(num));
			}
			if((radix & ((~radix) + 1)) == radix){
				SizeT exp = static_cast<SizeT>(std::ceil(std::log2(radix)));
				if(exp < ENTRY_SIZE){
					producer = new _SmallPower2Radix<Digit, BI>(num, exp);
				}
				if(exp == ENTRY_SIZE){
					producer = new _ExactDigitExtract<Digit, BI>(num.buf.data, num.buf.len);
				}
				if(exp > ENTRY_SIZE){
					producer = new _LargePower2Radix<Digit, BI>(num, exp);
				}
				assert(nullptr != producer);
			}
			if(nullptr == producer){
				producer = new _GenericRadix<Digit, BI>(std::move(num), radix);
			}
			return iterator(producer);
		}
		// signed radix
		iterator beginImpl(std::true_type) const{
			if(radix > 1){
				return beginImpl(std::false_type{});
			}
			
			if(radix < -1){
				// TODO: implement conversion to negative base
				assert(false);
			}
			
			throw std::domain_error("a radix less than 2 and greater than -2 is not accepted.");
			// errno = EDOM;
		}
		
		BI num;
		Digit radix;
	};// class RadixConvertEnumer
	
	template <typename Digit, class BI>
	inline decltype(auto) begin(const RadixConvertEnumer<Digit, BI> &c){
		return c.begin();
	}
	
	template <typename Digit, class BI>
	inline decltype(auto) end(const RadixConvertEnumer<Digit, BI> &c){
		return c.end();
	}
	
};// namespace bignum
#endif // _BIG_INT_OUTPUT_HPP_