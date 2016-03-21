#ifndef _BIG_NUM_HPP_
#error "This header must be included through BigNum.hpp"
#endif // _BIG_NUM_HPP_

#ifndef _BIG_INT_INPUT_HPP_
#define _BIG_INT_INPUT_HPP_

#include <deque>
#include <vector>
#include <type_traits>
#include <utility>
#include <iterator>
#include <cmath>
#include <string>
#include <locale>
#include <istream>

#include "../Libs/BigNumTypeTrait.hpp"
#include "../Libs/BigNumMemory.hpp"
#include "../Libs/BigNumGenerics.hpp"

namespace bignum{
	
	template <typename, class, typename>
	class _GenericAutomatic;
	template <typename, class, typename>
	class _DecimalAutomatic;
	template <typename, class, typename, class>
	class _SmallPower2RadixAutomatic;
	template <typename, class, typename, class>
	class _ExactDigitAutomatic;
	template <typename, class, typename, class>
	class _LargePower2RadixAutomatic;
	
	namespace _type{
		
		template <typename Digit, class BI>
		class DigitReceiver{
		public:
			//virtual DigitReceiver *clone() const = 0;
			virtual void _readDigit(Digit) = 0;
			virtual BI _finish() && = 0;
			virtual BI _finish() const & = 0;
			
			virtual ~DigitReceiver(){}
		};
		
		template <typename Digit, class BI, class Derived>
		class DigitReceiverCRTP:public DigitReceiver<Digit, BI>{
		public:
			/*virtual DigitReceiver *clone() const{
				return new Derived(static_cast<const Derived &>(*this));
			}*/
		};
		
		// uses an ostream-iterator-like prototype
		template <typename Digit, class BI>
		class DigitRecvIterator
			:public std::iterator<
				std::output_iterator_tag,	// iterator_category
				void,						// value_type
				void,						// difference_type
				void,						// pointer
				void>{						// reference
		public:
			using iterator_category = std::output_iterator_tag;
			using value_type = void;
			using difference_type = void;
			using pointer = void;
			using reference = void;
		private:
			using Receiver = DigitReceiver<Digit, BI>;
			
			using SizeT = typename BI::SizeT;
			using LogSizeT = typename BI::LogSizeT;
			
			constexpr static LogSizeT ENTRY_SIZE = BI::ENTRY_SIZE;
		public:
			template <class Derived>
			explicit DigitRecvIterator(Derived *_receiver)
				:receiver(_receiver){}
			
			DigitRecvIterator(const DigitRecvIterator &) = default;
			DigitRecvIterator(DigitRecvIterator &&) = default;
			
			DigitRecvIterator(Digit radix):receiver(nullptr){
				if(radix == 10){
					receiver = std::make_shared<_DecimalAutomatic<Digit, BI, SizeT>>();
					return ;
				}
				if((radix & ((~radix) + 1)) == radix){
					SizeT exp = static_cast<SizeT>(std::ceil(std::log2(radix)));
					if(exp < ENTRY_SIZE){
						receiver = std::make_shared<_SmallPower2RadixAutomatic<Digit, 
							BI, SizeT, std::input_iterator_tag>>(exp);
						return ;
					}
					if(exp == ENTRY_SIZE){
						receiver = std::make_shared<_ExactDigitAutomatic<Digit, 
							BI, SizeT, std::input_iterator_tag>>(exp);
						return ;
					}
					if(exp > ENTRY_SIZE){
						receiver = std::make_shared<_LargePower2RadixAutomatic<Digit, 
							BI, SizeT, std::input_iterator_tag>>(exp);
						return ;
					}
				}
				receiver = std::make_shared<_GenericAutomatic<Digit, BI, SizeT>>(radix);
			}
			
			DigitRecvIterator &operator=(const DigitRecvIterator &) = default;
			DigitRecvIterator &operator=(DigitRecvIterator &&) = default;
			
			~DigitRecvIterator() = default;
			
			DigitRecvIterator &operator=(Digit digit){
				receiver->_readDigit(digit);
				return *this;
			}
			
			DigitRecvIterator &operator*() noexcept{
				return *this;
			}
			
			DigitRecvIterator &operator++() noexcept{
				return *this;
			}
			
			DigitRecvIterator &operator++(int) noexcept{
				return *this;
			}
		private:
			std::shared_ptr<Receiver> receiver;
		};
		
	}; // namespace _type
	
	using _type::isSigned;
	
	template <typename Char, class Trait>
	class _ThousandSepParser{
	private:
		using Stream = std::basic_istream<Char, Trait>;
		using Tr = typename Stream::traits_type;
		using Int = typename Stream::int_type;
		
		using Finder = _utility::CharFindHelper<Char, Int, Tr>;
	public:
		_ThousandSepParser(const _ThousandSepParser &) = delete;
		_ThousandSepParser(_ThousandSepParser &&) = delete;
		
		explicit _ThousandSepParser(Stream &_is, const Char *_digit, const Char *_digitend, 
			const std::ctype<Char> &_cF, const std::numpunct<Char> &_npF)
			:is(_is), ctypeFacet(_cF), npFacet(_npF), digit(_digit), digitend(_digitend), 
			gp(npFacet.grouping()), firstGpLen(0), gpCnt(0), curGpIdx(gp.size()){}
		
		_ThousandSepParser &operator=(const _ThousandSepParser &) = delete;
		_ThousandSepParser &operator=(_ThousandSepParser &&) = delete;
		
		~_ThousandSepParser() = default;
		
		const Char *checkInput(Int _ch){
			if(Tr::eq_int_type(Tr::eof(), _ch)){
				is.setstate(std::ios_base::eofbit);
				return digitend;
			}
			
			const Char *ptr = Finder::find(digit, digitend, _ch);
			if((ptr >= digit) && (ptr < digitend)){
				++gpCnt;
				return ptr;
			}
			if(ctypeFacet.is(std::ctype<Char>::space, Tr::to_char_type(_ch))){
				if(!gp.empty()){
					if(firstGpLen == 0){
						if(firstGpLen > gp[gp.size() - 1]){
							is.setstate(std::ios_base::failbit);
							return nullptr;
						}
						return digitend;
					}
					bool flag = true;
					if(gpCnt != gp[gp.size() - 1]){
						flag = false;
					}
					if(gp.size() > 1){
						if(curGpIdx + 1 != gp.size() - 1){
							flag = false;
						}
					}
					if(!flag){
						is.setstate(std::ios_base::failbit);
						return nullptr;
					}
				}
				return digitend;
			}
			if(Tr::eq_int_type(Tr::to_int_type(npFacet.thousands_sep()), _ch)){
				if(gp.empty()){
					is.setstate(std::ios_base::failbit);
					return nullptr;
				}
				if(gpCnt == 0){
					// ",," or ",......"
					is.setstate(std::ios_base::failbit);
					return nullptr;
				}
				if(firstGpLen != 0){
					if(curGpIdx == gp.size()){
						if(firstGpLen > gpCnt){
							is.setstate(std::ios_base::failbit);
							return nullptr;
						}
						curGpIdx = gp.find(gpCnt);
						if(curGpIdx == std::string::npos){
							is.setstate(std::ios_base::failbit);
							return nullptr;
						}
					}
					else{
						if(curGpIdx == 0){
							bool flag = gp[0] == gpCnt;
							if(gp.size() > 1){
								if(gp[1] == gpCnt){
									flag = true;
									curGpIdx = 1;
								}
							}
							if(!flag){
								is.setstate(std::ios_base::failbit);
								return nullptr;
							}
						}
						else{
							if(curGpIdx + 1 == gp.size()){
								is.setstate(std::ios_base::failbit);
								return nullptr;
							}
							if(gp[curGpIdx + 1] != gpCnt){
								is.setstate(std::ios_base::failbit);
								return nullptr;
							}
							++curGpIdx;
						}
					}
				}
				else{
					firstGpLen = gpCnt;
					curGpIdx = 0;
				}
				return digitend;
			}
			is.setstate(std::ios_base::failbit);
			return nullptr;
		}
		
		const Char *tryNextDigit(){
			Int _ch = is.peek();
			const Char *digitPtr = checkInput(_ch);
			while((digitPtr < digit) || (digitPtr > digitend)){
				if(digitPtr == nullptr){
					return nullptr;
				}
				if(!is){
					return digitPtr;
				}
				is.ignore();
				_ch = is.peek();
				digitPtr = checkInput(_ch);
			}
			is.ignore();
			return digitPtr;
		}
	private:
		Stream &is;
		const std::ctype<Char> &ctypeFacet;
		const std::numpunct<Char> &npFacet;
		const Char *digit, *digitend;
		
		std::string gp;
		char firstGpLen, gpCnt;
		std::string::size_type curGpIdx;
	};
	
	using _type::DigitReceiverCRTP;
	
	template <typename Digit, class BI, typename Diff>
	class _GenericAutomatic
		:public DigitReceiverCRTP<Digit, BI, _GenericAutomatic<Digit, BI, Diff>>{
	private:
		using SizeT = Diff;
		using VSize = typename std::vector<BI>::size_type;
	public:
		explicit _GenericAutomatic(Digit _radix)
			:count(0), baseCap(0), radix(_radix){}
		
		_GenericAutomatic(const _GenericAutomatic &) = delete;
		_GenericAutomatic(_GenericAutomatic &&) = delete;
		
		_GenericAutomatic &operator=(const _GenericAutomatic &) = delete;
		_GenericAutomatic &operator=(_GenericAutomatic &&) = delete;
		
		~_GenericAutomatic() = default;
		
		void _readDigit(Digit digit){
			asset(digit < radix);
			
			numStack.emplace_back(digit);
			++count;
			SizeT _lowbit = count & ((~count) + 1);
			VSize i = 0;
			for(;(1 << i) < _lowbit;++i){
				assert(numStack.size() >= 2);
				BI top1 = std::move(numStack.back());
				numStack.pop_back();
				
				assert(i <= base.size());
				if(i == base.size()){
					if(i == 0){
						base.emplace_back(radix);
					}
					else{
						BI next = base.back();
						next.selfMultiply();
						base.push_back(std::move(next));
					}
				}
				numStack.back().multiplyMedium(base[i]);
				numStack.back().add(std::move(top1));
			}
			
			if(baseCap + 1 < i){
				baseCap = i - 1;
			}
		}
		
		BI _finish() &&{
			if(numStack.empty()){
				return BI(static_cast<Digit>(0));
			}
			
			BI res = std::move(numStack.front());
			numStack.pop_front();
			if(numStack.empty()){
				return res;
			}
			
			VSize i = baseCap;
			//assert((count & (1 << i)) != 0);
			SizeT _lowbit = count & ((~count) + 1);
			assert(_lowbit > 0);
			assert(_lowbit < count);
			assert((1 << i) < count);
			for(SizeT mask = (1 << i);mask >= _lowbit;mask >>= 1, --i){
				assert(!numStack.empty());
				if((mask & count) == 0){
					continue;
				}
				res.multiply(base[i]);
				res.add(std::move(numStack.front()));
				numStack.pop_front();
			}
			assert(numStack.empty());
			return res;
		}
		BI _finish() const &{
			if(numStack.empty()){
				return BI(static_cast<Digit>(0));
			}
			
			using ConstIter = typename std::deque<BI>::const_iterator;
			ConstIter iter = numStack.cbegin();
			BI res = *iter;
			++iter;
			if(iter == numStack.cend()){
				return res;
			}
			
			VSize i = baseCap;
			SizeT _lowbit = count & ((~count) + 1);
			assert(_lowbit > 0);
			assert(_lowbit < count);
			assert((1 << i) < count);
			for(SizeT mask = (1 << i);mask >= _lowbit;mask >>= 1, --i){
				assert(iter != numStack.end());
				if((mask & count) == 0){
					continue;
				}
				res.multiply(base[i]);
				res.add(*iter);
				++iter;
			}
			assert(iter == numStack.end());
			return res;
		}
	private:
		std::vector<BI> base;
		std::deque<BI> numStack;
		
		SizeT count;
		VSize baseCap;
		
		Digit radix;
	};
	
	template <typename Digit, class BI, typename Diff>
	class _DecimalAutomatic
		:public DigitReceiverCRTP<Digit, BI, _DecimalAutomatic<Digit, BI, Diff>>{
	private:
		using SizeT = Diff;
		using VSize = typename std::vector<BI>::size_type;
		
		constexpr static Digit radix = 10;
	public:
		explicit _DecimalAutomatic()
			:decimalBase(BI::getDecimalBase()), count(0), baseCap(0){}
		
		_DecimalAutomatic(const _DecimalAutomatic &) = delete;
		_DecimalAutomatic(_DecimalAutomatic &&) = delete;
		
		_DecimalAutomatic &operator=(const _DecimalAutomatic &) = delete;
		_DecimalAutomatic &operator=(_DecimalAutomatic &&) = delete;
		
		~_DecimalAutomatic() = default;
		
		void _readDigit(Digit digit){
			assert(digit < radix);
			
			numStack.emplace_back(digit);
			++count;
			SizeT _lowbit = count & ((~count) + 1);
			VSize i = 0;
			for(;(1 << i) < _lowbit;++i){
				assert(numStack.size() >= 2);
				BI top1 = std::move(numStack.back());
				numStack.pop_back();
				
				assert(i <= decimalBase.size());
				if(i == decimalBase.size()){
					if(i == 0){
						decimalBase.emplace_back(10);
					}
					else{
						BI next = decimalBase.back();
						next.selfMultiply();
						decimalBase.push_back(std::move(next));
					}
				}
				numStack.back().multiplyMedium(decimalBase[i]);
				//numStack.back().buf.addRaw(std::move(top1.buf));
				numStack.back().add(std::move(top1));
			}
			
			if(baseCap + 1 < i){
				baseCap = i - 1;
			}
		}
		
		BI _finish() &&{
			// elinimate "baseCap==0 && i==0" case
			if(numStack.empty()){
				return BI(static_cast<Digit>(0));
			}
			
			BI res = std::move(numStack.front());
			numStack.pop_front();
			if(numStack.empty()){
				return res;
			}
			
			VSize i = baseCap;
			//assert((count & (1 << i)) != 0);
			SizeT _lowbit = count & ((~count) + 1);
			assert(_lowbit > 0);
			assert(_lowbit < count);
			assert((1 << i) < count);
			for(SizeT mask = (1 << i);mask >= _lowbit;mask >>= 1, --i){
				assert(!numStack.empty());
				if((mask & count) == 0){
					continue;
				}
				res.multiply(decimalBase[i]);
				//res.multiplySmall(decimalBase[i]);
				res.add(std::move(numStack.front()));
				numStack.pop_front();
			}
			assert(numStack.empty());
			return res;
		}
		BI _finish() const &{
			if(numStack.empty()){
				return BI(static_cast<Digit>(0));
			}
			
			using ConstIter = typename std::deque<BI>::const_iterator;
			ConstIter iter = numStack.cbegin();
			BI res = *iter;
			++iter;
			if(iter == numStack.cend()){
				return res;
			}
			
			VSize i = baseCap;
			SizeT _lowbit = count & ((~count) + 1);
			assert(_lowbit > 0);
			assert(_lowbit < count);
			assert((1 << i) < count);
			for(SizeT mask = (1 << i);mask >= _lowbit;mask >>= 1, --i){
				assert(iter != numStack.end());
				if((mask & count) == 0){
					continue;
				}
				res.multiply(decimalBase[i]);
				res.add(*iter);
				++iter;
			}
			assert(iter == numStack.end());
			return res;
		}
	private:
		// not thread-safe
		std::vector<BI> &decimalBase;
		std::deque<BI> numStack;
		
		SizeT count;
		VSize baseCap;
	};
	
	template <typename, class, typename, class>
	class _SmallPower2RadixAutomatic;
	
	template <typename Digit, class BI, typename Diff>
	class _SmallPower2RadixAutomatic<Digit, BI, Diff, std::random_access_iterator_tag>
		:public DigitReceiverCRTP<Digit, BI, 
			_SmallPower2RadixAutomatic<Digit, BI, Diff, std::random_access_iterator_tag>>{
	private:
		using SizeT = typename BI::SizeT;
		using LogSizeT = typename BI::LogSizeT;
		
		static constexpr LogSizeT ENTRY_SIZE = BI::ENTRY_SIZE;
	public:
		explicit _SmallPower2RadixAutomatic(Diff _dist, SizeT _exp)
			:num(typename BI::NullTag{}), dist(_dist), exp(_exp){
			assert(dist > 0);
			assert(exp < ENTRY_SIZE);
			
			SizeT _len = (exp * _dist + ENTRY_SIZE - 1) / ENTRY_SIZE;
			num.buf.setLen(_len);
			num.allocator.allocate(num.buf.data, static_cast<std::size_t>(num.buf.cap));
			SizeT i = 0;
			try{
				for(;i < num.buf.len;++i){
					num.allocator.construct(num.buf.data + i, 0);
				}
			}
			catch(...){
				_utility::destroyAll(num.buf.data, num.buf.data + i, num.allocator);
				num.allocator.deallocate(num.buf.data, static_cast<std::size_t>(num.buf.cap));
				num.buf.data = nullptr;
				num.buf.zeroLen();
				throw;
			}
		}
		
		_SmallPower2RadixAutomatic(const _SmallPower2RadixAutomatic &) = delete;
		_SmallPower2RadixAutomatic(_SmallPower2RadixAutomatic &&) = delete;
		
		_SmallPower2RadixAutomatic &operator=(const _SmallPower2RadixAutomatic &) = delete;
		_SmallPower2RadixAutomatic &operator=(_SmallPower2RadixAutomatic &&) = delete;
		
		~_SmallPower2RadixAutomatic() = default;
		
		void _setDigit(Digit digit, Diff count){
			assert(digit < (1 << exp));
			
			assert(count < dist);
			SizeT curI = ((dist - 1 - count) * exp) / ENTRY_SIZE;
			LogSizeT curBit = ((dist - 1 - count) * exp) % ENTRY_SIZE;
			
			if(curBit + exp <= ENTRY_SIZE){
				num.buf.data[curI] |= digit << curBit;
			}
			else{
				assert(curI + 1 < num.buf.len);
				num.buf.data[curI] |= (digit & ((1 << (ENTRY_SIZE - curBit)) - 1)) << curBit;
				num.buf.data[curI + 1] |= digit >> (ENTRY_SIZE - curBit);
			}
		}
		
		void _readDigit(Digit digit){
			assert(digit < (1 << exp));
			_setDigit(digit, count);
			++count;
		}
		
		BI _finish() &&{
			return std::move(num);
		}
		BI _finish() &{
			return num;
		}
	private:
		Diff dist;
		SizeT exp;
				
		Diff count;
		BI num;
	};
	
	template <typename Digit, class BI, typename Diff>
	class _SmallPower2RadixAutomatic<Digit, BI, Diff, std::input_iterator_tag>
		:public DigitReceiverCRTP<Digit, BI, 
			_SmallPower2RadixAutomatic<Digit, BI, Diff, std::input_iterator_tag>>{
	private:
		using SizeT = typename BI::SizeT;
		using LogSizeT = typename BI::LogSizeT;
		using DB = typename BI::DigitBuffer;
		using Ele = typename BI::Ele;
		
		static constexpr LogSizeT ENTRY_SIZE = BI::ENTRY_SIZE;
	public:
		explicit _SmallPower2RadixAutomatic(SizeT _exp)
			:exp(_exp), num(), curI(0), curBit(ENTRY_SIZE){
			assert(exp < ENTRY_SIZE);
		}
		
		_SmallPower2RadixAutomatic(const _SmallPower2RadixAutomatic &) = delete;
		_SmallPower2RadixAutomatic(_SmallPower2RadixAutomatic &&) = delete;
		
		_SmallPower2RadixAutomatic &operator=(const _SmallPower2RadixAutomatic &) = delete;
		_SmallPower2RadixAutomatic &operator=(_SmallPower2RadixAutomatic &&) = delete;
		
		~_SmallPower2RadixAutomatic() = default;
		
		void _readDigit(Digit digit){
			assert(digit < (1 << exp));
			
			if(curBit >= exp){
				num.buf.data[curI] |= digit << (curBit - exp);
				curBit -= exp;
				return ;
			}
			
			if((curI == 0) && (curBit < exp)){
				DB tmp(&num.allocator, nullptr);
				tmp.setLen(num.buf.len * 2);
				tmp.data = num.allocator.allocate(static_cast<std::size_t>(tmp.cap));
				SizeT i = tmp.len;
				try{
					for(;i > num.buf.len;--i){
						num.allocator.construct(tmp.data + i - 1, num.buf.data[i - num.buf.len - 1]);
					}
					for(i = num.buf.len;i > 0;--i){
						num.allocator.construct(tmp.data + i - 1, 0);
					}
				}
				catch(...){
					_utility::destroyAll(tmp.data + i, tmp.data + tmp.len, num.allocator);
					num.allocator.deallocate(tmp.data, static_cast<std::size_t>(tmp.cap));
					tmp.data = nullptr;
					tmp.zeroLen();
					throw;
				}
				
				curI = num.buf.len;
				
				_utility::destroyAll(num.buf.data, num.buf.data + num.buf.len, num.allocator);
				num.allocator.deallocate(num.buf.data, static_cast<std::size_t>(num.buf.cap));
				num.buf.data = tmp.data;
				num.buf.len = tmp.len;
				num.buf.cap = tmp.cap;
				tmp.data = nullptr;
				tmp.zeroLen();
			}
			
			assert(curI > 0);
			num.buf.data[curI] |= digit >> (exp - curBit);
			num.buf.data[curI - 1] |= (digit & ((1 << (exp - curBit)) - 1)) << (ENTRY_SIZE + curBit - exp);
			curBit = ENTRY_SIZE + curBit - exp;
			--curI;
			
			return ;
		}
		
		BI _finish() &&{
			assert(num.buf.len == num.buf.cap);
			//assert((curI + curBit / ENTRY_SIZE) * 2 < num.buf.len);
			if((curI == 0) && (curBit == 0)){
				return std::move(num);
			}
			for(SizeT i = 0;i + curI < num.buf.len - 1;++i){
				num.buf.data[i] = (num.buf.data[i + curI] >> curBit) | ((num.buf.data[i + curI + 1] & ((1 << curBit) - 1)) << (ENTRY_SIZE - curBit));
			}
			num.buf.data[num.buf.len - curI - 1] = num.buf.data[num.buf.len - 1] >> curBit;
			
			_utility::destroyAll(num.buf.data + num.buf.len - curI, num.buf.data + num.buf.len, num.allocator);
			num.buf.len -= curI;
			
			return std::move(num);
		}
		BI _finish() const &{
			assert(num.buf.len == num.buf.cap);
			if((curI == 0) && (curBit == 0)){
				return num;
			}
			
			Ele leftover = num.buf.data[num.buf.len - 1] >> curBit;
			assert((leftover == 0) || (num.buf.len > curI));
			assert((leftover > 0) || (num.buf.len > curI + 1));
			SizeT _len = (leftover > 0)? (num.buf.len - curI): (num.buf.len - curI - 1);
			
			BI res(typename BI::NullTag{});
			res.buf.setLen(_len);
			res.buf.data = res.allocator.allocate(static_cast<std::size_t>(res.buf.cap));
			SizeT i = 0;
			try{
				for(;i + curI < num.buf.len - 1;++i){
					res.allocator.construct(res.buf.data + i, (num.buf.data[i + curI] >> curBit) | ((num.buf.data[i + curI + 1] & ((1 << curBit) - 1)) << (ENTRY_SIZE - curBit)));
				}
				if(leftover > 0){
					i = res.buf.len - 1;
					res.allocator.construct(res.buf.data + i, leftover);
				}
			}
			catch(...){
				_utility::destroyAll(res.buf.data, res.buf.data + i, res.allocator);
				res.allocator.deallocate(res.buf.data, static_cast<std::size_t>(res.buf.cap));
				res.buf.data = nullptr;
				res.buf.zeroLen();
				throw ;
			}
			
			return res;
		}
	private:
		SizeT exp;
		BI num;
		
		SizeT curI;
		LogSizeT curBit;
	};
	
	template <typename, class, typename, class>
	class _LargePower2RadixAutomatic;
	
	template <typename Digit, class BI, typename Diff>
	class _LargePower2RadixAutomatic<Digit, BI, Diff, std::random_access_iterator_tag>
		:public DigitReceiverCRTP<Digit, BI, 
			_LargePower2RadixAutomatic<Digit, BI, Diff, std::random_access_iterator_tag>>{
	private:
		using SizeT = typename BI::SizeT;
		using LogSizeT = typename BI::LogSizeT;
		using Ele = typename BI::Ele;
		
		static constexpr LogSizeT ENTRY_SIZE = BI::ENTRY_SIZE;
	public:
		explicit _LargePower2RadixAutomatic(Diff _dist, SizeT _exp)
			:dist(_dist), exp(_exp), count(0), num(typename BI::NullTag{}){
			assert(dist > 0);
			assert(exp > ENTRY_SIZE);
			
			SizeT _len = (exp * _dist + ENTRY_SIZE - 1) / ENTRY_SIZE;
			num.buf.setLen(_len);
			num.buf.data = num.allocator.allocate(static_cast<std::size_t>(num.buf.cap));
			SizeT i = 0;
			try{
				for(;i < num.buf.len;++i){
					num.allocator.construct(num.buf.data + i, 0);
				}
			}
			catch(...){
				_utility::destroyAll(num.buf.data, num.buf.data + i, num.allocator);
				num.allocator.deallocate(num.buf.data, static_cast<std::size_t>(num.buf.cap));
				num.buf.data = nullptr;
				num.buf.zeroLen();
				throw;
			}
		}
		
		_LargePower2RadixAutomatic(const _LargePower2RadixAutomatic &) = delete;
		_LargePower2RadixAutomatic(_LargePower2RadixAutomatic &&) = delete;
		
		_LargePower2RadixAutomatic &operator=(const _LargePower2RadixAutomatic &) = delete;
		_LargePower2RadixAutomatic &operator=(_LargePower2RadixAutomatic &&) = delete;
		
		~_LargePower2RadixAutomatic() = default;
		
		void _setDigit(Digit digit, Diff count){
			assert(digit < (1 << exp));
			
			assert(count < dist);
			SizeT curI = ((dist - 1 - count) * exp) / ENTRY_SIZE;
			LogSizeT curBit = ((dist - 1 - count) * exp) % ENTRY_SIZE;
			
			assert(curI + 1 < num.buf.len);
			SizeT stepI = (exp + curBit - ENTRY_SIZE) / ENTRY_SIZE;
			LogSizeT leftStep = (exp + curBit - ENTRY_SIZE) % ENTRY_SIZE;
			num.buf.data[curI] |= static_cast<Ele>((digit & ((1 << (ENTRY_SIZE - curBit)) - 1)) << curBit);
			digit >>= ENTRY_SIZE - curBit;
			for(SizeT i = 1;i <= stepI;++i){
				num.buf.data[curI + i] |= static_cast<Ele>(digit & ((1 << ENTRY_SIZE) - 1));
				digit >>= ENTRY_SIZE;
			}
			assert(digit < (1 << leftStep));
			if(leftStep > 0){
				num.buf.data[curI + stepI + 1] |= static_cast<Ele>(digit);
			}
		}
		
		void _readDigit(Digit digit){
			assert(digit < (1 << exp));
			_setDigit(digit, count);
			++count;
		}
		
		BI _finish() &&{
			return std::move(num);
		}
		BI _finish() &{
			return num;
		}
	private:
		Diff dist;
		SizeT exp;
				
		Diff count;
		BI num;
	};
	
	template <typename Digit, class BI, typename Diff>
	class _LargePower2RadixAutomatic<Digit, BI, Diff, std::input_iterator_tag>
		:public DigitReceiverCRTP<Digit, BI, 
			_LargePower2RadixAutomatic<Digit, BI, Diff, std::input_iterator_tag>>{
	private:
		using SizeT = typename BI::SizeT;
		using LogSizeT = typename BI::LogSizeT;
		using DB = typename BI::DigitBuffer;
		using Ele = typename BI::Ele;
		
		static constexpr LogSizeT ENTRY_SIZE = BI::ENTRY_SIZE;
	public:
		explicit _LargePower2RadixAutomatic(SizeT _exp)
			:exp(_exp), num(), curI(0), curBit(ENTRY_SIZE){
			assert(exp > ENTRY_SIZE);
		}
		
		_LargePower2RadixAutomatic(const _LargePower2RadixAutomatic &) = delete;
		_LargePower2RadixAutomatic(_LargePower2RadixAutomatic &&) = delete;
		
		_LargePower2RadixAutomatic &operator=(const _LargePower2RadixAutomatic &) = delete;
		_LargePower2RadixAutomatic &operator=(_LargePower2RadixAutomatic &&) = delete;
		
		~_LargePower2RadixAutomatic() = default;
		
		void _readDigit(Digit digit){
			assert(digit < (1 << exp));
			
			SizeT resLen = num.buf.len, space = 0;
			while((curI + space) * ENTRY_SIZE + curBit < exp){
				space += resLen;
				resLen <<= 1;
			}
			
			if(space != 0){
				DB tmp(&num.allocator, nullptr);
				tmp.setLen(resLen);
				tmp.data = num.allocator.allocate(static_cast<std::size_t>(tmp.cap));
				SizeT i = tmp.len;
				try{
					for(;i > space;--i){
						num.allocator.allocate(tmp.data + i - 1, num.buf.data[i - num.buf.len - 1]);
					}
					for(i = space;i > 0;--i){
						num.allocator.allocate(tmp.data + i - 1, 0);
					}
				}
				catch(...){
					_utility::destroyAll(tmp.data + i, tmp.data + tmp.len, num.allocator);
					num.allocator.deallocate(tmp.data, static_cast<std::size_t>(tmp.cap));
					tmp.data = nullptr;
					tmp.zeroLen();
					throw;
				}
				
				curI += space;
				
				_utility::destroyAll(num.buf.data, num.buf.data + num.buf.len, num.allocator);
				num.allocator.deallocate(num.buf.data, static_cast<std::size_t>(num.buf.cap));
				num.buf.data = tmp.data;
				num.buf.len = tmp.len;
				num.buf.cap = tmp.cap;
				tmp.data = nullptr;
				tmp.zeroLen();
			}
			
			SizeT stepI = (exp - curBit) / ENTRY_SIZE;
			LogSizeT leftStep = (exp - curBit) % ENTRY_SIZE;
			if(leftStep > 0){
				assert(curI > stepI);
				num.buf.data[curI - stepI - 1] |= static_cast<Ele>(digit & ((1 << leftStep) - 1)) << (ENTRY_SIZE - leftStep);
				digit >>= leftStep;
			}
			else{
				assert(curI >= stepI);
			}
			for(SizeT i = stepI;i > 0;--i){
				num.buf.data[curI - stepI] |= static_cast<Ele>(digit & ((1 << ENTRY_SIZE) - 1));
				digit >>= ENTRY_SIZE;
			}
			assert(digit < (1 << curBit));
			num.buf.data[curI] |= static_cast<Ele>(digit);
			
			if(leftStep > 0){
				curI -= stepI + 1;
				curBit = ENTRY_SIZE - leftStep;
			}
			else{
				curI -=stepI;
				curBit = 0;
			}
		}
		
		BI _finish() &&{
			assert(num.buf.len == num.buf.cap);
			assert((curI + curBit / ENTRY_SIZE) * 2 < num.buf.len);
			if((curI == 0) && (curBit == 0)){
				return std::move(num);
			}
			for(SizeT i = 0;i + curI < num.buf.len - 1;++i){
				num.buf.data[i] = (num.buf.data[i + curI] >> curBit) | ((num.buf.data[i + curI + 1] & ((1 << curBit) - 1)) << (ENTRY_SIZE - curBit));
			}
			num.buf.data[num.buf.len - curI - 1] = num.buf.data[num.buf.len - 1] >> curBit;
			
			_utility::destroyAll(num.buf.data + num.buf.len - curI, num.buf.data + num.buf.len, num.allocator);
			num.buf.len -= curI;
			
			return std::move(num);
		}
		BI _finish() const &{
			assert(num.buf.len == num.buf.cap);
			assert((curI + curBit / ENTRY_SIZE) * 2 < num.buf.len);
			if((curI == 0) && (curBit == 0)){
				return num;
			}
			
			Ele leftover = num.buf.data[num.buf.len - 1] >> curBit;
			assert((leftover == 0) || (num.buf.len > curI));
			assert((leftover > 0) || (num.buf.len > curI + 1));
			SizeT _len = (leftover > 0)? (num.buf.len - curI): (num.buf.len - curI - 1);
			
			BI res(typename BI::NullTag{});
			res.buf.setLen(_len);
			res.buf.data = res.allocator.allocate(static_cast<std::size_t>(res.buf.cap));
			SizeT i = 0;
			try{
				for(;i + curI < num.buf.len - 1;++i){
					res.allocator.construct(res.buf.data + i, (num.buf.data[i + curI] >> curBit) | ((num.buf.data[i + curI + 1] & ((1 << curBit) - 1)) << (ENTRY_SIZE - curBit)));
				}
				if(leftover > 0){
					i = res.buf.len - 1;
					res.allocator.construct(res.buf.data + i, leftover);
				}
			}
			catch(...){
				_utility::destroyAll(res.buf.data, res.buf.data + i, res.allocator);
				res.allocator.deallocate(res.buf.data, static_cast<std::size_t>(res.buf.cap));
				res.buf.data = nullptr;
				res.buf.zeroLen();
				throw ;
			}
			return res;
		}
	private:
		SizeT exp;
		BI num;
		
		SizeT curI;
		LogSizeT curBit;
	};
	
	template <typename, class, typename, class>
	class _ExactDigitAutomatic;
	
	template <typename Digit, class BI, typename Diff>
	class _ExactDigitAutomatic<Digit, BI, Diff, std::random_access_iterator_tag>
		:public DigitReceiverCRTP<Digit, BI, 
			_ExactDigitAutomatic<Digit, BI, Diff, std::random_access_iterator_tag>>{
	private:
		using SizeT = typename BI::SizeT;
		using Ele = typename BI::Ele;
		using LogSizeT = typename BI::LogSizeT;
		
		static constexpr LogSizeT ENTRY_SIZE = BI::ENTRY_SIZE;
	public:
		explicit _ExactDigitAutomatic(Diff _dist)
			:dist(_dist), count(0), num(typename BI::NullTag{}){
			assert(_dist > 0);
			
			num.buf.setLen(dist);
			num.buf.data = num.allocator.allocate(static_cast<std::size_t>(num.buf.cap));
			SizeT i = 0;
			try{
				for(;i < num.buf.len;++i){
					num.allocator.construct(num.buf.data + i, 0);
				}
			}
			catch(...){
				_utility::destroyAll(num.buf.data, num.buf.data + i, num.allocator);
				num.allocator.deallocate(num.buf.data, static_cast<std::size_t>(num.buf.cap));
				num.buf.data = nullptr;
				num.buf.zeroLen();
				throw;
			}
		}
		
		_ExactDigitAutomatic(const _ExactDigitAutomatic &) = delete;
		_ExactDigitAutomatic(_ExactDigitAutomatic &&) = delete;
		
		_ExactDigitAutomatic &operator=(const _ExactDigitAutomatic &) = delete;
		_ExactDigitAutomatic &operator=(_ExactDigitAutomatic &&) = delete;
		
		~_ExactDigitAutomatic() = default;
		
		void _setDigit(Digit digit, Diff count){
			assert(digit < (1 << ENTRY_SIZE));
			num.buf.data[count] = digit;
		}
		
		void _readDigit(Digit digit){
			assert(digit < (1 << ENTRY_SIZE));
			_setDigit(digit, count);
			++count;
		}
		
		BI _finish() &&{
			return std::move(num);
		}
		BI _finish() &{
			return num;
		}
	private:
		Diff dist, count;
		BI num;
	};
	
	template <typename Digit, class BI, typename Diff>
	class _ExactDigitAutomatic<Digit, BI, Diff, std::input_iterator_tag>
		:public DigitReceiverCRTP<Digit, BI, 
			_ExactDigitAutomatic<Digit, BI, Diff, std::input_iterator_tag>>{
	private:
		using SizeT = typename BI::SizeT;
		using Ele = typename BI::Ele;
		using DB = typename BI::DigitBuffer;
		using LogSizeT = typename BI::LogSizeT;
		
		static constexpr LogSizeT ENTRY_SIZE = BI::ENTRY_SIZE;
	public:
		explicit _ExactDigitAutomatic()
			:num(0), curI(1){}
		
		_ExactDigitAutomatic(const _ExactDigitAutomatic &) = delete;
		_ExactDigitAutomatic(_ExactDigitAutomatic &&) = delete;
		
		_ExactDigitAutomatic &operator=(const _ExactDigitAutomatic &) = delete;
		_ExactDigitAutomatic &operator=(_ExactDigitAutomatic &&) = delete;
		
		~_ExactDigitAutomatic() = default;
		
		void _readDigit(Digit digit){
			assert(digit < (1 << ENTRY_SIZE));
			
			if(curI == 0){
				DB tmp(&num.allocator, nullptr);
				tmp.setLen(num.buf.len * 2);
				tmp.data = num.allocator.allocate(static_cast<std::size_t>(tmp.cap));
				SizeT i = tmp.len;
				try{
					for(;i > num.buf.len;--i){
						num.allocator.allocate(tmp.data + i - 1, num.buf.data[i - num.buf.len - 1]);
					}
					for(i = num.buf.len;i > 0;--i){
						num.allocator.allocate(tmp.data + i - 1, 0);
					}
				}
				catch(...){
					_utility::destroyAll(tmp.data + i, tmp.data + tmp.len, num.allocator);
					num.allocator.deallocate(tmp.data, static_cast<std::size_t>(tmp.cap));
					tmp.data = nullptr;
					tmp.zeroLen();
					throw;
				}
				
				curI = num.buf.len;
				
				_utility::destroyAll(num.buf.data, num.buf.data + num.buf.len, num.allocator);
				num.allocator.deallocate(num.buf.data, static_cast<std::size_t>(num.buf.cap));
				num.buf.data = tmp.data;
				num.buf.len = tmp.len;
				num.buf.cap = tmp.cap;
				tmp.data = nullptr;
				tmp.zeroLen();
			}
			
			assert(curI > 0);
			num.buf.data[curI - 1] = digit;
			--curI;
		}
		
		BI _finish() &&{
			assert(num.buf.len == num.buf.cap);
			if(curI == 0){
				return std::move(num);
			}
			for(SizeT i = 0;i + curI < num.buf.len;++i){
				num.buf.data[i] = num.buf.data[i + curI];
			}
			
			_utility::destroyAll(num.buf.data + num.buf.len - curI, num.buf.data + num.buf.len, num.allocator);
			num.buf.len -= curI;
			
			return std::move(num);
		}
		BI _finish() const &{
			assert(num.buf.len == num.buf.cap);
			if(curI == 0){
				return num;
			}
			
			BI res(typename BI::NullTag{});
			res.buf.setLen(num.buf.len - curI);
			res.buf.data = res.allocator.allocate(static_cast<std::size_t>(res.buf.cap));
			SizeT i = 0;
			try{
				for(;i + curI < num.buf.len;++i){
					res.buf.data[i] = num.buf.data[i + curI];
				}
			}
			catch(...){
				_utility::destroyAll(res.buf.data, res.buf.data + i, res.allocator);
				res.allocator.deallocate(res.buf.data, static_cast<std::size_t>(res.buf.cap));
				res.buf.data = nullptr;
				res.buf.zeroLen();
				throw ;
			}
			
			return res;
		}
	private:
		BI num;
		SizeT curI;
	};
	
	template <typename Digit, class BI>
	class RadixConvertRecver{
	public:
		//using iterator = _type::DigitRecvIterator<Digit, BI>;
		using Unsigned = typename std::make_unsigned<Digit>::type;
		using SizeT = typename BI::SizeT;
		using LogSizeT = typename BI::LogSizeT;
		
		static constexpr LogSizeT ENTRY_SIZE = BI::ENTRY_SIZE;
	public:
		explicit RadixConvertRecver(Digit _radix)
			:radix(_radix){}
		
		RadixConvertRecver(const RadixConvertRecver &) = default;
		RadixConvertRecver(RadixConvertRecver &&) = default;
		
		RadixConvertRecver &operator=(const RadixConvertRecver &) = default;
		RadixConvertRecver &operator=(RadixConvertRecver &&) = default;
		
		~RadixConvertRecver() = default;
		
		template <typename Iter>
		BI readDigits(Iter _begin, Iter _end){
			if(radix == 10){
				return readDigitsDecimal(_begin, _end);
			}
			if((radix & ((~radix) + 1)) == radix){
				SizeT exp = static_cast<SizeT>(std::ceil(std::log2(radix)));
				using iterTag = typename std::iterator_traits<Iter>::iterator_category;
				if(exp < ENTRY_SIZE){
					return readDigitsSmallPowerOf2(exp, _begin, _end, iterTag{});
				}
				if(exp == ENTRY_SIZE){
					return readDigitsExact(exp, _begin, _end, iterTag{});
				}
				if(exp > ENTRY_SIZE){
					return readDigitsLargePowerOf2(exp, _begin, _end, iterTag{});
				}
			}
			
			return readDigitsGeneric(radix, _begin, _end);
		}
		
		/*iterator getIterator() const{
			
		}*/
	private:
		template <typename Iter>
		BI readDigitsDecimal(Iter _begin, Iter _end){
			using Diff = typename std::iterator_traits<Iter>::difference_type;
			_DecimalAutomatic<Unsigned, BI, Diff> automatic;
			for(Iter iter = _begin;iter != _end;++iter){
				automatic._readDigit(static_cast<Unsigned>(*iter));
			}
			return std::move(automatic)._finish();
		}
		
		template <typename Iter>
		BI readDigitsSmallPowerOf2(SizeT exp, Iter _begin, Iter _end, std::random_access_iterator_tag){
			using Diff = typename std::iterator_traits<Iter>::difference_type;
			Diff dist = std::distance(_begin, _end);
			_SmallPower2RadixAutomatic<Unsigned, BI, Diff, std::random_access_iterator_tag> automatic(dist, exp);
			for(Iter iter = _begin;iter != _end;++iter){
				automatic._readDigit(static_cast<Unsigned>(*iter));
			}
			return std::move(automatic)._finish();
		}
		template <typename Iter>
		BI readDigitsSmallPowerOf2(SizeT exp, Iter _begin, Iter _end, std::input_iterator_tag){
			using Diff = typename std::iterator_traits<Iter>::difference_type;
			_SmallPower2RadixAutomatic<Unsigned, BI, Diff, std::input_iterator_tag> automatic(exp);
			for(Iter iter = _begin;iter != end;++iter){
				automatic._readDigit(static_cast<Unsigned>(*iter));
			}
			return std::move(automatic)._finish();
		}
		
		template <typename Iter>
		BI readDigitsExact(Iter _begin, Iter _end, std::random_access_iterator_tag){
			using Diff = typename std::iterator_traits<Iter>::difference_type;
			Diff dist = std::distance(_begin, _end);
			_ExactDigitAutomatic<Unsigned, BI, Diff, std::random_access_iterator_tag> automatic(dist);
			for(Iter iter = _begin;iter != _end;++iter){
				automatic._readDigit(static_cast<Unsigned>(*iter));
			}
			return std::move(automatic)._finish();
		}
		template <typename Iter>
		BI readDigitsExact(Iter _begin, Iter _end, std::input_iterator_tag){
			using Diff = typename std::iterator_traits<Iter>::difference_type;
			_ExactDigitAutomatic<Unsigned, BI, Diff, std::input_iterator_tag> automatic;
			for(Iter iter = _begin;iter != _end;++iter){
				automatic._readDigit(static_cast<Unsigned>(*iter));
			}
			return std::move(automatic)._finish();
		}
		
		template <typename Iter>
		BI readDigitsLargePowerOf2(SizeT exp, Iter _begin, Iter _end, std::random_access_iterator_tag){
			using Diff = typename std::iterator_traits<Iter>::difference_type;
			Diff dist = std::distance(_begin, _end);
			_LargePower2RadixAutomatic<Unsigned, BI, Diff,std::random_access_iterator_tag> automatic(dist, exp);
			for(Iter iter = _begin;iter != _end;++iter){
				automatic._readDigit(static_cast<Unsigned>(*iter));
			}
			return std::move(automatic)._finish();
		}
		template <typename Iter>
		BI readDigitsLargePowerOf2(SizeT exp, Iter _begin, Iter _end, std::input_iterator_tag){
			using Diff = typename std::iterator_traits<Iter>::difference_type;
			_LargePower2RadixAutomatic<Unsigned, BI, Diff,std::random_access_iterator_tag> automatic(exp);
			for(Iter iter = _begin;iter != _end;++iter){
				automatic._readDigit(static_cast<Unsigned>(*iter));
			}
			return std::move(automatic)._finish();
		}
		
		template <typename Iter>
		BI readDigitsGeneric(Digit radix, Iter _begin, Iter _end){
			using Diff = typename std::iterator_traits<Iter>::difference_type;
			_GenericAutomatic<Unsigned, BI, Diff> automatic(radix);
			for(Iter iter = _begin;iter != _end;++iter){
				automatic._readDigit(static_cast<Unsigned>(*iter));
			}
			return std::move(automatic)._finish();
		}
		
		Digit radix;
	};
	
	namespace _type{
		
		template <class Automatic, char ...Args>
		inline void inputDigitImpl(Automatic &automatic, StaticList<std::integral_constant<char, Args>...>){
			PackExpandHelper<int> helper({0, (automatic._readDigit(Args), 0)...});
		}
		
		template <typename>
		struct NotSep;
		template <char ch>
		struct NotSep<std::integral_constant<char, ch>>
			:public std::true_type{};
		template <>
		struct NotSep<std::integral_constant<char, '\''>>
			:public std::false_type{};
		
		template <typename>
		struct IsBin;
		template <char Ch>
		struct IsBin<std::integral_constant<char, Ch>>
			:public std::integral_constant<bool, (Ch > '0' && Ch <= '1')>{};
		
		template <typename>
		struct IsOct;
		template <char Ch>
		struct IsOct<std::integral_constant<char, Ch>>
			:public std::integral_constant<bool, (Ch > '0' && Ch <= '7')>{};
		
		template <typename>
		struct IsDec;
		template <char Ch>
		struct IsDec<std::integral_constant<char, Ch>>
			:public std::integral_constant<bool, (Ch > '0' && Ch <= '9')>{};
		
		template <typename>
		struct IsHex;
		template <char Ch>
		struct IsHex<std::integral_constant<char, Ch>>
			:public std::integral_constant<bool, ((Ch > '0' && Ch <= '9') || 
				(Ch > 'a' && Ch <= 'f') || (Ch > 'A' && Ch <= 'F'))>{};
		
		template <typename>
		struct ToBin;
		template <char Ch>
		struct ToBin<std::integral_constant<char, Ch>>{
			using type = std::integral_constant<char, Ch - '0'>;
		};
		
		template <typename>
		struct ToOct;
		template <char Ch>
		struct ToOct<std::integral_constant<char, Ch>>{
			using type = std::integral_constant<char, Ch - '0'>;
		};
		
		template <typename>
		struct ToDec;
		template <char Ch>
		struct ToDec<std::integral_constant<char, Ch>>{
			using type = std::integral_constant<char, Ch - '0'>;
		};
		
		template <typename>
		struct ToHex;
		template <char Ch>
		struct ToHex<std::integral_constant<char, Ch>>{
			using type = std::integral_constant<char, (Ch > '0' && Ch <= '9')? (Ch - '0'): 
				((Ch > 'a' && Ch <= 'f')? (Ch - 'a'): (Ch - 'A'))>;
		};
		
		template <class, class>
		struct LiteralParser{};
		// hex
		template <class BI, char ...Args>
		struct LiteralParser<BI, StaticList<std::integral_constant<char, '0'>, 
			std::integral_constant<char, 'x'>, 
			std::integral_constant<char, Args>...>>{
		private:
			using CharSeq = typename Filter<NotSep, StaticList<std::integral_constant<char, Args>...>>::type;
			using DigitSeq = typename Map<ToHex, typename Check<IsHex, CharSeq>::type>::type;
			
			using Automatic = typename CompareCond<typename BI::LogSizeT, 4, BI::ENTRY_SIZE, 
				_LargePower2RadixAutomatic<char, BI, std::size_t, std::random_access_iterator_tag>, 
				_ExactDigitAutomatic<char, BI, std::size_t, std::random_access_iterator_tag>,
				_SmallPower2RadixAutomatic<char, BI, std::size_t, std::random_access_iterator_tag>>::type;
		public:
			inline static BI getBI(){
				Automatic automatic(sizeof...(Args), 4);
				inputDigitImpl(automatic, DigitSeq{});
				return std::move(automatic)._finish();
			}
		};
		template <class BI, char ...Args>
		struct LiteralParser<BI, StaticList<std::integral_constant<char, '0'>, 
			std::integral_constant<char, 'X'>, 
			std::integral_constant<char, Args>...>>
			:public LiteralParser<BI, StaticList<std::integral_constant<char, '0'>, 
				std::integral_constant<char, 'x'>, 
				std::integral_constant<char, Args>...>>{};
		// oct
		template <class BI, char ...Args>
		struct LiteralParser<BI, StaticList<std::integral_constant<char, '0'>, 
			std::integral_constant<char, Args>...>>{
		private:
			using CharSeq = typename Filter<NotSep, StaticList<std::integral_constant<char, Args>...>>::type;
			using DigitSeq = typename Map<ToOct, typename Check<IsOct, CharSeq>::type>::type;
			
			using Automatic = typename CompareCond<typename BI::LogSizeT, 3, BI::ENTRY_SIZE, 
				_LargePower2RadixAutomatic<char, BI, std::size_t, std::random_access_iterator_tag>, 
				_ExactDigitAutomatic<char, BI, std::size_t, std::random_access_iterator_tag>,
				_SmallPower2RadixAutomatic<char, BI, std::size_t, std::random_access_iterator_tag>>::type;
		public:
			inline static BI getBI(){
				Automatic automatic(sizeof...(Args), 3);
				inputDigitImpl(automatic, DigitSeq{});
				return std::move(automatic)._finish();
			}
		};
		// bin
		template <class BI, char ...Args>
		struct LiteralParser<BI, StaticList<std::integral_constant<char, '0'>, 
			std::integral_constant<char, 'b'>, 
			std::integral_constant<char, Args>...>>{
		private:
			using CharSeq = typename Filter<NotSep, StaticList<std::integral_constant<char, Args>...>>::type;
			using DigitSeq = typename Map<ToBin, typename Check<IsBin, CharSeq>::type>::type;
			
			using Automatic = typename CompareCond<typename BI::LogSizeT, 1, BI::ENTRY_SIZE, 
				_LargePower2RadixAutomatic<char, BI, std::size_t, std::random_access_iterator_tag>, 
				_ExactDigitAutomatic<char, BI, std::size_t, std::random_access_iterator_tag>,
				_SmallPower2RadixAutomatic<char, BI, std::size_t, std::random_access_iterator_tag>>::type;
		public:
			inline static BI getBI(){
				Automatic automatic(sizeof...(Args), 1);
				inputDigitImpl(automatic, DigitSeq{});
				return std::move(automatic)._finish();
			}
		};
		template <class BI, char ...Args>
		struct LiteralParser<BI, StaticList<std::integral_constant<char, '0'>, 
			std::integral_constant<char, 'B'>, 
			std::integral_constant<char, Args>...>>
			:public LiteralParser<BI, StaticList<std::integral_constant<char, '0'>, 
				std::integral_constant<char, 'b'>, 
				std::integral_constant<char, Args>...>>{};
		// dec
		template <class BI, char ...Args>
		struct LiteralParser<BI, StaticList<std::integral_constant<char, Args>...>>{
			using CharSeq = typename Filter<NotSep, StaticList<std::integral_constant<char, Args>...>>::type;
			using DigitSeq = typename Map<ToDec, typename Check<IsDec, CharSeq>::type>::type;
			
			using Automatic = _DecimalAutomatic<char, BI, std::size_t>;
			inline static BI getBI(){
				Automatic automatic;
				inputDigitImpl(automatic, DigitSeq{});
				return std::move(automatic)._finish();
			}
		};
		
	}; // namespace _type
	
}; // namespace bignum
#endif // _BIG_INT_INPUT_HPP_