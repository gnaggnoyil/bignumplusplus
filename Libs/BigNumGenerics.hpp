#ifndef _BIG_NUM_HPP_
#error "This header must be included through BigNum.hpp"
#endif // _BIG_NUM_HPP_

#ifndef _BIG_NUM_GENERICS_HPP_
#define _BIG_NUM_GENERICS_HPP_

#include <type_traits>
#include <utility>
#include <functional>
#include <cassert>
#include <memory>
#include <algorithm>

#include "BigNumTypeTrait.hpp"

namespace bignum{
	
	namespace _utility{
		
		template <typename Char, typename Int, class Trait>
		class CharFindHelper{
		public:
			static const Char *find(const Char *begin, const Char *end, Int ch) noexcept{
				return std::find_if(begin, end, [&ch](const Char &ele) -> bool{
					return Trait::eq_int_type(ch, Trait::to_int_type(ele));
				});
			}
		};
		
	}; // namespace _utility
	
	namespace _type{
		
		template <typename Callable, typename SFINAE = void>
		struct CallableHolder;
		
		template <class Functor>
		struct CallableHolder<Functor, 
			typename std::enable_if<std::is_class<Functor>::value>::type>
			:private Functor{
		public:
			using Functor::operator();
			
			template <class FunctorRef, 
				typename std::enable_if<isRLRef<Functor, FunctorRef &&>::value>::type * = nullptr>
			explicit constexpr CallableHolder(FunctorRef &&_func)
				:Functor(std::forward<FunctorRef>(_func)){}
		};
		
		template <typename Ret, typename... Args>
		struct CallableHolder<Ret (*)(Args...)>{
		private:
			using FuncPtr = Ret (*)(Args...);
		public:
			explicit constexpr CallableHolder(FuncPtr _func)
				:func(_func){}
				
			Ret operator()(Args... args) const{
				return (*func)(std::forward<Args>(args)...);
			}
		private:
			FuncPtr func;
		};
		
		template <typename Ret, class T>
		struct CallableHolder<Ret T::*, 
			typename std::enable_if<std::is_class<T>::value>::type>
			:private CallableHolder<decltype(std::mem_fn(std::declval<Ret T::*>()))>{
		private:
			using MemPtr = Ret T::*;
			using base = CallableHolder<decltype(std::mem_fn(std::declval<MemPtr>()))>;
		public:
			explicit constexpr CallableHolder(MemPtr _func)
				:base(std::mem_fn(_func)){}
				
			using base::operator();
		};
		
		template <typename... Callable>
		struct CallableOverloader;
		
		template <typename Callable>
		struct CallableOverloader<Callable>:private CallableHolder<Callable>{
		private:
			using base = CallableHolder<Callable>;
		public:
			template <typename CallableRef, 
				typename std::enable_if<isRLRef<Callable, CallableRef &&>::value>::type * = nullptr>
			constexpr CallableOverloader(CallableRef &&_func)
				:base(std::forward<CallableRef>(_func)){}
			
			using base::operator();
		};
		
		template <typename Callable1, typename Callable2, typename... Callables>
		struct CallableOverloader<Callable1, Callable2, Callables...>
			:private CallableHolder<Callable1>, private CallableOverloader<Callable2, Callables...>{
		public:
			template <typename CallableRef1, typename CallableRef2, typename... CallableRefs,
				typename std::enable_if<conj<isRLRef<Callable1, CallableRef1 &&>, 
					isRLRef<Callable2, CallableRef2 &&>, 
					isRLRef<Callables, CallableRefs &&>...>::value>::type * = nullptr>
			constexpr CallableOverloader(CallableRef1 &&cr1, CallableRef2 &&cr2, CallableRefs && ...crs)
				:CallableHolder<Callable1>(std::forward<CallableRef1>(cr1)),
					CallableOverloader<Callable2, Callables...>(std::forward<CallableRef2>(cr2), 
						std::forward<CallableRefs>(crs)...){}
			
			using CallableHolder<Callable1>::operator();
			using CallableOverloader<Callable2, Callables...>::operator();
		};
		
		template <typename... Callables>
		inline CallableOverloader<Callables...> makeOverload(Callables... callables){
			return CallableOverloader<Callables...>(std::forward<Callables>(callables)...);
		}
		
#define _BIG_NUM_GENERIC_LITERAL_(type, literal) \
	::bignum::_type::makeOverload([](char) -> decltype(auto){return literal;}, \
		[](wchar_t) -> decltype(auto){return L ## literal;}, \
		[](char16_t) -> decltype(auto){return u ## literal;}, \
		[](char32_t) -> decltype(auto){return U ## literal;})(type())
		
	};  // namespace _type
	
}; // namespace bignum

#endif // _BIG_NUM_GENERICS_HPP_