#ifndef _BIG_NUM_HPP_
#error "This header must be included through BigNum.hpp"
#endif // _BIG_NUM_HPP_

#ifndef _BIG_NUM_MEMORY_HPP_
#define _BIG_NUM_MEMORY_HPP_

#include <type_traits>
#include <utility>

#include "BigNumTypeTrait.hpp"

namespace bignum{
	
	namespace _utility{
		
		using _type::isRLRef;
		
		template <class>
		struct destroyHelper{};
		// case for trivially destructible types
		template <>
		struct destroyHelper<std::true_type>{
			template <typename ForIter, class Alloc>
			static void destroy(ForIter, ForIter, Alloc &){
				// do nothing
			}
		};
		// case for non trivially destructible types
		template <>
		struct destroyHelper<std::false_type>{
			template <typename ForIter, class Alloc>
			static void destroy(ForIter _start, ForIter _last, Alloc &_alloc){
				using valueT = typename std::iterator_traits<ForIter>::value_type;
				for(;_start != _last;++_start){
					_alloc.template destroy<valueT>(_alloc.address(*_start));
				}
			}
		};
		
		template <typename ForIter, class Alloc>
		void destroyAll(ForIter _start, ForIter _last, Alloc &_alloc){
			using valueT = typename std::iterator_traits<ForIter>::value_type;
			destroyHelper<typename std::is_trivially_destructible<valueT>::type>::destroy(_start, _last, _alloc);
		}
		
		template <typename Callable, typename SFINAE1 = void, typename SFINAE2 = void>
		class FinalizerImpl;
		
		template <class Functor>
		class FinalizerImpl<Functor,  
			typename std::enable_if<std::is_class<Functor>::value>::type, 
			typename std::result_of<Functor()>::type>
			:private Functor{
		private:
			using base = Functor;
		public:
			template <typename FunctorRef, 
				typename std::enable_if<isRLRef<Functor, FunctorRef>::value>::type * = nullptr>
			explicit FinalizerImpl(FunctorRef &&_func)
				:Functor(std::forward<FunctorRef>(_func)){}
			
			~FinalizerImpl(){
				base::operator()();
			}
		};
		
		template <typename Ret>
		class FinalizerImpl<Ret (*)()>{
		private:
			using FuncPtr = Ret (*)();
		public:
			explicit FinalizerImpl(FuncPtr _func)
				:func(_func){}
			
			~FinalizerImpl(){
				(*func)();
			}
		private:
			FuncPtr func;
		};
		
		template <typename Callable>
		class Finalizer:public FinalizerImpl<Callable>{
		private:
			using base = FinalizerImpl<Callable>;
		public:
			Finalizer(const Finalizer &) = delete;
			Finalizer(Finalizer &&) = default;
			
			template <typename CallableRef>
			explicit Finalizer(CallableRef &&_func)
				:base(std::forward<CallableRef>(_func)){}
			
			Finalizer &operator=(const Finalizer &) = delete;
			Finalizer &operator=(Finalizer &&) = default;
		};
		
		template <typename Callable>
		static inline Finalizer<Callable> makeFinalizer(Callable &&_func){
			return Finalizer<Callable>(std::forward<Callable>(_func));
		}
		
#define _BIG_NUM_CONTACT_(x, y) x ## y
#define _BIG_NUM_CONTACT2_(x, y) _BIG_NUM_CONTACT_(x, y)
#define _BIG_NUM_ADD_FINALIZER_(x) \
	auto _BIG_NUM_CONTACT2_(_finalizer, __LINE__) = bignum::_utility::makeFinalizer(x);
#define _BIG_NUM_ADD_FINALIZER_CODE_(x) _BIG_NUM_ADD_FINALIZER_([&](){x});
		
	}; // namepsace _utility
	
}; // namespace bignum
#endif // _BIG_NUM_MEMORY_HPP_