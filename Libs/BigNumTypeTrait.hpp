#ifndef _BIG_NUM_HPP_
#error "This header must be included through BigNum.hpp"
#endif // _BIG_NUM_HPP_

#ifndef _BIG_NUM_TYPE_TRAIT_HPP_
#define _BIG_NUM_TYPE_TRAIT_HPP_

#include <limits>
#include <cstdint>
#include <type_traits>
#include <initializer_list>

namespace bignum{
	
	namespace _type{
		
		template <typename T>
		struct isSigned
			:public std::integral_constant<bool, std::numeric_limits<T>::is_signed>{};
		
		template <typename T, typename U>
		struct isRLRef
			:public std::is_same<T, typename std::remove_cv<typename std::remove_reference<U>::type>::type>{};
		
		/*template <typename, typename>
		struct isRLRef{
			static constexpr bool value = false;
		};
		template <typename T>
		struct isRLRef<T, const T &>{
			static constexpr bool value = true;
		};
		template <typename T>
		struct isRLRef<T, T &>{
			static constexpr bool value = true;
		};
		template <typename T>
		struct isRLRef<T, const T &&>{
			static constexpr bool value - true;
		};
		template <typename T>
		struct isRLRef<T, T &&>{
			static constexpr bool value = true;
		};*/
		
		template <typename T>
		struct squareType{};
		template <>
		struct squareType<int8_t>{
			using type = int16_t;
		};
		template <>
		struct squareType<int16_t>{
			using type = int32_t;
		};
		template <>
		struct squareType<int32_t>{
			using type = int64_t;
		};
		template <>
		struct squareType<uint8_t>{
			using type = uint16_t;
		};
		template <>
		struct squareType<uint16_t>{
			using type = uint32_t;
		};
		template <>
		struct squareType<uint32_t>{
			using type = uint64_t;
		};
		
		template <class T1, class T2, typename SFINAE1 = void, typename SFINAE2 = void>
		struct conj2;
		
		template <class T1, class T2>
		struct conj2<T1, T2, 
			typename std::enable_if<std::is_base_of<std::true_type, T1>::value>::type,
			typename std::enable_if<std::is_base_of<std::true_type, T2>::value>::type>
			:public std::true_type{};
		
		template <class T1, class T2>
		struct conj2<T1, T2, 
			typename std::enable_if<std::is_base_of<std::true_type, T1>::value>::type,
			typename std::enable_if<std::is_base_of<std::false_type, T2>::value>::type>
			:public std::false_type{};
		
		template <class T1, class T2>
		struct conj2<T1, T2, 
			typename std::enable_if<std::is_base_of<std::false_type, T1>::value>::type,
			typename std::enable_if<std::is_base_of<std::true_type, T2>::value>::type>
			:public std::false_type{};
		
		template <class T1, class T2>
		struct conj2<T1, T2, 
			typename std::enable_if<std::is_base_of<std::false_type, T1>::value>::type,
			typename std::enable_if<std::is_base_of<std::false_type, T2>::value>::type>
			:public std::false_type{};
		
		template <class... Ts>
		struct conj;
		
		template <class T>
		struct conj<T>:public std::is_base_of<std::true_type, T>::type{};
		
		template <class T1, class T2, class... Ts>
		struct conj<T1, T2, Ts...>:public conj2<T1, conj<T2, Ts...>>{};
		
		template <typename Int, Int M, Int N, typename G, typename E, typename L>
		struct CompareCond{
			using type = typename std::conditional<(M > N), G, 
				typename std::conditional<M == N, E, L>::type>::type;
		};
		
		template <typename ...Args>
		struct StaticList{};
		
		struct NothingType;
		
		template <typename T>
		struct SFINAEWrapper{
			using type = T;
		};
		template<>
		struct SFINAEWrapper<NothingType>{};
		
		template <typename, class>
		struct PushFront;
		template <typename Arg, typename ...Args>
		struct PushFront<Arg, StaticList<Args...>>{
			using type = StaticList<Arg, Args...>;
		};
		template <typename Arg>
		struct PushFront<Arg, StaticList<>>{
			using type = StaticList<Arg>;
		};
		template <typename Arg>
		struct PushFront<Arg, NothingType>{
			using type = NothingType;
		};
		
		template <template <typename> class, class>
		struct Filter;
		template <template <typename> class Predicate, typename Arg, typename ...Args>
		struct Filter<Predicate, StaticList<Arg, Args...>>{
			using type = typename std::conditional<Predicate<Arg>::value, 
				typename PushFront<Arg, typename Filter<Predicate, StaticList<Args...>>::type>::type, 
				typename Filter<Predicate, StaticList<Args...>>::type>::type;
		};
		template <template <typename> class Predicate>
		struct Filter<Predicate, StaticList<>>{
			using type = StaticList<>;
		};
		
		// SFINAE-friendly Check
		template <template <typename> class, class>
		struct CheckImpl;
		template <template <typename> class Predicate, typename Arg, typename ...Args>
		struct CheckImpl<Predicate, StaticList<Arg, Args...>>{
			using type = typename std::conditional<Predicate<Arg>::value, 
				typename PushFront<Arg, typename CheckImpl<Predicate, StaticList<Args...>>::type>::type, 
				NothingType>::type;
		};
		template <template <typename> class Predicate>
		struct CheckImpl<Predicate, StaticList<>>{
			using type = StaticList<>;
		};
		
		template <template <typename> class Predicate, class List>
		struct Check:public SFINAEWrapper<typename CheckImpl<Predicate, List>::type>{};
		
		template <template <typename> class, class>
		struct Map;
		template <template <typename> class Predicate, typename Arg, typename ...Args>
		struct Map<Predicate, StaticList<Arg, Args...>>{
			using type = typename PushFront<typename Predicate<Arg>::type, 
				typename Map<Predicate, StaticList<Args...>>::type>::type;
		};
		template <template <typename> class Predicate>
		struct Map<Predicate, StaticList<>>{
			using type = StaticList<>;
		};
		
		template <typename T>
		class PackExpandHelper{
		public:
			constexpr explicit PackExpandHelper(std::initializer_list<T>){}
		};
		
	}; // namespace bignum::_type
	
}; // namespace bignum
#endif // _BIG_NUM_TYPE_TRAIT_HPP_