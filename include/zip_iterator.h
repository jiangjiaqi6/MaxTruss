#pragma once

#include <iterator>
#include <limits>
#include <tuple>
#include <type_traits>

template <size_t N>
struct tuple_util {
  template <typename TupleType, typename DifferenceType>
  static void increment(TupleType& it, DifferenceType forward) {
    std::get<N - 1>(it) += forward;
    tuple_util<N - 1>::increment(it, forward);
  }
  template <typename TupleType, typename DifferenceType>
  static bool check_sync(const TupleType& it1, const TupleType& it2, DifferenceType val) {
    if (std::get<N - 1>(it1) - std::get<N - 1>(it2) != val) return false;
    return tuple_util<N - 1>::check_sync(it1, it2, val);
  }
};

template <>
struct tuple_util<0> {
  template <typename TupleType, typename DifferenceType>
  static void increment(TupleType&, DifferenceType) {}
  template <typename TupleType, typename DifferenceType>
  static bool check_sync(const TupleType&, const TupleType&, DifferenceType) {
    return true;
  }
};

template <typename TupleReturnType>
struct make_references {
  template <typename TupleType, std::size_t... Is>
  TupleReturnType operator()(const TupleType& t, std::index_sequence<Is...>) {
    return std::tie(*std::get<Is>(t)...);
  }
};

// A simple wrapper over a tuple of references.
// The class is designed to hold a temporary tuple of reference
// after dereferencing a zip_iterator; in particular, it is needed
// to swap these rvalue tuples. Any other usage is not supported.
template <typename... T>
struct tuplewrapper : public std::tuple<typename std::enable_if<std::is_reference<T>::value, T&&>::type...> {
  // In the context of this class, T is a reference, so T&& is a "forwarding reference"
  typedef std::tuple<T&&...> base_type;
  // Construct from the result of std::tie
  tuplewrapper(const base_type& in) : base_type(in) {}

  // Assign any tuple convertible to std::tuple<T&&...>: *it = a_tuple;
  template <typename... U>
  tuplewrapper& operator=(const std::tuple<U...>& other) {
    base_type::operator=(other);
    return *this;
  }
  // Swap rvalue tuples: swap(*it1,*it2);
  friend void swap(tuplewrapper&& a, tuplewrapper&& b) { std::swap<T&&...>(a, b); }
};

template <typename... Types>
class zip_iterator {
  static const std::size_t num_types = sizeof...(Types);
  typedef std::tuple<Types...> it_types;

 public:
  typedef typename std::make_signed<std::size_t>::type difference_type;
  typedef std::tuple<typename std::iterator_traits<Types>::value_type...> value_type;
  typedef tuplewrapper<typename std::iterator_traits<Types>::reference...> reference;
  typedef std::tuple<typename std::iterator_traits<Types>::pointer...> pointer;
  typedef std::random_access_iterator_tag iterator_category;

  zip_iterator() : my_it() {}
  explicit zip_iterator(Types... args) : my_it(std::make_tuple(args...)) {}
  zip_iterator(const zip_iterator& input) : my_it(input.my_it) {}
  zip_iterator& operator=(const zip_iterator& input) {
    my_it = input.my_it;
    return *this;
  }

  reference operator*() const { return make_references<reference>()(my_it, std::make_index_sequence<num_types>()); }
  reference operator[](difference_type i) const { return *(*this + i); }

  difference_type operator-(const zip_iterator& it) const { return std::get<0>(my_it) - std::get<0>(it.my_it); }

  zip_iterator& operator+=(difference_type forward) {
    tuple_util<num_types>::increment(my_it, forward);
    return *this;
  }
  zip_iterator& operator-=(difference_type backward) { return *this += -backward; }
  zip_iterator& operator++() { return *this += 1; }
  zip_iterator& operator--() { return *this -= 1; }

  zip_iterator operator++(int) {
    zip_iterator it(*this);
    ++(*this);
    return it;
  }
  zip_iterator operator--(int) {
    zip_iterator it(*this);
    --(*this);
    return it;
  }

  zip_iterator operator-(difference_type backward) const {
    zip_iterator it(*this);
    return it -= backward;
  }
  zip_iterator operator+(difference_type forward) const {
    zip_iterator it(*this);
    return it += forward;
  }
  friend zip_iterator operator+(difference_type forward, const zip_iterator& it) { return it + forward; }

  bool operator==(const zip_iterator& it) const { return *this - it == 0; }
  it_types base() const { return my_it; }

  bool operator!=(const zip_iterator& it) const { return !(*this == it); }
  bool operator<(const zip_iterator& it) const { return *this - it < 0; }
  bool operator>(const zip_iterator& it) const { return it < *this; }
  bool operator<=(const zip_iterator& it) const { return !(*this > it); }
  bool operator>=(const zip_iterator& it) const { return !(*this < it); }

 private:
  it_types my_it;
};

template <typename... T>
zip_iterator<T...> make_zip_iterator(T... args) {
  return zip_iterator<T...>(args...);
}
