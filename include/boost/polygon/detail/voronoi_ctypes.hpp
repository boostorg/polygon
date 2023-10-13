// Boost.Polygon library detail/voronoi_ctypes.hpp header file

//          Copyright Andrii Sydorchuk 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_POLYGON_DETAIL_VORONOI_CTYPES
#define BOOST_POLYGON_DETAIL_VORONOI_CTYPES

#include <boost/cstdint.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>

#if BOOST_VORONOI_USE_GMP
#include <gmp.h>
#if GMP_NAIL_BITS != 0
  #error boost::voronoi expects GMP_LIMB_BITS == 0
#endif
#if GMP_NUMB_BITS != 32 && GMP_NUMB_BITS != 64
  #error boost::voronoi expects GMP_NUMB_BITS == 32 or 64
#endif
#endif /* BOOST_VORONOI_USE_GMP */

namespace boost {
namespace polygon {
namespace detail {

typedef boost::int32_t int32;
typedef boost::int64_t int64;
typedef boost::uint32_t uint32;
typedef boost::uint64_t uint64;
typedef double fpt64;

// If two floating-point numbers in the same format are ordered (x < y),
// then they are ordered the same way when their bits are reinterpreted as
// sign-magnitude integers. Values are considered to be almost equal if
// their integer bits reinterpretations differ in not more than maxUlps units.
template <typename _fpt>
struct ulp_comparison;

template <>
struct ulp_comparison<fpt64> {
  enum Result {
    LESS = -1,
    EQUAL = 0,
    MORE = 1
  };

  Result operator()(fpt64 a, fpt64 b, unsigned int maxUlps) const {
    uint64 ll_a, ll_b;

    // Reinterpret double bits as 64-bit signed integer.
    std::memcpy(&ll_a, &a, sizeof(fpt64));
    std::memcpy(&ll_b, &b, sizeof(fpt64));

    // Positive 0.0 is integer zero. Negative 0.0 is 0x8000000000000000.
    // Map negative zero to an integer zero representation - making it
    // identical to positive zero - the smallest negative number is
    // represented by negative one, and downwards from there.
    if (ll_a < 0x8000000000000000ULL)
      ll_a = 0x8000000000000000ULL - ll_a;
    if (ll_b < 0x8000000000000000ULL)
      ll_b = 0x8000000000000000ULL - ll_b;

    // Compare 64-bit signed integer representations of input values.
    // Difference in 1 Ulp is equivalent to a relative error of between
    // 1/4,000,000,000,000,000 and 1/8,000,000,000,000,000.
    if (ll_a > ll_b)
      return (ll_a - ll_b <= maxUlps) ? EQUAL : LESS;
    return (ll_b - ll_a <= maxUlps) ? EQUAL : MORE;
  }
};

template <typename _fpt>
struct extened_exponent_fpt_traits;

template <>
struct extened_exponent_fpt_traits<fpt64> {
 public:
  typedef int exp_type;
  enum {
    MAX_SIGNIFICANT_EXP_DIF = 54
  };
};

// Floating point type wrapper. Allows to extend exponent boundaries to the
// integer type range. This class does not handle division by zero, subnormal
// numbers or NaNs.
template <typename _fpt, typename _traits = extened_exponent_fpt_traits<_fpt> >
class extended_exponent_fpt {
 public:
  typedef _fpt fpt_type;
  typedef typename _traits::exp_type exp_type;

  explicit extended_exponent_fpt(fpt_type val) {
    val_ = std::frexp(val, &exp_);
  }

  extended_exponent_fpt(fpt_type val, exp_type exp) {
    val_ = std::frexp(val, &exp_);
    exp_ += exp;
  }

  bool is_pos() const {
    return val_ > 0;
  }

  bool is_neg() const {
    return val_ < 0;
  }

  bool is_zero() const {
    return val_ == 0;
  }

  extended_exponent_fpt operator-() const {
    return extended_exponent_fpt(-val_, exp_);
  }

  extended_exponent_fpt operator+(const extended_exponent_fpt& that) const {
    if (this->val_ == 0.0 ||
        that.exp_ > this->exp_ + _traits::MAX_SIGNIFICANT_EXP_DIF) {
      return that;
    }
    if (that.val_ == 0.0 ||
        this->exp_ > that.exp_ + _traits::MAX_SIGNIFICANT_EXP_DIF) {
      return *this;
    }
    if (this->exp_ >= that.exp_) {
      exp_type exp_dif = this->exp_ - that.exp_;
      fpt_type val = std::ldexp(this->val_, exp_dif) + that.val_;
      return extended_exponent_fpt(val, that.exp_);
    } else {
      exp_type exp_dif = that.exp_ - this->exp_;
      fpt_type val = std::ldexp(that.val_, exp_dif) + this->val_;
      return extended_exponent_fpt(val, this->exp_);
    }
  }

  extended_exponent_fpt operator-(const extended_exponent_fpt& that) const {
    if (this->val_ == 0.0 ||
        that.exp_ > this->exp_ + _traits::MAX_SIGNIFICANT_EXP_DIF) {
      return extended_exponent_fpt(-that.val_, that.exp_);
    }
    if (that.val_ == 0.0 ||
        this->exp_ > that.exp_ + _traits::MAX_SIGNIFICANT_EXP_DIF) {
      return *this;
    }
    if (this->exp_ >= that.exp_) {
      exp_type exp_dif = this->exp_ - that.exp_;
      fpt_type val = std::ldexp(this->val_, exp_dif) - that.val_;
      return extended_exponent_fpt(val, that.exp_);
    } else {
      exp_type exp_dif = that.exp_ - this->exp_;
      fpt_type val = std::ldexp(-that.val_, exp_dif) + this->val_;
      return extended_exponent_fpt(val, this->exp_);
    }
  }

  extended_exponent_fpt operator*(const extended_exponent_fpt& that) const {
    fpt_type val = this->val_ * that.val_;
    exp_type exp = this->exp_ + that.exp_;
    return extended_exponent_fpt(val, exp);
  }

  extended_exponent_fpt operator/(const extended_exponent_fpt& that) const {
    fpt_type val = this->val_ / that.val_;
    exp_type exp = this->exp_ - that.exp_;
    return extended_exponent_fpt(val, exp);
  }

  extended_exponent_fpt& operator+=(const extended_exponent_fpt& that) {
    return *this = *this + that;
  }

  extended_exponent_fpt& operator-=(const extended_exponent_fpt& that) {
    return *this = *this - that;
  }

  extended_exponent_fpt& operator*=(const extended_exponent_fpt& that) {
    return *this = *this * that;
  }

  extended_exponent_fpt& operator/=(const extended_exponent_fpt& that) {
    return *this = *this / that;
  }

  extended_exponent_fpt sqr() const {
    return extended_exponent_fpt(this->val_ * this->val_, this->exp_ * 2);
  }

  extended_exponent_fpt sqrt() const {
    fpt_type val = val_;
    exp_type exp = exp_;
    if (exp & 1) {
      val *= 2.0;
      --exp;
    }
    return extended_exponent_fpt(std::sqrt(val), exp >> 1);
  }

  // Add to the exponent.
  extended_exponent_fpt addexp(int n) const {
    return extended_exponent_fpt(this->val_, this->exp_ + n);
  }

  fpt_type d() const {
    return std::ldexp(val_, exp_);
  }

 private:
  fpt_type val_;
  exp_type exp_;
};
typedef extended_exponent_fpt<double> efpt64;

template <typename _fpt>
extended_exponent_fpt<_fpt> sqr(const extended_exponent_fpt<_fpt>& that) {
  return that.sqr();
}

template <typename _fpt>
extended_exponent_fpt<_fpt> get_sqrt(const extended_exponent_fpt<_fpt>& that) {
  return that.sqrt();
}

template <typename _fpt>
bool is_pos(const extended_exponent_fpt<_fpt>& that) {
  return that.is_pos();
}

template <typename _fpt>
bool is_neg(const extended_exponent_fpt<_fpt>& that) {
  return that.is_neg();
}

template <typename _fpt>
bool is_zero(const extended_exponent_fpt<_fpt>& that) {
  return that.is_zero();
}

#if BOOST_VORONOI_USE_GMP
  #if GMP_NUMB_BITS == 32
    #define BOOST_VORONOI_64_T 0
  #elif GMP_NUMB_BITS == 64
    #define BOOST_VORONOI_64_T 1
  #else
    #error BOOST_VORONOI_USE_GMP: mp_limb_t is expected to be either 32 or 64bit integer type.
  #endif
#else
  #if INTPTR_MAX == INT64_MAX
    // Multi-precision using 64bit ints on a 64bit compiler.
    #define BOOST_VORONOI_64_T 1
  #else
    // Multi-precision using 32bit ints otherwise.
    #define BOOST_VORONOI_64_T 0
  #endif
#endif /* BOOST_VORONOI_USE_GMP */

// Very efficient stack allocated big integer class.
// Supports next set of arithmetic operations: +, -, *.
template<std::size_t AN>
class extended_int {
 public:

  static constexpr std::size_t N = AN;

  using chunk_type  = 
#if BOOST_VORONOI_USE_GMP
    mp_limb_t
#elif BOOST_VORONOI_64_T
    uint64
#else
    uint32
#endif
    ;

  using chunk_type2 = 
#if BOOST_VORONOI_64_T
    #ifdef _MSC_VER
        // Visual studio does not support 128bit type natively. boost::multiprecision::uint128_t utilizes 32bit limbs for multiplication,
        // which is not the best we can do on a 64bit system. Give the caller an option to replace the 128bit uint with its own optimized solution.
        #ifdef BOOST_VORONOI_UINT128T
            BOOST_VORONOI_UINT128T
        #else
            boost::multiprecision::uint128_t
        #endif
    #else
        unsigned __int128
    #endif
#else
    uint64
#endif
    ;

  static constexpr int chunk_shift =
#if BOOST_VORONOI_USE_GMP
    GMP_NUMB_BITS
#elif BOOST_VORONOI_64_T
    64
#else
    32
#endif
    ;

  // Exponent conversion between this type and the extended floating point type.
  static constexpr int chunks_to_fptexp =
#if BOOST_VORONOI_64_T
    10
#else
    5
#endif
    ;

  extended_int() {}

  extended_int(int32 that) {
    if (that > 0) {
      this->chunks_[0] = that;
      this->count_ = 1;
    } else if (that < 0) {
      this->chunks_[0] = -that;
      this->count_ = -1;
    } else {
      this->count_ = 0;
    }
  }

  extended_int(int64 that) {
    if (that > 0) {
#if BOOST_VORONOI_64_T
      this->chunks_[0] = that;
      this->count_ = 1;
#else
      this->chunks_[0] = static_cast<chunk_type>(that);
      this->chunks_[1] = that >> 32;
      this->count_ = this->chunks_[1] ? 2 : 1;
#endif
    } else if (that < 0) {
#if BOOST_VORONOI_64_T
      this->chunks_[0] = -that;
      this->count_ = -1;
#else
      that = -that;
      this->chunks_[0] = static_cast<chunk_type>(that);
      this->chunks_[1] = that >> 32;
      this->count_ = this->chunks_[1] ? -2 : -1;
#endif
    } else {
      this->count_ = 0;
    }
  }

  template<std::size_t M>
  extended_int(const extended_int<M>& that) {
    this->count_ = that.count();
    std::memcpy(this->chunks_, that.chunks(), that.size() * sizeof(chunk_type));
  }

  extended_int& operator=(int32 that) {
    if (that > 0) {
      this->chunks_[0] = that;
      this->count_ = 1;
    } else if (that < 0) {
      this->chunks_[0] = -that;
      this->count_ = -1;
    } else {
      this->count_ = 0;
    }
    return *this;
  }

  extended_int& operator=(int64 that) {
    if (that > 0) {
#if BOOST_VORONOI_64_T
      this->chunks_[0] = that;
      this->count_ = 1;
#else
      this->chunks_[0] = static_cast<chunk_type>(that);
      this->chunks_[1] = that >> 32;
      this->count_ = this->chunks_[1] ? 2 : 1;
#endif
    } else if (that < 0) {
#if BOOST_VORONOI_64_T
      this->chunks_[0] = -that;
      this->count_ = -1;
#else
      that = -that;
      this->chunks_[0] = static_cast<chunk_type>(that);
      this->chunks_[1] = that >> 32;
      this->count_ = this->chunks_[1] ? -2 : -1;
#endif
    } else {
      this->count_ = 0;
    }
    return *this;
  }

  template<std::size_t M>
  extended_int& operator=(const extended_int<M>& that) {
    this->count_ = that.count();
    std::memcpy(this->chunks_, that.chunks(), that.size() * sizeof(chunk_type));
    return *this;
  }

  bool is_pos() const {
    return this->count_ > 0;
  }

  bool is_neg() const {
    return this->count_ < 0;
  }

  bool is_zero() const {
    return this->count_ == 0;
  }

  bool operator==(const extended_int& that) const {
    if (this->count_ != that.count())
      return false;
    for (std::size_t i = 0; i < this->size(); ++i)
      if (this->chunks_[i] != that.chunks()[i])
        return false;
    return true;
  }

  bool operator!=(const extended_int& that) const {
    return !(*this == that);
  }

  bool operator<(const extended_int& that) const {
    if (this->count_ != that.count())
      return this->count_ < that.count();
    std::size_t i = this->size();
    if (!i)
      return false;
    do {
      --i;
      if (this->chunks_[i] != that.chunks()[i])
        return (this->chunks_[i] < that.chunks()[i]) ^ (this->count_ < 0);
    } while (i);
    return false;
  }

  bool operator>(const extended_int& that) const {
    return that < *this;
  }

  bool operator<=(const extended_int& that) const {
    return !(that < *this);
  }

  bool operator>=(const extended_int& that) const {
    return !(*this < that);
  }

  extended_int operator-() const {
    extended_int ret_val = *this;
    ret_val.neg();
    return ret_val;
  }

  void neg() {
    this->count_ = -this->count_;
  }

  extended_int operator+(const extended_int& that) const {
    extended_int ret_val;
    ret_val.add(*this, that);
    return ret_val;
  }

  void add(const extended_int& e1, const extended_int& e2) {
    if (!e1.count()) {
      *this = e2;
      return;
    }
    if (!e2.count()) {
      *this = e1;
      return;
    }
    if ((e1.count() > 0) ^ (e2.count() > 0)) {
      dif(e1.chunks(), e1.size(), e2.chunks(), e2.size());
    } else {
      add(e1.chunks(), e1.size(), e2.chunks(), e2.size());
    }
    if (e1.count() < 0)
      this->count_ = -this->count_;
  }

  extended_int operator-(const extended_int& that) const {
    extended_int ret_val;
    ret_val.dif(*this, that);
    return ret_val;
  }

  void dif(const extended_int& e1, const extended_int& e2) {
    if (!e1.count()) {
      *this = e2;
      this->count_ = -this->count_;
      return;
    }
    if (!e2.count()) {
      *this = e1;
      return;
    }
    if ((e1.count() > 0) ^ (e2.count() > 0)) {
      add(e1.chunks(), e1.size(), e2.chunks(), e2.size());
    } else {
      dif(e1.chunks(), e1.size(), e2.chunks(), e2.size());
    }
    if (e1.count() < 0)
      this->count_ = -this->count_;
  }

  extended_int operator*(int32 that) const {
    extended_int temp(that);
    return (*this) * temp;
  }

  extended_int operator*(int64 that) const {
    extended_int temp(that);
    return (*this) * temp;
  }

  extended_int operator*(const extended_int& that) const {
    extended_int ret_val;
    ret_val.mul(*this, that);
    return ret_val;
  }

  extended_int sqr() const {
    extended_int ret_val;
    if (! this->count_)
      ret_val.count_ = 0;
    else
      ret_val.mksqr(this->chunks_, this->size());
    return ret_val;
  }

  // Square of a positive non-zero value with more than half the limbs utilized.
  // Such a number does not fit the result, thus it is right shifted.
  // The caller shall be aware of the right shift performed by this function
  // on the result.
  extended_int sqrext() const {
    extended_int ret_val;
    assert(this->size() * 2 > N);
    ret_val.mksqrext(this->chunks_, this->size());
    return ret_val;
  }

  void mul(const extended_int& e1, const extended_int& e2) {
    if (!e1.count() || !e2.count()) {
      this->count_ = 0;
      return;
    }
    mul(e1.chunks(), e1.size(), e2.chunks(), e2.size());
    if ((e1.count() > 0) ^ (e2.count() > 0))
      this->count_ = -this->count_;
  }

  // Multiple of two nonzero values, result of which does not fit the resulting type,
  // thus it is right shifted.
  // The caller shall be aware of the right shift performed by this function
  // on the result.
  extended_int mulext(const extended_int& rhs) const {
    assert(this->size() > 0);
    assert(rhs.size() > 0);
    assert(this->size() + rhs.size() > N);
    extended_int ret_val;
    ret_val.mulext_(this->chunks(), this->size(), rhs.chunks(), rhs.size());
    if ((this->count() > 0) ^ (rhs.count() > 0))
      ret_val.count_ = -ret_val.count_;
    return ret_val;
  }

  // Shift the limbs of this object by nshift positions right.
  void rshiftme(std::size_t nshift) {
    std::size_t m = this->size();
    if (nshift >= m) {
      // Underflow to zero.
      this->count_ = 0;
    } else {
      m -= nshift;
      for (std::size_t i = 0; i < m; ++ i)
        this->chunks_[i] = this->chunks_[i + nshift];
      this->count_ = this->count_ > 0 ? m : - m;
    }
  }

  const chunk_type* chunks() const {
    return chunks_;
  }

  int32 count() const {
    return count_;
  }

  std::size_t size() const {
    return (std::abs)(count_);
  }

  std::pair<fpt64, int> p() const {
    std::pair<fpt64, int> ret_val(0, 0);
    std::size_t sz = this->size();
    if (!sz) {
      return ret_val;
    } else {
#if BOOST_VORONOI_64_T
      if (sz == 1) {
        auto   l0 = this->chunks_[0];
        uint32 l, h;
        memcpy(&l, &l0, sizeof(uint32));
        memcpy(&h, reinterpret_cast<const char*>(&l0) + sizeof(uint32), sizeof(uint32));
        if (h == 0) {
          ret_val.first = static_cast<fpt64>(static_cast<uint32>(l));
        } else {
          ret_val.first = static_cast<fpt64>(h) *
                          static_cast<fpt64>(0x100000000LL) +
                          static_cast<fpt64>(l);
        }
      } else {
        assert(sz >= 2);
        auto last = this->chunks_[sz - 1];
        auto prev = this->chunks_[sz - 2];
        uint32 l, h;
        memcpy(&l, &last, sizeof(uint32));
        memcpy(&h, reinterpret_cast<const char*>(&last) + sizeof(uint32), sizeof(uint32));
        sz *= 2;
        if (h == 0) {
          -- sz;
          ret_val.first  = static_cast<fpt64>(l);
          memcpy(&l, &prev, sizeof(uint32));
          memcpy(&h, reinterpret_cast<const char*>(&prev) + sizeof(uint32), sizeof(uint32));          
          ret_val.first *= static_cast<fpt64>(0x100000000LL);
          ret_val.first += static_cast<fpt64>(h);
          ret_val.first *= static_cast<fpt64>(0x100000000LL);
          ret_val.first += static_cast<fpt64>(l);
        } else {
          ret_val.first  = static_cast<fpt64>(h);
          ret_val.first *= static_cast<fpt64>(0x100000000LL);
          ret_val.first += static_cast<fpt64>(l);
          memcpy(&h, reinterpret_cast<const char*>(&prev) + sizeof(uint32), sizeof(uint32));
          ret_val.first *= static_cast<fpt64>(0x100000000LL);
          ret_val.first += static_cast<fpt64>(h);
        }
        ret_val.second = static_cast<int>((sz - 3) << (chunks_to_fptexp / 2));
      }
#else
      if (sz == 1) {
        ret_val.first = static_cast<fpt64>(this->chunks_[0]);
      } else if (sz == 2) {
        ret_val.first = static_cast<fpt64>(this->chunks_[1]) *
                        static_cast<fpt64>(0x100000000LL) +
                        static_cast<fpt64>(this->chunks_[0]);
      } else {
        for (std::size_t i = 1; i <= 3; ++i) {
          ret_val.first *= static_cast<fpt64>(0x100000000LL);
          ret_val.first += static_cast<fpt64>(this->chunks_[sz - i]);
        }
        ret_val.second = static_cast<int>((sz - 3) << chunks_to_fptexp);
      }
#endif
    }
    if (this->count_ < 0)
      ret_val.first = -ret_val.first;
    return ret_val;
  }

  fpt64 d() const {
    std::pair<fpt64, int> p = this->p();
    return std::ldexp(p.first, p.second);
  }

 private:

  void add(const chunk_type* c1, std::size_t sz1,
           const chunk_type* c2, std::size_t sz2) {
    if (sz1 < sz2) {
      add(c2, sz2, c1, sz1);
      return;
    }

    assert(sz1 >= sz2);
    assert(sz2 > 0);
    this->count_ = static_cast<int32>(sz1);
#if BOOST_VORONOI_USE_GMP
    mp_limb_t temp = mpn_add(this->chunks_, c1, sz1, c2, sz2);
#else
    chunk_type2 temp = 0;
    for (std::size_t i = 0; i < sz2; ++i) {
      temp += static_cast<chunk_type2>(c1[i]) + static_cast<chunk_type2>(c2[i]);
      this->chunks_[i] = static_cast<chunk_type>(temp);
      temp >>= chunk_shift;
    }
    for (std::size_t i = sz2; i < sz1; ++i) {
      temp += static_cast<chunk_type2>(c1[i]);
      this->chunks_[i] = static_cast<chunk_type>(temp);
      temp >>= chunk_shift;
    }
#endif
    // One should never trim the most significant limb.
    assert(temp == 0 || this->count_ < N);
    if (temp && (this->count_ != N)) {
      this->chunks_[this->count_] = static_cast<chunk_type>(temp);
      ++this->count_;
    }
  }

  void dif(const chunk_type* c1, std::size_t sz1,
           const chunk_type* c2, std::size_t sz2,
           bool rec = false) {
    if (sz1 < sz2) {
      dif(c2, sz2, c1, sz1, true);
      this->count_ = -this->count_;
      return;
    }

    if ((sz1 == sz2) && !rec) {
      do {
        --sz1;
        if (c1[sz1] < c2[sz1]) {
          ++sz1;
          dif(c2, sz1, c1, sz1, true);
          this->count_ = -this->count_;
          return;
        } else if (c1[sz1] > c2[sz1]) {
          ++sz1;
          break;
        }
      } while (sz1);
      if (!sz1) {
        this->count_ = 0;
        return;
      }
      sz2 = sz1;
    }
#if BOOST_VORONOI_USE_GMP
    assert(sz1 > sz2 || (sz1 == sz2 && mpn_cmp(c1, c2, sz1) >= 0));
    mp_limb_t borrow = mpn_sub(this->chunks_, c1, sz1, c2, sz2);
    assert(! borrow);
    this->count_ = sz1;
    if (! this->chunks_[this->count_ - 1])
      -- this->count_;
#else
    this->count_ = static_cast<int32>(sz1-1);
    bool flag = false;
    for (std::size_t i = 0; i < sz2; ++i) {
      this->chunks_[i] = c1[i] - c2[i] - (flag?1:0);
      flag = (c1[i] < c2[i]) || ((c1[i] == c2[i]) && flag);
    }
    for (std::size_t i = sz2; i < sz1; ++i) {
      this->chunks_[i] = c1[i] - (flag?1:0);
      flag = !c1[i] && flag;
    }
    if (this->chunks_[this->count_])
      ++this->count_;
#endif
  }

  // multiple of two positive nonzero values
  void mul(const chunk_type* c1, std::size_t sz1,
           const chunk_type* c2, std::size_t sz2) {
    assert(sz1 > 0);
    assert(sz2 > 0);
    assert(this->chunks_ != c1);
    assert(this->chunks_ != c2);
    assert(sz1 + sz2 <= N);
#if BOOST_VORONOI_USE_GMP
    this->count_ = sz1 + sz2;
    auto msb = (sz1 >= sz2) ? mpn_mul(this->chunks_, c1, sz1, c2, sz2) : mpn_mul(this->chunks_, c2, sz2, c1, sz1);
    assert(msb == this->chunks_[this->count_ - 1]);
    if (msb == 0)
      -- this->count_;
#else
    chunk_type2 cur = 0, nxt, tmp;
    // This stinks. One only calculates maximum N of the most significant limbs
    // without letting the caller know.
    this->count_ = sz1 + sz2 - 1;
    for (std::size_t shift = 0; shift < static_cast<std::size_t>(this->count_); ++shift) {
      nxt = 0;
      for (std::size_t first = 0; first <= shift; ++first) {
        if (first >= sz1)
          break;
        std::size_t second = shift - first;
        if (second >= sz2)
          continue;
        tmp = static_cast<chunk_type2>(c1[first]) * static_cast<chunk_type2>(c2[second]);
        cur += static_cast<chunk_type>(tmp);
        nxt += tmp >> chunk_shift;
      }
      this->chunks_[shift] = static_cast<chunk_type>(cur);
      cur = nxt + (cur >> chunk_shift);
    }
    // One should never trim the most significant limb.
    assert(cur == 0 || this->count_ < N);
    if (cur && (this->count_ != N)) {
      this->chunks_[this->count_] = static_cast<chunk_type>(cur);
      ++this->count_;
    }
#endif
  }

  // Multiple of two positive nonzero values, result of which does not fit the resulting type,
  // thus it is right shifted.
  void mulext_(const chunk_type* c1, std::size_t sz1,
               const chunk_type* c2, std::size_t sz2) {
    assert(sz1 > 0 && sz2 > 0 && sz1 + sz2 > N);
    chunk_type2 cur = 0, nxt, tmp;
    // This stinks. One only calculates maximum N of the most significant limbs
    // without letting the caller know.
    this->count_ = N-1;
    for (std::size_t shift = 0; shift < static_cast<std::size_t>(this->count_); ++shift) {
      nxt = 0;
      for (std::size_t first = 0; first <= shift; ++first) {
        if (first >= sz1)
          break;
        std::size_t second = shift - first;
        if (second >= sz2)
          continue;
        tmp = static_cast<chunk_type2>(c1[first]) * static_cast<chunk_type2>(c2[second]);
        cur += static_cast<chunk_type>(tmp);
        nxt += tmp >> chunk_shift;
      }
      this->chunks_[shift] = static_cast<chunk_type>(cur);
      cur = nxt + (cur >> chunk_shift);
    }
    // One should never trim the most significant limb.
    assert(cur == 0 || this->count_ < N);
    if (cur && (this->count_ != N)) {
      this->chunks_[this->count_] = static_cast<chunk_type>(cur);
      ++this->count_;
    }
  }

  // square of a positive non-zero value
  void mksqr(const chunk_type* c, std::size_t sz) {
    assert(2 * sz <= N);
#if BOOST_VORONOI_USE_GMP
    assert(sz > 0);
    assert(this->chunks_ != c);
    this->count_ = sz * 2;
    mpn_sqr(this->chunks_, c, sz);
    if (this->chunks_[this->count_ - 1] == 0)
      -- this->count_;
#else
    chunk_type2 cur = 0, nxt, tmp;
    this->count_ = sz * 2 - 1;
    for (std::size_t shift = 0; shift < static_cast<std::size_t>(this->count_); ++shift) {
      nxt = 0;
      for (std::size_t first = 0; first <= shift; ++first) {
        if (first >= sz)
          break;
        std::size_t second = shift - first;
        if (second >= sz)
          continue;
        tmp = static_cast<chunk_type2>(c[first]) * static_cast<chunk_type2>(c[second]);
        cur += static_cast<chunk_type>(tmp);
        nxt += tmp >> chunk_shift;
      }
      this->chunks_[shift] = static_cast<chunk_type>(cur);
      cur = nxt + (cur >> chunk_shift);
    }
    // One should never trim the most significant limb.
    assert(cur == 0 || this->count_ < N);
    if (cur && (this->count_ != N)) {
      this->chunks_[this->count_] = static_cast<chunk_type>(cur);
      ++this->count_;
    }
#endif
  }

  // Square of a positive non-zero value with more than half the limbs utilized.
  // Such a number does not fit the result, thus it is right shifted.
  void mksqrext(const chunk_type* c, std::size_t sz) {
    chunk_type2 cur = 0, nxt, tmp;
    this->count_ = N-1;
    for (std::size_t shift = 0; shift < static_cast<std::size_t>(this->count_); ++shift) {
      nxt = 0;
      for (std::size_t first = 0; first <= shift; ++first) {
        if (first >= sz)
          break;
        std::size_t second = shift - first;
        if (second >= sz)
          continue;
        tmp = static_cast<chunk_type2>(c[first]) * static_cast<chunk_type2>(c[second]);
        cur += static_cast<chunk_type>(tmp);
        nxt += tmp >> chunk_shift;
      }
      this->chunks_[shift] = static_cast<chunk_type>(cur);
      cur = nxt + (cur >> chunk_shift);
    }
    // One should never trim the most significant limb.
    assert(cur == 0 || this->count_ < N);
    if (cur && (this->count_ != N)) {
      this->chunks_[this->count_] = static_cast<chunk_type>(cur);
      ++this->count_;
    }
  }

  chunk_type chunks_[N];
  int32 count_;
};

template <std::size_t N>
extended_int<N> sqr(const extended_int<N>& that) {
  return that.sqr();
}

template <std::size_t N>
bool is_pos(const extended_int<N>& that) {
  return that.is_pos();
}

template <std::size_t N>
bool is_neg(const extended_int<N>& that) {
  return that.is_neg();
}

template <std::size_t N>
bool is_zero(const extended_int<N>& that) {
  return that.is_zero();
}

struct type_converter_fpt {
  template <typename T>
  fpt64 operator()(const T& that) const {
    return static_cast<fpt64>(that);
  }

  template <std::size_t N>
  fpt64 operator()(const extended_int<N>& that) const {
    return that.d();
  }

  fpt64 operator()(const extended_exponent_fpt<fpt64>& that) const {
    return that.d();
  }
};

struct type_converter_efpt {
  template <std::size_t N>
  extended_exponent_fpt<fpt64> operator()(const extended_int<N>& that) const {
    std::pair<fpt64, int> p = that.p();
    return extended_exponent_fpt<fpt64>(p.first, p.second);
  }
};

// Voronoi coordinate type traits make it possible to extend algorithm
// input coordinate range to any user provided integer type and algorithm
// output coordinate range to any ieee-754 like floating point type.
template <typename T>
struct voronoi_ctype_traits;

template <>
struct voronoi_ctype_traits<int32> {
  typedef int32 int_type;
  typedef int64 int_x2_type;
  typedef uint64 uint_x2_type;
#if BOOST_VORONOI_64_T
  // using uint64 
  typedef extended_int<33> big_int_type;
#else
  // using uint32
  typedef extended_int<64> big_int_type;
#endif
  typedef fpt64 fpt_type;
  typedef extended_exponent_fpt<fpt_type> efpt_type;
  typedef ulp_comparison<fpt_type> ulp_cmp_type;
  typedef type_converter_fpt to_fpt_converter_type;
  typedef type_converter_efpt to_efpt_converter_type;
};
}  // detail
}  // polygon
}  // boost

#endif  // BOOST_POLYGON_DETAIL_VORONOI_CTYPES
