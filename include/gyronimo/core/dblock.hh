// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// @dblock.hh

#ifndef GYRONIMO_DBLOCK
#define GYRONIMO_DBLOCK

#include <ranges>

namespace gyronimo {

//! Interface to a *read-only* contiguous range of doubles.
/*!
    Useful to seamlessly interact with 3rd-party client code (c, fortran, etc.)
    that assume an array to be contiguously stored in memory and its name to
    hold a pointer to the first element. General `gyronimo` code should use only
    this abstract interface, actual translation is achived by deriving adapters
    for each specific case (c/fortran arrays, gsl_block, whatever).
*/
class dblock {
 public:
  virtual ~dblock() {};
  virtual const double* begin() const = 0;
  virtual const double* end() const = 0;
  const double* data() const {return this->begin();};
  size_t size() const {return this->end() - this->begin();};
};

//! Adapter for any contiguous STL range of doubles.
template<typename Range> requires
  std::ranges::contiguous_range<Range> &&
  std::same_as<std::ranges::range_value_t<Range>, double>
class dblock_adapter : public dblock {
 public:
  dblock_adapter(const Range& obj) : obj_ref_(obj) {};
  dblock_adapter(Range&& temp_obj)
      : moved_in_obj_(temp_obj), obj_ref_(moved_in_obj_) {};
  virtual const double* begin() const override {
      return std::ranges::data(obj_ref_);};
  virtual const double* end() const override {
      return std::ranges::data(obj_ref_) + std::ranges::size(obj_ref_);};
 private:
  const Range& obj_ref_;
  Range moved_in_obj_;
};

} // end namespace gyronimo.

#endif // GYRONIMO_DBLOCK
