// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// ::gyronimo:: is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ::gyronimo:: is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with ::gyronimo::.  If not, see <https://www.gnu.org/licenses/>.

// @dblock.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_DBLOCK
#define GYRONIMO_DBLOCK

#include <ranges>

namespace gyronimo {

//! Interface to a **read-only** contiguous range of doubles.
/*!
    Useful to seamlessly interact with 3rd-party client code (c, fortran, etc.)
    or libraries that provide their own data structures to store arrays of
    doubles contiguously in memory. General-purpose `gyronimo` code **should use
    only** this abstract interface. The actual translation to each particular
    data type is achieved by deriving adapter classes for each specific case
    (c/fortran arrays, gsl_block, whatever). The class supports std::ranges
    semantics.
*/
class dblock {
 public:
  typedef size_t size_type;
  typedef double value_type;
  typedef const double* iterator;
  typedef const double& reference;
  typedef const double* const_iterator;
  typedef const double& const_reference;
  virtual ~dblock() {};
  virtual iterator begin() const = 0;
  virtual iterator end() const = 0;
  iterator data() const {return this->begin();};
  size_type size() const {return this->end() - this->begin();};
  value_type front() const {return this->data()[0];};
  value_type back() const {return this->data()[this->size() - 1];};
  value_type operator[](size_t k) const {return this->data()[k];};
};

//! Templated adapter for any contiguous STL range of doubles.
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
