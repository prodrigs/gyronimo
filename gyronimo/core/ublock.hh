// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2023 Paulo Rodrigues.

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

// @ublock.hh, this file is part of ::gyronimo::

#ifndef GYRONIMO_UBLOCK
#define GYRONIMO_UBLOCK

#include <execution>
#include <memory>

namespace gyronimo {

//! Continuous block of uninitialised memory to store objects of type `T`.
/*!
    Provides a convenient way to allocate a block of continous memory for
    objects that are not `DefaultConstructible` and thus not allowed as template
    arguments for common containers like `std::vector`. It would be possible to
    define an empty `std::vector<T>` and then construct and push individual
    objects, but this task could not be done in parallel. `ublock<T>` allows
    such parallel initialisation and provides also parallel destruction.
*/
template<typename T>
class ublock {
 public:
  ublock(size_t size)
      : size_(size),
        sptr_(std::allocator<T>().allocate(size), pdeleter(size)) {};
  T* data() { return sptr_.get(); };
  const T* data() const { return sptr_.get(); };
  size_t size() const { return size_; };
 private:
  const size_t size_;
  const std::shared_ptr<T> sptr_;
  struct pdeleter {
    const size_t size_;
    pdeleter(size_t size) : size_(size) {};
    void operator()(T* ptr) const {
      std::destroy(std::execution::par, ptr, ptr + size_);
      std::allocator<T>().deallocate(ptr, size_);
    };
  };
};

}  // end namespace gyronimo.

#endif  // GYRONIMO_UBLOCK
