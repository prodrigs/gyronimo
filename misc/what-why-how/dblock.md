gyronimo::dblock
================

What is it?
-----------

A `dblock` is an abstraction of the *access* interface to a very simple
kind of container: a contiguous collection (i.e., a block) of `double`
numbers. In short, it provides access to the block's begin and end,
together to its size. It is not supposed to work as a container
structure by itself (i.e., supplying the usual container functionality),
as there are many good alternatives to that purpose. Instead, it is to
be employed as a 'lingua franca' when working with disparate types of
double containers, in order to facilitate the exchange of the contained
data and the writing of general, container-agnostic code.

Why is it needed for?
---------------------

Many c++ third-party libraries with which `gyronimo` interacts have
their own implementation of the concept "a vector/array of doubles".
Possible candidates are, for instance, the ubiquitous STL containers
`std::array<double>`, `std::vector<double>`, and
`std::valarray<double>`. Besides these, the
[GSL](https://www.gnu.org/software/gsl) library adds the type
`gsl_block` to the party, while [Eigen](https://eigen.tuxfamily.org)
provides the `MatrixVd` and `VectorXd` types. As one would easily
expect, such libraries have many of their routine arguments tied to the
particular kind of types they happen to define, e.g.,

```
gsl_block A, B;
...
double C = gsl_my_favourite_function(A, B)
```

Such proliferation of particular types creates a problem for `gyronimo`
developers that strive to write general code intended to be called by a
number of users, which may have potentially different ideas about the
most convenient container to store and handle their data. Two possible
alternatives to deal with this problem are:

- Pick a specific type (either an existing or a newly developed one) and
  force everyone (developers and users alike) to use it;
- Develop several versions of the same "general" algorithm and rely on
  the function overloading mechanism to select the correct one when
  compiling;

However, both alternatives are bad design choices. The first one would
tie the library to some specific container, whose functionality would
likely be found to be insufficient soon after being chosen. Any request
by some user to extend the container's functionality would have to be
checked by a developer and all relevant code would have to be changed
and compiled from scratch. The alternative to this waste of development
effort would be to have the specific container choice written on stone
and prevent/discourage any user to try something different, which would
severely limit `gyronimo` flexibility. Moreover, choosing a specific
container would most likely demand costly copy operations to move in
data from other data structures. The second option is not much better
because having a set of overloaded algorithms like

```
double my_general_gyronimo_algorithm(double* A);
double my_general_gyronimo_algorithm(gsl_block* A);
double my_general_gyronimo_algorithm(std`vector<double>& A);
double my_general_gyronimo_algorithm(std`valarray<double>& A);
...
```

would also demand the algorithm's developer to add a newly written
version of the code for each container type requested by users, followed
by compiling the relevant parts of the library.

How does it work?
-----------------

The solution followed in `gyronimo` tries to somehow split the workload
between developers and users, with the heavier share falling towards the
latter (if a user wants to enjoy some exotic functionality, he has to do
a bit more than just demanding its implementation by developers or
complaining about it...). The main idea is to define a very basic
abstract *interface* to a generic continuous container of doubles (e.g.,
first element, size, etc.) that every `gyronimo` developer knows about
and to write every general algorithms in terms of it. It is then the
responsibility of each user to *derive* from such abstract interface an
adaptor to the container of his particular choice and to implement the
very few (only two, actually) virtual functions that are needed. By
passing a reference to their specific adaptor into general algorithms
that take as argument a reference to the abstract interface, the user
will take advantage of the virtual function mechanism in order to get
the library code to seamlessly work with the container type of his
choice. A quite significant detail is that `gyronimo` developers are not
required to get involved into the business. Moreover, because we are
talking about an interface and not a newly defined container structure
by itself, no expensive copies are needed.

Among other things, the abstract class `dblock` defines the minimal
interface that `gyronimo` developers need to know about and demands
every derived class (or adaptors) to implement two simple virtual
functions:

```
class dblock {
 public:
  ...
  virtual const double* begin() const = 0;
  virtual const double* end() const = 0;
  ...
};
```

As an example, suppose some user has a fancy about using a `gsl_block`
(interested readers may check the relevant documentation available
[here](https://www.gnu.org/software/gsl/doc/html/vectors.html)) as a
container to store a continuous collection of doubles that needs to be
passed into some general `gyronimo` algorithm. Defining an adaptor to a
`gsl_block` is as simple as writing

```
class adapter_gsl_block : public dblock {
 public:
  adapter_gsl_block(gsl_block* block) : ptr_(block) {};
  virtual const double* begin() const override {
      return ptr_->data;};
  virtual const double* end() const override {
      return ptr_->data + ptr_->size;};
 private:
  const gsl_block* ptr_;
};
```

After writing the relatively small piece of code above, the user is able
to call general algorithms as

```
// algorithm declaration:
double my_general_gyronimo_algorithm(const dblock& A);

// allocates memory:
gsl_block* b = gsl_block_alloc(100);

// do something with b:
...

// builds the adaptor and calls some algorithm:
adapter_gsl_block adapter(b);
double c = my_general_gyronimo_algorithm(adapter);

// deletes data and frees memory.
gsl_block_free(b);
```

Well, it couldn't be more straightforward...
