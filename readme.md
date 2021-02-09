::gyronimo:: - gyromotion for the people, by the people -
=========================================================
*An object-oriented library for gyromotion applications in plasma physics.*  

Philosophy and purpose:
-----------------------
Have you ever had a bright and promising idea about gyromotion in plasmas that
just faded away the moment you realised the amount of non-trivial, tedious,
unrewarding, non-physics details you would have to implement before you could
get a simple glimpse over the results?

`gyronimo` was designed to address this problem. It provides a library of
objects that take care of most tedious tasks in gyromotion simulations, allowing
developers to quickly implement and test their ideas in order to understand if
further work and optimisation really worth the trouble. Built over
object-oriented programming concepts like inheritance, polymorphism, and
specialisation, `gyronimo` provides abstract algorithms to many tasks that
developers may take advantage of, either by using ready-to-use objects or by
derive new ones adapted to their own needs, in order to quickly finish a
"working" workflow. In a later stage, derived objects and procedures may be
specialised and further optimised for performance, without the need to
significantly change the logic structure of the initial code.

In brief, the goal of `gyronimo` is to keep maintenance efforts to a minimum,
code readability to its best, and performance to the level that is really needed
by each specific gyromotion application.

"What's in a name?", wondered Juliet.
-------------------------------------
Quite obviously, `gyronimo` is an acronym for GYRO-MOtion. To those with a sharp
eye and an inquisitive mind, the author must confess that he has absolutely No
Idea on how the `ni` managed to sneak into the middle.

Once, in the wild west of yore, Apache warriors were said to charge their foes
shouting "Geronimo!", the nickname of their fierce leader; Other sources say the
same cry was uttered instead by the fleeing Mexican soldiers, invoking their
saint protector "Santo Geronimo"; Whichever the actual cause might have been,
library users are expected to occasionally shout an enthusiastic "gyronimo!"
while happily coding their favorite gyromotion application.

Installation and dependencies:
------------------------------
The library requires a c++ compiler supporting the **c++20** standard (e.g.,
gcc-10.1.0). HTML documentation can be extracted from the source code by
doxygen.

- Edit the makefile in the folder `build` and adjust the variables (e.g.,
  compiler, compiling options, and installation folders) to your system's
  requirements;
- Run `make lib` to generate the static library `libgyronimo.a`; copy the static
  library and the `include` folder to wherever your system needs them;
- Run `make dox` to extract HTML documentation into the folder `gyronimo_dox`;
  move it to wherever you find it more convenient;
- Run `make clean` to cleanup all the mess.

Applications developed using `gyronimo` may, eventually, have to be compiled
and/or linked against the following libraries:

- GNU Scientific Library [GSL](https://www.gnu.org/software/gsl).
- odeint, an ODE c++ library distributed with [boost](https://www.boost.org).
