# Markov Categories in C++

This is a companion implementation for the paper *Hidden Markov Models and the Bayes Filter
in Categorical Probability* by Fritz et al., demonstrating how Markov categories can be represented as a C++ abstract class interface.

### Structure of the code

* `mcat.h` defines the abstract interface for Markov categories. Any algorithm given as a string diagram should operate only on pointers to the classes `Obj` and `Mor` defined herein.
* `finstoch.h` defines an implementation of the above interface for the category *FinStoch*. If you wish to run a given abstract algorithm on finite probability spaces, you should execute it on instances of these classes (that you created explicitly, to make them contain interesting probabilistic maps). The implementation uses the `armadillo` library, which implements a subset of the MATLAB API in C++, for linear algebra.
* `bayes.cpp` contains an implementation of the instantiated Bayes filter as an abstract algorithm on Markov categories (see the function `ibf`). Moreover, it contains (in `main`) some simple code that tests this implementation and outputs diagnostic input on a simple HMM with two states, of which one is absorbing, and a noisy observation channel.
* `test.cpp` contains miscellaneous stress and functionality tests of the API without particular implications or structure.

### Adding more categories

To add support for an additional category, it is sufficient to derive classes for objects (from `Obj`) and morphisms (from `Mor`) implementing the pure virtual functions in `mcat.h`, in particular

* For any object,
  * `Obj::id`, returning the identity morphism
  * `Obj::copy_n`, returning an n-branch "copy" morphism
  * `Obj::term`, returning the unique morphism to the terminal object
  * `Obj::tensor_length`, returning the number of basic objects an object is tensored from
  * `Obj::is_equal`, performing an equality test on two objects
  * `Obj::tensor`, tensoring two objects
  * (`Obj::debug_print`, printing any debug information you would like about the object)
  * You may also provide an optimised implementation of `Obj::select`, returning the projection morphism to the ith factor in a tensor.
* For any morphism,
  * `Mor::dom`, returning its domain object
  * `Mor::cod`, returning its codomain object
  * `Mor::cond`, conditioning on the last `n` factors of the codomain
  * `Mor::compose`, composing with another morphism
  * `Mor::tensor`, tensoring two morphisms
  * (`Mor::debug_print`, printing any debug information you would like about the morphism)

### Caveats

The consumer of the Markov category API is responsible for freeing any memory (with `operator delete`), but the example code is not diligent about doing this and therefore leaks memory.
