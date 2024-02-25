# Markov Categories in C++

This is a companion implementation for the paper *Hidden Markov Models and the Bayes Filter
in Categorical Probability* by Fritz et al., demonstrating how Markov categories can be represented as a C++ abstract class interface.

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


