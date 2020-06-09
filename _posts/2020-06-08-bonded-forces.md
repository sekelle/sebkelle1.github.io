---
layout: post
title: Bonded forces in MD simulations
subtitle: Rewiring the dataflow in Gromac's implementation of bonded forces
#gh-repo: daattali/beautiful-jekyll
#gh-badge: [star, fork, follow]
tags: [PRACE, NBLIB]
comments: true
---

While the main purpose of NBLIB was always to expose non-bonded force calculation capabilites,
we still needed a way to obtain bonded forces for the systems set up with our newly developed
topology functionality.

## Making bonded forces available to NBLIB

At the heart of it, calculating forces for different types of bonds is quite simple:
for an array of input 3D coordinates, for each bond you need a pair of indices to refer to
the two participating particles plus a handful of parameters like spring constants and equilibrium distances
to calculate the bond forces. And it this works just the same for forces due to angles and dihedrals.
The commonality between these types of forces is that they act on predefined sets of coordinates. That's why
in GROMACS, they are called "listed" interactions - there's a list which contains the indices along with the
parameters for each bond, angle or dihedral. This list is normally constant during a simulation run and is thus
part of the system topology.

NBLIB implements system topology and listed interactions are naturally a part of it. The data format looks
like this:

```c++
template <class Interaction>
struct InteractionData
{
    std::vector<Index<Interaction>> indices;
    std::vector<Interaction>        parameters;
};
```

For each type of interaction, we store the interaction indices plus the interaction parameters.
An example for `Interaction` would be `HarmonicBond`, which looks like this:

```c++
struct HarmonicBond : public std::tuple<real, real>
{
    real& forceConstant() { return std::get<0>(*this); }
    real& equilDistance() { return std::get<1>(*this); }
};
```
The `Index` traits class deduces to `std::array<int, 3>`, because for each harmonic bond, we need two `int`s
for the coordinate indices and a third `int` to look up the bond parameters in the `parameters` vector.
For bonds and dihedrals, the `Index` trait adds an additional one or two `int`s to hold the additional coordinate indices.
Finally, each type of interaction is stored in a `std::tuple`, such that the complete definition of all listed forces
in the NBLIB topology looks like this:

```c++
using ListedInteractions = std::tuple<InteractionData<HarmonicBond>, ..., InteractionData<HarmonicAngle>, ...>;
```

## Listed interactions in Gromacs

How is the same data stored in Gromacs? The format looks like this:
```c++
struct InteractionDefinitions
{
    std::vector<t_iparams> iparams;
    std::array<std::vector<int>, F_NRE> il;
};
```
This includes all interaction types, i.e. `t_iparams` is a union type which can hold the parameters of any type.
Thus `iparams` corresponds to the concatenation of the NBLIB `parameters` vectors for all interaction types.
The other member called `il` contains the indices for each interaction type, where `F_NRE` is the number of interaction types
that Gromacs supports. Just like NBLIB, Gromacs stores a vector of indices for each type. The only difference
is that the `value_type` is always `int` for Gromacs, while for NBLIB, it's `std::array<int, N+1>` for an `N-`center
interaction type. In summary this means that in Gromacs too, a harmonic bond is represented as three integers, two for the
coordinates and a third to look up the bond parameters in a look-up table. But the fundamental difference is:
NBLIB has a different _(C++)-type_ for each interaction type while Gromacs uses a single (C++)-type for all kinds of interactions.


## To convert or not to convert?

It seemed like a straight forward thing to implement. We were going to write a translation function
to convert from `ListedInteractions` to Gromac's `InteractionDefinitions`, call `calc_listed(const InteractionDefinitions& idef)`
and voilà, there's our listed forces.

Pretty soon, however, we began to have some doubts. The fact that such a conversion would rely on many
implementation details, such as the specific positions of each interaction type in the `il` array, for instance, was one concern.
By design, such a layout translation would be inherently fragile with respect to future changes in the underlying code
and create a future maintenance burden. Another issue was that `calc_listed`, the entry point to listed forces in Gromacs,
had accumulated a lot of "baggage" over the years. For example, the argument signature of `calc_listed` contains a translation
table to convert from (node)-local to global atom index. This table gets passed all the way down to the kernel level where forces
for individual bonds are calculated. Its only purpose is to translate the local index to the global one in the error message
for an obscure case that nobody remembers about. Another example is the use of multi-purpose data structures, such as `t_forcerec`.
Yes, `calc_listed` computes the listed forces, but where do these forces actually end up? Turns out there's one of those `t_forcerecs`
in the argument list which holds a sub-object called `bonded_threading` that (among other things) has the force buffer we're looking for.
The `t_forcerec` object also holds data and parameters related to non-bonded force computations and various hardware related information
such as resources allocated for calculations on GPUs. In short, `t_forcerec` is used everywhere where forces are computed, which,
in an MD code, is basically _everywhere_.

This means that if we wanted to obtain listed forces in NBLIB through `calc_listed`, the only realistic option would
be to create a `t_forcerec` object and fill it with dummy values except for the part actually relevant for bonded forces.
At this point, we where asking ourselves:
>_What if we rewired the dataflow to feed the NBLIB interaction data directly to the force kernels in Gromacs?_

In fact, we quickly realized that we'd be doing ourselves a big disservice if we just threw away the interaction types
by which we'd already neatly grouped all the listed interactions. C++ compilers have become powerful code generators,
but they operate on types. If we don't have different types, we forego the possibility of leveraging the compiler
to generate separate code paths optimizied and tailored to each listed interaction type.


## A type-aware approach to bonded forces


In essence, the implementation of `calc_listed(const InteractionDefinitions& idef, ...)` in Gromacs looks like this:

```c++
void calc_listed(const InteractionDefinitions& idef, ...)
{
    // manage timing and multi-threading 

    for (int ftype = 0; ftype < F_NRE; ++type)
    {
        // branch out and descend stack for 2 intermediate functions based on
        // the type of interaction that ftype corresponds to
        // then call a function from a pointer table

        bondFunction* bonded = bondedInteractionFunctions[ftype]; 

        // compute all forces for ftype
        bonded(idef.iparams, idef.il[ftype], ...);
    }

    // reduce thread output
}
```



<!--
**Here is some bold text**

## Here is a secondary heading

Here's a code chunk:

~~~
var foo = function(x) {
  return(x + 5);
}
foo(3)
~~~

And here is the same code yet again but with line numbers:

{% highlight javascript linenos %}
var foo = function(x) {
  return(x + 5);
}
foo(3)
{% endhighlight %}

## Boxes
You can add notification, warning and error boxes like this:

### Notification

{: .box-note}
**Note:** This is a notification box.

### Warning

{: .box-warning}
**Warning:** This is a warning box.

### Error

{: .box-error}
**Error:** This is an error box.
-->
