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
to generate separate code paths optimized and tailored to each listed interaction type.


## A type-aware approach to listed forces


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

Since `ftype` is a runtime value, the whole interaction type dispatch that's controlled through it is managed
by branches and function pointer tables, i.e. the implementation above contains a lot of control flow logic
of the form `if (isPairInteraction(ftype)) {...}`. If on the one hand each interaction type is its own C++ type, we additionally
have overloads available to control flow or template instantiations to create separate code paths that don't
require their flow to be controlled anymore.

On the other hand, we can't just loop over different C++ types. This is illustrated by the fact that we
can't, for instance, write a `for` loop over the `ListedInteractions` tuple from above:

```c++
for (int i = 0; i < interactions.size(); ++i)
{
    // error: i is not a compile time constant
    calc_one_type(std::get<i>(interactions));
}
```

Back in the 90s, the only option would have been to manually unroll the tuple loop.
Given that there's several dozen different types of interactions implemented in Gromacs, this would have been
quite cumbersome and the union data type and function pointer table approach for `InteractionDefinitions` was the only
practicable solution.

Luckily, things have changed. Consider the following code.

```c++
template<class Buffer>
auto reduceListedForces(const ListedInteractions& interactions,
                        const std::vector<gmx::RVec>& x,
                        Buffer* forces)
{
    std::array<real, std::tuple_size<ListedInteractions>::value> energies;

    // calculate one bond type
    auto computeForceType = [forces, &x, &energies](const auto& ielem) {
        real energy = computeForces(ielem.indices, ielem.parameters, x, forces);
        energies[FindIndex<std::decay_t<decltype(ilem)>, ListedInteractions>{}] = energy;
    };

    // calculate all bond types, returns a tuple with the energies for each type
    for_each_tuple(computeForceType, interactions);

    return energies;
}
```

With the help of a generic lambda and C++17's `std::apply` in the one-liner `for_each_tuple`, we can
generate the loop over the different types in the tuple quite effortlessly. While `reduceListedForces`
implements a loop over the interaction types, the next layer, `computeForces` implements a loop over all
interactions of a given type:


```c++
template <class Index, class InteractionType, class Buffer>
real computeForces(const std::vector<Index>& indices,
                const std::vector<InteractionType>& iParams,
                const std::vector<gmx::RVec>& x,
                Buffer* forces)
{
    real Epot = 0.0;

    for (const auto& index : indices)
    {
        Epot += dispatchInteraction(index, iParams, x, forces);
    }

    return Epot;
}
```

We're now down to the level of individual bonds, angles and dihedrals. At this
point, the next steps depend on the actual type of the interaction.
But instead of dispatching each harmonic bond, cubic bond, harmonic angle and so on
to their seperate paths just yet, we just differentiate based on the number of interaction
centers for now. Through overload resolution, the appropriate version `dispatchInteraction`
gets called now, such as this one for the case of 2-center interactions:

```c++
template <class Buffer, class TwoCenterType>
std::enable_if_t<IsTwoCenter<TwoCenterType>::value, real>
dispatchInteraction(const InteractionIndex<TwoCenterType>& index,
                    const std::vector<TwoCenterType>& bondInstances,
                    const std::vector<gmx::RVec>& x,
                    Buffer* forces,
                    PbcHolder pbc)
{
    int i = std::get<0>(index);
    int j = std::get<1>(index);
    const gmx::RVec& x1 = x[i];
    const gmx::RVec& x2 = x[j];
    const TwoCenterType& bond = bondInstances[std::get<2>(index)];

    gmx::RVec dx;
    // calculate x1 - x2 modulo pbc
    pbc.dxAiuc(x1, x2, dx);
    real dr2 = dot(dx, dx);
    real dr  = std::sqrt(dr2);

    real force = 0.0, energy = 0.0;
    std::tie(force, energy) = bondKernel(dr, bond);

    // avoid division by 0
    if (dr2 != 0.0)
    {
        force /= dr;
        detail::spreadTwoCenterForces(force, dx, &(*forces)[i], &(*forces)[j]);
    }

    return energy;
}
```

We can see that coordinate retrieval, computation of the scalar distance and the spreading
of the scalar part of the force to the two centers is actually shared between all the different
types of 2-center interactions. The only remaining thing to do now is to call the actual kernel
to compute the force. Since `bond` has a distinct type, we can again use overload resolution:

```c++
template <class T>
auto bondKernel(T dr, const HarmonicBond& bond)
{
    return harmonicScalarForce(bond.forceConstant(), bond.equilDistance(), dr);
}
```

and call the actual kernel, which in its simplest form for a harmonic bond looks like this:

```c++
template <class T>
std::tuple<T, T> harmonicScalarForce(T k, T x0, T x)
{
    real dx  = x - x0;
    real dx2 = dx * dx;

    real force = -k * dx;
    real epot = 0.5 * k * dx2;

    return std::make_tuple(force, epot);

    /* That was 6 flops */
}
```

That's it! From here on, we just had to add a separate dispatch for the 3- to 5-center interactions
and add the type-aware wrappers for all the different kernels implemented in Gromacs.

What have we achieved?
We've avoided future maintenance headaches by not writing a fragile translation layer, routed
the listed force data flow through type-aware channels with maximum code reuse and control logic
already baked into the type dispatch, and all of that while leveraging the code from Gromacs that actually matters:
the force kernels that implement all the physics!

## The fun's only starting now

We're not done yet, howere. The best part is yet to come, stay tuned.

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
