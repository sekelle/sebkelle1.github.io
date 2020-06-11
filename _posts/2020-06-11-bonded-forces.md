---
layout: post
title: Bonded forces in GROMACS and NB-LIB
subtitle: Rewiring the dataflow in Gromac's implementation of bonded forces
#gh-repo: daattali/beautiful-jekyll
#gh-badge: [star, fork, follow]
tags: [PRACE, NB-LIB]
comments: true
---

The main purpose of [NB-LIB](https://gitlab.com/gromacs/nb-lib) was always to expose non-bonded force calculation capabilites.
The "NB",  after all, stands for "non-bonded".
That didn't mean however that we, the NB-LIB team got to ignore the bonded part.
We still needed a way to obtain bonded forces for the systems set up with our newly developed topology functionality.
One way or another.

## Making bonded forces available to NB-LIB

At the heart of it, calculating forces for different types of bonds is quite simple:
for an array of input 3D coordinates, for each bond you need a pair of indices to refer to
the two participating particles plus a handful of parameters like spring constants and equilibrium distances
to calculate the bond forces. And this works just the same for forces due to angles and dihedrals.
The commonality between these types of forces is that they act on predefined sets of coordinates. That's why
in [GROMACS](http://www.gromacs.org/), they are called "listed" interactions - there's a list which contains the indices along with the
parameters for each bond, angle or dihedral. This list is normally constant during a simulation run and is thus
considered part of the system topology.

NB-LIB implements a system topology and listed interactions are naturally a part of it.
The data format looks like this:

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
in the NB-LIB topology looks like this:

```c++
using ListedInteractions = std::tuple<InteractionData<HarmonicBond>, ..., InteractionData<HarmonicAngle>, ...>;
```
&nbsp;

## Listed interactions in GROMACS

How is the same data stored in GROMACS? The format looks like this:
```c++
struct InteractionDefinitions
{
    std::vector<t_iparams> iparams;
    std::array<std::vector<int>, F_NRE> il;
};
```
This includes all interaction types, i.e. `t_iparams` is a union type which can hold the parameters of any type.
Thus `iparams` corresponds to the concatenation of the NB-LIB `parameters` vectors for all interaction types.
The other member called `il` contains the indices for each interaction type, where `F_NRE` is the number of interaction types
that GROMACS supports. Just like NB-LIB, GROMACS stores a vector of indices for each type. The only difference
is that the `value_type` is always `int` for GROMACS, while for NB-LIB, it's `std::array<int, N+1>` for an `N-`center
interaction type. In summary this means that in GROMACS too, a harmonic bond is represented as three integers, two for the
coordinates and a third to look up the bond parameters in a look-up table. But the fundamental difference is:
NB-LIB has a different _(C++)-type_ for each interaction type while GROMACS uses a single (C++)-type for all kinds of interactions.


## To convert or not to convert?

It seemed like a straight forward thing to implement. We were going to write a translation function
to convert from `ListedInteractions` to Gromac's `InteractionDefinitions`, call `calc_listed(const InteractionDefinitions& idef)`
and voilà, there's our listed forces.

Pretty soon, however, we began to have some doubts. The fact that such a conversion would rely on many
implementation details, such as the specific positions of each interaction type in the `il` array, for instance, was one concern.
By design, such a layout translation would be inherently fragile with respect to future changes in the underlying code
and create a future maintenance burden. Another issue was that `calc_listed`, the entry point to listed forces in GROMACS,
had accumulated a lot of "baggage" over the years. For example, the argument signature of `calc_listed` contains a translation
table to convert from (node)-local to global atom index. This table gets passed all the way down to the kernel level where forces
for individual bonds are calculated. Its only purpose is to translate the local index to the global one in the
[error message](https://gitlab.com/gromacs/gromacs/-/blob/master/src/gromacs/listed_forces/pairs.cpp#L77)
for the case that the distance of a listed 1-4 interaction using tabulation exceeds the table limit.
But only if either energies, virials or free energy are computed. Yet, the lookup table has to be part of the function signature
for all the other interaction types as well due to the function pointer table (see below) requirint a uniform signature.
Another example is the use of multi-purpose data structures, such as `t_forcerec`.
Yes, `calc_listed` computes the listed forces, but where do these forces actually end up? Turns out there's one of those `t_forcerecs`
in the argument list which holds a sub-object called `bonded_threading` that (among other things) has the force buffer we're looking for.
The `t_forcerec` object also holds data and parameters related to non-bonded force computations and various hardware related information
such as resources allocated for calculations on GPUs. In short, `t_forcerec` is used everywhere where forces are computed, which,
in an MD code, is basically _everywhere_.

This means that if we wanted to obtain listed forces in NB-LIB through `calc_listed`, the only realistic option would
be to create a `t_forcerec` object and fill it with dummy values except for the part actually relevant for bonded forces.
At this point, we where asking ourselves:
>_What if we rewired the dataflow to feed the NB-LIB interaction data directly to the force kernels in GROMACS?_

In fact, we quickly realized that we'd be doing ourselves a big disservice if we just threw away the interaction types
by which we'd already neatly grouped all the listed interactions. C++ compilers have become powerful code generators,
but they operate on types. If we don't have different types, we forego the possibility of leveraging the compiler
to generate separate code paths optimized and tailored to each listed interaction type.


## A type-aware approach to listed forces


In essence, the implementation of `calc_listed(const InteractionDefinitions& idef, ...)` in GROMACS looks like this:

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
of the form `if (someProperty(ftype)) {...}`. If on the one hand each interaction type is its own C++ type, we additionally
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

About 20 years ago when GROMACS was first written, the only practicable option would have been to manually unroll the tuple loop.
Given that there's several dozen different types of interactions implemented in GROMACS, this would have been
quite cumbersome and the union data type and function pointer table approach for `InteractionDefinitions` was the only
feasible solution.

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

That's it! Not a single `if(ftype)` required.
From here on, we just had to add a separate dispatch for the 3- to 5-center interactions
and add the type-aware wrappers for all the different kernels implemented in GROMACS.

What have we achieved? We've
* avoided future maintenance headaches by not writing a fragile translation layer
* routed the listed force data flow through type-aware channels with
* maximum code reuse and
* control logic already baked into the type dispatch

And all of the above while leveraging the code from GROMACS that actually matters:
the force kernels that implement all the physics!


## Extra toppings

We could grab a cup of coffee now and focus again on the non-bonded part, which is NB-LIB's main objective, and we should!
But, what if, like at the time of writing this, it's a rainy holiday and we knew nothing better to do than to continue experimenting a bit?
Should we vectorize the code? And then there's also an idea that Berk Hess had, but never got around to try out.
He said that, instead of processing all types of interactions separately,
when an interaction involves coordinates `i` and `j`, one should also evaluate all the other
types of interactions that involve these two coordinates. Seems like we have just the right tools now to give that a shot.

### A quick roofline analysis

To address the first question, consider a simple harmonic bond for instance.
In the [GROMACS source](https://gitlab.com/gromacs/gromacs/-/blob/master/src/gromacs/listed_forces/bonded.cpp#L450),
we can see that each
of them costs 59 flops, which drops to 46 if we use the version above that doesn't do free energy perturbation.
How many bytes of data do we have to access for that? Assuming single precision,
* 6 `float`s = 24 bytes to load the two coordinates
* 6 `float`s = 24 bytes to store the two forces
* 3 `int`s = 12 bytes for the `i`, `j`, and parameter index
* 2 `float`s = 8 bytes for the two bond parameters

Total: 68 bytes. A [Skylake core](https://en.wikichip.org/wiki/intel/microarchitectures/skylake_(server)) for example
has two AVX-512 units. Thus, it can perform up to 4 flops/cycle for non-vectorized code and up to 64 flops/cycle
if both units perform 512-bit FMA instructions in single precision, or 10 Gflop/s in scalar mode and 160 Gflop/s
vectorized, assuming a clock of 2.5 GHz. How much bandwidth would we need to sustain these processing speeds?
According to the [roofline model](https://dl.acm.org/doi/10.1145/1498765.1498785), performance P = &beta; B,
where &beta; is the operational intensity in flops/byte and B denotes the bandwidth. For our example,
&beta; = 0.68, which means that we need a bandwidth of about 1.5 GB/s for every Gflop/s, so 15 GB/s for the scalar case
and 240 GB/s for AVX-512. 

The bottom line of these metrics is that if the number of atoms per CPU core is so small that the coordinate and force
buffers fit into the L2 cache, it makes sense to vectorize the code. If the bond force computations
are served from main memory, memory bandwith is the bottleneck, and if the data fits into L2, bandwidth is again the
bottleneck if the code is vectorized (L2 data bandwidth on Skylake is 64 B/cycle).

### Synthetic aggregate interaction types

Let's consider and angle between particles at positions `i`, `j` and `k`. Then there's almost certainly going
to be two bonds `i-j` and `j-k` as well. After computing the angle forces when the three coordinates are already
fetched from memory, we could also calculate those two bond forces, reuse the input coordinates and store back the sum of
all the forces. If we managed to do that, the two bonds would essentially be getting a free ride in a bandwidth limited scenario.

One way of implementing this would be through an additional, synthetic aggregate interaction type:
```c++
template<class Angle, class Bond>
struct ThreeCenterAggregate
{
    Angle angle;
    Bond  bonds[2];
};
```

We can just append this type to the `ListedInteractions` tuple and send it down the pipeline like all the other interaction types.
At the `dispatchInteraction` level we can then plug together the angle and bond components.
The routines for that we already have lying around and we can select the right kernels based on the `Angle` and `Bond` types.

At setup time, we will have to go through all bonds and see whether we can find an appropriate angle, at a cost of `NlogM`,
where `N` is the number of bonds and `M` the number of angles.
This means that interactions are subtracted from the `Angle` and `Bond` lists and added to the one for `ThreeCenterAggregate<Angle, Bond>`.
The best place to do this is in a `ListedCalculator` object. Its constructor receives the `ListedInteractions` from the topology,
and then creates the aggregate transformations. By doing that, we can make aggregates an implementation detail of the `ListedCalculator`.
Both the topology above and the whole `reduceListedForces` pipeline below `ListedCalculator` are oblivious to whether it
chose to populate any aggregates or not.


## Some remainders

There are a couple additional aspects to implementing listed forces, such as multi-threading, which we've already implemented,
and calculating virials and free energy perturbations (FEP).


### Free energies

The calculation of free energies, or more precisely, the derivative of the potential energy with respect to &lambda;, which interpolates
between two system states A and B, deserves a mention here.
In GROMACS, each interaction type that supports FEP actually has two sets of parameters, one for the A state and an other one for the B state.
The additional memory for the B state is always there in GROMACS, whether the user wants to do FEP or not. 

Instead of adding a second B state to each interaction, one could also envisage just having two separate interactions for the A and B state.
The advantage of that approach would be that the B state could be absent if the user doesn't want to do FEP. Also we'd be saving a lot of
code, because we could just build two separate topologies A and B. They would then supply each a `ListedInteractions` tuple
and we could add an overload to `reduceListedForces` that accepts a second B tuple and the &lambda; parameters.
The only additional ingredients needed are force kernels that take the second parameter set and &lambda; as inputs, so the signature
looks like this:

```c++
template <class T>
std::tuple<T, T, T> harmonicScalarForce(T kA, T kB, T xA, T xB, T x, T lambda);
// return is <F, V, dV/dl>
```

Of course, these kernels already exist in GROMACS.

With that I'd like to conclude and thank the whole NB-LIB team for the productive last two weeks.
