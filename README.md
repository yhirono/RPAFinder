# RPAFinder

**RPAFinder** is a _Mathematica_ package that offers useful functions for studying the impact of parameter variations and explore the intricate relationships within a chemical reaction network effectively. 
It facilitates the determination of the dependencies of steady-state concentrations and reaction rates on all the reaction parameters (and the values of conserved quantitites if any) in a model-independent manner.
This allows for systematic discovery of Robust Perfect Adaptation (RPA) in a reaction network.
The determination is achieved through the enumeration of all the labeled buffering structures for a given reaction network. 
Once we have the list of labeled buffering structures, one can, for example, do the following:

- Find all the concentrations and reactions affected by the change of a chosen reaction parameter. 
- Find all the reactions that affect the steady-state value of the concentration of a chosen species. 

See "demo.nb" for demonstration. More examples can be found in "examples.nb".


# Usage

## Installation

Download the package file "RPAFinder.wl". 
One can load the package from a _Mathematica_ notebook by

```mathematica
<< "RPAFinder.wl"
```

## Finding labeled buffering structures

First, prepare a list of reactions and create a reaction system object. For example, let us consider a reaction network with three species $\{v_1,v_2,v_3\}$ and the following five reactions:

- $e_1: \emptyset \to v_1$

- $e_2: v_1 \to v_2$

- $e_3:v_2 \to v_3$

- $e_4:v_3 \to v_1$

- $e_5:v_2 \to \emptyset$

To prepare a reaction system for this network, call the method `makeReactionSystem[list]`

```mathematica 
system = makeReactionSystem[ 
	{
		{0,  v1},
		{v1, v2},
		{v2, v3},
		{v3, v1},
		{v2, 0}
	}
];
```
The argument of this function is a list of reactions. For $\emptyset$, use 0. 

A labeled buffering structure is expressed by a quadruple, 
$( P, V_P, E_P, \mathcal E_P)$, where 

- $P$ denotes a set of parameters that are perturbed
- $V_P$ denotes the set of species affected by the perturbation of the parameters in $P$ 
- $E_P$ denotes the set of reactions affected by the perturbation of the parameters in $P$
- $\mathcal E_P$ denote a set of added reactions to make the subnetwork $(V_P,E_P \cup \mathcal E_P)$ output-complete 

Labeled buffering structures for a given reaction system can be found by the following function, 

```mathematica 
enumerateLabeledBufferingStructures[system]
```

The output is 

```mathematica 
{
	{{1}, {v1, v2, v3}, {1, 2, 3, 4, 5}, {}}, 
	{{2}, {v1}, {}, {2}}, 
	{{3}, {v1, v3}, {2, 3, 4}, {}}, 
	{{4}, {v3}, {}, {4}}, 
	{{5}, {v1, v2, v3}, {2, 3, 4}, {5}}
}
```

We can read off the dependencies of concentrations and reaction rates on all the system parameters. For example, from the third labeled buffering structure, we can see that the perturbation of the parameter of reaction $e_3$ affects the concentrations of $v_1$ and $v_3$ and the rates of reactions $e_2,e_3,e_4$. It does not affect the concentration of $v_2$ and reaction rates of $e_1$ and $e_5$, which means that they exhibit RPA with respect to this parameter. 

Note that the indices specifying the parameters is arranged in the order of $(\vec k, \vec \ell)$, where $\vec k$ are rate parameters and $\vec \ell$ are the values of conserved quantities. 
As a basis of conserved quantities, those return by `NullSpace@Transpose@s` is used where `s` is the stoichiometric matrix of the reaction system. 

Instead, if we use `enumerateLabeledBufferingStructuresByIndex[system]`, the species are denoted by indices. 


# Reference

- Y. Hirono, A. Gupta, M. Khammash, "_Complete characterization of robust perfect adaptation in biochemical reaction networks_," [arXiv:????]. 


# Contact

If you have any questions or suggestions, feel free to drop an email to [Yuji Hirono](mailto:yuji.hirono@gmail.com).


# License

RPAFinder is licensed under the [MIT License](https://opensource.org/license/mit/). See [LICENSE](LICENSE) for details.
