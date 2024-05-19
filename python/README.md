# rpa_finder


`rpa_finder` is a python package for the systematic discovery of Robust Perfect Adaptation (RPA) in chemical reaction networks.

The package allows for finding all the RPA properties implemented in a
deterministic chemical reaction system through the enumeration of *labeled buffering structures*, that are in one-to-one correspondence with elementary RPA properties of the system. 
In other words, it facilitates the determination of the dependencies of steady-state concentrations and reaction rates on all the reaction parameters (and the values of conserved quantitites if any) in a model-independent manner.

Once we have the list of labeled buffering structures, one can, for example, do the following:

- Finds all the concentrations and reactions affected by the change of a chosen reaction parameter. 
- Finds all the reactions that affect the steady-state value of the concentration of a chosen species. 
- Finds integral control realizing each RPA property represented by a (labeled) buffering structure.

## Installation

The package can be installed via pip as

```bash
pip install rpa_finder
```


## Usage

We here briefly describe basic usage. For more examples, see jupyter notebooks in `examples` directory.

### Finding labeled buffering structures

First, prepare a list of reactions and create a reaction system object. For example, let us consider a reaction network with three species $\{v_1,v_2,v_3\}$ and the following five reactions:

- $e_1: \emptyset \to v_1$

- $e_2: v_1 \to v_2$

- $e_3:v_2 \to v_3$

- $e_4:v_3 \to v_1$

- $e_5:v_2 \to \emptyset$


To prepare a reaction system for this network, let us first define a reaction network in a text format as 

```python 
network1 = """
"", "v1"
"v1", "v2"
"v2", "v3"
"v3", "v1"
"v2", ""
"""
```

Note that $\emptyset$ is represented as an empty string `""` in the text format.

We can create a `ReactionSystem` object corresponding to the network. 
Call the method `enumerate_labeled_buffering_structures` to enumerate labeled buffering structures for this network: 

```python 
r_system1 = ReactionSystem(network1)
```


A labeled buffering structure is expressed by a quadruple, 
$( P, V_P, E_P, \mathcal E_P)$, where 

- $P$ denotes a set of parameters that are perturbed
- $V_P$ denotes the set of species affected by the perturbation of the parameters in $P$ 
- $E_P$ denotes the set of reactions affected by the perturbation of the parameters in $P$
- $\mathcal E_P$ denote a set of added reactions to make the subnetwork $(V_P,E_P \cup \mathcal E_P)$ output-complete 

Labeled buffering structures for a given reaction system can be found by the following function, 

```python
lbs_list = r_system1.enumerate_labeled_buffering_structures()
```

The output is 

```python
[
    [[0], [0, 1, 2], [0, 1, 2, 3, 4], []],
    [[1], [0], [], [1]],
    [[2], [0, 2], [1, 2, 3], []],
    [[3], [2], [], [3]],
    [[4], [0, 1, 2], [1, 2, 3], [4]]
]
```

The species and reactions are indicated by indices. Note that, unlike the *Mathematica* version, indicees start with zero. 
To use names for species, one can use `lbs_to_name`:

```
lbs_name = [ r_system1.lbs_to_name(l) for l in lbs_list]
```

The content of `lbs_name` is

```
[
    [[0], ['v1', 'v2', 'v3'], [0, 1, 2, 3, 4], []],
    [[1], ['v1'], [], [1]],
    [[2], ['v1', 'v3'], [1, 2, 3], []],
    [[3], ['v3'], [], [3]],
    [[4], ['v1', 'v2', 'v3'], [1, 2, 3], [4]]
]
```

We can read off the dependencies of concentrations and reaction rates on all the system parameters. For example, from the third labeled buffering structure, we can see that the perturbation of the parameter of reaction $e_3$ (recall that the index  of reactions in the code starts with zero, while the index in the original reaction list starts with one) affects the concentrations of $v_1$ and $v_3$ and the rates of reactions $e_2,e_3,e_4$. It does not affect the concentration of $v_2$ and reaction rates of $e_1$ and $e_5$, which means that they exhibit RPA with respect to this parameter. 

### Finding integral control

For every labeled buffering structure, one can find integral control realizing the RPA property. For example, let us find integrator equations for `lbs_list[1]`:

```python

integrators = r_system1.find_integrators_from_lbs(
    lbs_list[1], symbol_mode='name'
    )

for i in range(4):
    display( 
        Markdown(
            "$\\frac{d}{dt}" + integrators[2*i] 
            + "=" 
            + integrators[2*i+1] + "$"
            ) 
        )

```

The output is 

$$\frac{d}{dt}\left[\begin{matrix}\end{matrix}\right]=\left[\begin{matrix}\end{matrix}\right]$$

$$\frac{d}{dt}\left[\begin{matrix}v_{1} + v_{2} \\ v_{3}\end{matrix}\right]=\left[\begin{matrix}r_{1} - r_{3} + r_{4} - r_{5} \\ r_{3} - r_{4}\end{matrix}\right]$$

$$\frac{d}{dt}\left[\begin{matrix}- v_{1}\end{matrix}\right]=\left[\begin{matrix}- r_{1} + r_{2} - r_{4}\end{matrix}\right]$$

$$\frac{d}{dt}\left[\begin{matrix}v_{2}\end{matrix}\right]=\left[\begin{matrix}r_{2} - r_{3} - r_{5}\end{matrix}\right]$$

In general, there are four sets of equations constituting integrators. 


Note that the indices specifying the parameters is arranged in the order of $(\vec k, \vec \ell)$, where $\vec k$ are rate parameters and $\vec \ell$ are the values of conserved quantities. 
As a basis of conserved quantities, those return by `scipy.linalg.null_space(s)` is used, where `s` is the stoichiometric matrix of the reaction system. 



## Testing
Tests can be run by `python -m unittest discover -s tests`

## Reference

- Y. Hirono, A. Gupta, M. Khammash, "_Complete characterization of robust perfect adaptation in biochemical reaction networks_," [arXiv:2307.07444](https://arxiv.org/abs/2307.07444). 

## Contact

If you have any questions or suggestions, feel free to drop an email to [Yuji Hirono](mailto:yuji.hirono@gmail.com).


## License

RPAFinder is licensed under the [MIT License](https://opensource.org/license/mit/). See [LICENSE](LICENSE) for details.
