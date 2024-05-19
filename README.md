# RPAFinder



**RPAFinder** is a package for the systematic discovery of Robust Perfect Adaptation (RPA) in chemical reaction networks. This utility is implemented in both Python and Mathematica, offering flexibility for users in different computational environments.

The package allows for finding all the RPA properties implemented in a
deterministic chemical reaction system through the enumeration of *labeled buffering structures*, that are in one-to-one correspondence with elementary RPA properties of the system. 
In other words, it facilitates the determination of the dependencies of steady-state concentrations and reaction rates on all the reaction parameters (and the values of conserved quantitites if any) in a model-independent manner.

Once we have the list of labeled buffering structures, one can, for example, do the following:

- Finds all the concentrations and reactions affected by the change of a chosen reaction parameter. 
- Finds all the reactions that affect the steady-state value of the concentration of a chosen species. 
- Finds integral control realizing each RPA property represented by a (labeled) buffering structure.

For more details, please refer to the `README` file of each version:
- [Python version](python/README.md)
- [Mathematica version](mathematica/README.md)



## Reference

- Y. Hirono, A. Gupta, M. Khammash, "_Complete characterization of robust perfect adaptation in biochemical reaction networks_," [arXiv:2307.07444](https://arxiv.org/abs/2307.07444). 

## Contact

If you have any questions or suggestions, feel free to drop an email to [Yuji Hirono](mailto:yuji.hirono@gmail.com).


## License

RPAFinder is licensed under the [MIT License](https://opensource.org/license/mit/). See [LICENSE](LICENSE) for details.
