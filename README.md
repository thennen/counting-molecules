# counting-molecules

This project is a set of functions to automate the counting and categorization of molecules, specifically tailored for data generated from low temperature scanning probe microscopes.

The functionality of this code was used to generate Fig. 3 of [this publication.](https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.201812334)

## Getting Started

This is a glorified script that makes use of existing python libraries.  Written in python 3+, if unsure, start by installing [Anaconda](https://www.anaconda.com/download).

### Dependencies

* [Python 3+](https://www.anaconda.com/download)
* [SciPy](https://www.scipy.org/)
* [NumPy](http://www.numpy.org/)
* [Matplotlib](https://matplotlib.org/)
* [matplotlib-scalebar](https://pypi.org/project/matplotlib-scalebar/)
* [Scikit-image](https://scikit-image.org/)
* [Scikit-learn](https://scikit-learn.org/stable/)
* [Nanonispy](https://github.com/underchemist/nanonispy)
* [Mahotas](https://mahotas.readthedocs.io/en/latest/)

### Installing

Clone this repository, navigate to said directory, and run:

```
pip install .
```

### Examples

`Helicene_example`, `APT_example` and `APT_2_example` are three example scripts that generate the figures in [the ArXiv article to be posted](arxiv.org/abs/).

Steps to reproduce:

1. Download the data files from this [figshare repo](https://doi.org/10.6084/m9.figshare.19217556) into the examples folder.
2. Run the individual script files in the examples folder to generate each figure.



## Versioning

We use github, see the [tags on this repository](https://github.com/thennen/counting-molecules/tags)

## Authors

* **Jack Hellerstedt** [jhellerstedt](https://github.com/jhellerstedt)
* **Tyler Hennen** [THennen](https://github.com/thennen)

## License

This project is licensed under the MIT License.

## Acknowledgements

This project was inspired by the need to analyze the data that eventually resulted in [this publication](https://onlinelibrary.wiley.com/doi/abs/10.1002/anie.201812334).
We also made use of the data available from [this work.](https://www.nature.com/articles/nchem.2662)
This happened at the [Nanosurf lab.](https://nanosurf.fzu.cz/)
