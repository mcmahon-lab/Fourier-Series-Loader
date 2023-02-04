# Fourier Series Loader 

This repository is a user-friendly implementation of the Fourier Series Loader (FSL) method introduced in Moosa, M. & Watts, T.W. *et. al.* (2023)<sup>[1](#how-to-cite-this-code)</sup>. The FSL is a method for preparing quantum states that encode multi-dimensional Fourier series using linear-depth quantum circuits. This method allows for the highly accurate loading of arbitrary complex-valued functions on a quantum computer.

We present four examples to illustrate the use of the FSL method in this repository. In each of these examples, we perform a quantum simulation on Qiskit of loading functions into a quantum state. For ease of use, we have presented these examples as step-by-step tutorials on how to use the FSL method. Hence, interested users can easily reuse the codes provided in this repository for their own use. In this repository, we provide examples of loading periodic functions, non-periodic functions, piece-wise continuous functions, and function of two variables.

Note that this repository does not discuss how the FSL method performs on a near-term noisy quantum computer. For this, we refer the readers to our paper<sup>[1](#how-to-cite-this-code)</sup>, where we have presented the results of the experiments performed on Quantinuum H1 generation quantum computers.

# Organization of the repository

This repository consists of four notebooks, each containing a different example of a quantum simulation of loading a function on a quantum state using the FSL method: 

- **Example_1.ipynb** focuses on loading periodic functions into a quantum state.
- **Example_2.ipynb** focuses on loading non-periodic functions into a quantum state.
- **Example_3.ipynb** focuses on loading piece-wise continuous functions into a quantum state.
- **Example_4.ipynb** focuses on loading functions of two variables into a quantum state.

This repository also contains two code files:

- **uniformly_controlled_rotations.py** contains codes to implement the uniformly controlled rotations, which were introduced in [M&ouml;tt&ouml;nen *et al.* (2005)](https://arxiv.org/abs/quant-ph/0407010) as a method to prepare any n-qubit state. We use the uniformly controlled rotations as a subroutine of the FSL method. (See our paper<sup>[1](#how-to-cite-this-code)</sup> or the example notebooks for further details.) 

- **supplementary.py** includes codes for various functions that are needed for the simulations performed in the examples presented in this repository.

# Getting started 

Each notebook in this repository contains a self-contained example of the implementation of the FSL method. Therefore, users can use them in any order based on their needs. Nevertheless, we recommend users start with notebook *Example_1.ipynb* as it contains a detailed discussion about each step of the FSL method. 

# Requirements

All of the quantum simulations of the FSL method in this repository are performed using [Qiskit](https://qiskit.org/) (version 0.41.0).

The code in this repository also made use of the following Python packages:

* [NumPy (version 1.23.1)](https://pypi.org/project/numpy/1.23.1/) 
* [Matplotlib (version 3.5.2)](https://pypi.org/project/matplotlib/3.5.2/)
* [Graycode (version 1.0.5)](https://pypi.org/project/graycode/)

# How to cite this code

<a id="how-to-cite-this-code"></a> If you use Fourier Series Loader (FSL) method in your work, please consider citing the following paper:

> Mudassir Moosa, Thomas W. Watts, Yiyou Chen, Abhijat Sarma, and Peter L. McMahon (2023). Linear-depth quantum circuits for loading Fourier approximations of arbitrary functions. *Manuscript in preparation.*  


# License

The code in this repository is released under the following license:

- [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/)

A copy of this license is given in this repository as [license.txt](https://github.com/mcmahon-lab/Fourier-Series-Loader/blob/master/license.txt)
