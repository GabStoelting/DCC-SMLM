<!-- PROJECT LOGO -->
<br />
<div align="center">
  </a>

  <h3 align="center">DCC-SMLM</h3>

  <p align="center">
    A python library to help in determining oligomeric states of plasma membrane protein complexes using single molecule localization microscopy (SMLM; e.g. PALM or STORM) images
    <br />
    ·
    <a href="https://github.com/GabStoelting/DCC-SMLM/issues">Report Bug</a>
    ·
    <a href="https://github.com/GabStoelting/DCC-SMLM/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#citations">Citations</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

This project is an optimized implementation of the DCC-SMLM algorithm published in Tan et al. <br> The general idea is
that the oligomeric state of protein complexes can be determined from counting the number of colocalized signals when
the proteins of interest are labeled with two independent fluorophores (e.g. fluorescent proteins). <br>
As of now, work on this project is still ongoing and proper documentation is currently still missing.
However, the accompanying example notebooks should be sufficient to understand and implement the algorithm using
your own data. Please contact us for further help and feature requests via the GitHub issue tracker above.
<p align="right">(<a href="#top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

This library is written in python and will likely work with any version of python > 3.6. I recommend using <a href="https://www.anaconda.com/products/individual">Anaconda</a> or <a href="https://docs.conda.io/en/latest/miniconda.html">Miniconda</a> but this is no hard requirement. However, for the installation instructions below, I'll assume that you are able to use the conda package manager.

### Prerequisites

This library and the example notebooks require python (likely>3.6), jupyter notebook, numpy, pandas, matplotlib and seaborn libraries. 

### Installation

Please follow these steps:

1. Clone the repository using command line tools (below) or by using a GUI version of Git
   ```sh
   git clone https://github.com/GabStoelting/DCC-SMLM
   ```
2. Install the required libraries via conda if not installed already
   ```sh
   conda install numpy
   conda install notebook
   conda install pandas
   conda install matplotlib
   conda install seaborn
   ```

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

At first, you may want to look at the example analysis notebook files using example data from our projects. Later, you can copy the DCCSMLM.py file to your
data folders and implement your own analyses. Start an anaconda prompt or terminal session and start jupyter notebook

```sh
jupyter notebook
```

We recommend to explore the notebooks in the following order:
1. CA_Example.ipynb - This notebooks demonstrates how to extract the chromatic aberration from recordings of 
fluorescent micro beads.
2. SciH5_Example.ipynb or SMAP_Example.ipynb - These notebooks show a typical processing pathway for localizations 
extracted using SNSMIL or SMAP software. It includes correction of the chromatic aberration (determined as shown in 
CA_Example.ipynb) and sample drift. Furthermore it extracts the colocalization ratio from a rectangular region of
interest.
3. Calibration_Example.ipynb - In order to determine the oligomeric state of unknown proteins, it is necessary
to first determine the colocalization ratios of known protein complexes. An example of this procedure is contained within this file.
4. Example_POI_analysis.ipynb - This notebook demonstrates how to analyze colocalization ratios from several recordings combined.
Both, the coefficient of mismatch as well as a Kolmogorov-Smirnov test are performed for statistics.

In general, you may use these notebooks also for data analysis. In this case, we recommend to create a copy of the SciH5_Example.ipynb
or SMAP_Example.ipynb for every recording and combine the output into CSV files that can be used by the other analysis notebooks.
<p align="right">(<a href="#top">back to top</a>)</p>




<!-- ROADMAP -->
## Roadmap

- [ ] Implement background correction
- [ ] Create documentation

See the [open issues](https://github.com/GabStoelting/DCC-SMLM/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- CITATIONS -->
## Citations

This tool has been used in the following publications:
<ol>
  <li>
    <a href="https://www.biorxiv.org/content/10.1101/2021.10.05.463155v4">
      Tan, H.L. et al. Determination of protein stoichiometries via dual-color colocalization with single molecule localization microscopy. BioRxiv (2022).
    </a>
  </li>
</ol>



<!-- LICENSE -->
## License

Distributed under the GPL 3.0 License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Gabriel Stölting - gabriel.stoelting@bih-charite.de
Project Link: [https://github.com/GabStoelting/DCC-SMLM](https://github.com/GabStoelting/DCC-SMLM)

<p align="right">(<a href="#top">back to top</a>)</p>
