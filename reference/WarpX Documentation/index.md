# WarpX

WarpX is an advanced **Particle-In-Cell** code.

It supports many features including:

> - Multiple types of field solvers (incl. for [Maxwell’s equations](theory/models_algorithms/electromagnetic_pic.md#theory-em-pic), Poisson’s equation, and [Ampere’s law coupled with Ohm’s law](theory/models_algorithms/kinetic_fluid_hybrid_model.md#theory-kinetic-fluid-hybrid-model))
> - Various grid geometries (1D/2D/3D Cartesian, cylindrical, spherical)
> - Multi-physics packages (incl. ionization, atomic, fusion and collisional physics, as well as quantum electrodynamics)
> - Advanced numerical methods (incl. explicit and implicit time advance, mesh refinement, boosted-frame simulations, embedded boundaries, pseudo-spectral solvers)

For details on these features, see the [theory section](theory/intro.md#theory).
WarpX has been applied to a wide variety of science projects, see [highlights](highlights.md#highlights).

In addition, WarpX is a *highly-parallel and highly-optimized code*:

> - Can run on multi-core CPUs as well as NVIDIA, AMD or Intel GPUs
> - Scales to the world’s largest supercomputers and includes load balancing capabilities. WarpX was awarded the [2022 ACM Gordon Bell Prize](https://www.exascaleproject.org/ecp-supported-collaborative-teams-win-the-2022-acm-gordon-bell-prize-and-special-prize/).
> - Multi-platform code that can run on Linux, macOS and Windows.
> - Can be run and [extended via its Python interface](usage/workflows/python_extend.md#usage-python-extend), e.g., to couple to other codes or AI/ML frameworks.

<a id="contact"></a>

## Contact us

The [WarpX GitHub repository](https://github.com/BLAST-WarpX/warpx) is the main communication platform:

> - If you are new to WarpX or have a question, we encourage you to visit our [discussions page](https://github.com/BLAST-WarpX/warpx/discussions) and connect with the community. This page is also a great place to browse answers to previously asked questions, post new ones, get help with installation, exchange ideas, and share feedback.
> - You can also explore the icons in the upper right corner of the [WarpX GitHub repository](https://github.com/BLAST-WarpX/warpx) (e.g., `Watch`, `Star`, etc.): feel free to watch the repository if you want to receive updates, or to star the repository to support the project.
> - For bug reports, feature requests, or installation issues, you can also open a new [issue](https://github.com/BLAST-WarpX/warpx/issues).
<style>
/\* front page: hide chapter titles
 \* needed for consistent HTML-PDF-EPUB chapters
 \*/
section#installation,
section#usage,
section#tutorials,
section#theory,
section#data-analysis,
section#development,
section#maintenance,
section#epilogue {
    display:none;
}
</style>

# Installation

<!-- install/changelog -->
<!-- install/upgrade -->

# Usage

# Tutorials

# Data Analysis

# Theory

# Development

<!-- good to have in the future: -->
<!-- developers/repostructure -->

# Maintenance

# Epilogue
