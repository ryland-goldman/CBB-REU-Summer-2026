<a id="dataanalysis-visualization-paraview"></a>

# 3D Visualization: ParaView

WarpX results can be visualized by ParaView, an open source visualization and analysis software.
ParaView can be downloaded and installed from httpshttps://www.paraview.org.
Use the latest version for best results.

## Tutorials

ParaView is a powerful, general parallel rendering program.
If this is your first time using ParaView, consider starting with a tutorial.

* [https://www.paraview.org/Wiki/The_ParaView_Tutorial](https://www.paraview.org/Wiki/The_ParaView_Tutorial)
* [https://www.youtube.com/results?search_query=paraview+introduction](https://www.youtube.com/results?search_query=paraview+introduction)
* [https://www.youtube.com/results?search_query=paraview+tutorial](https://www.youtube.com/results?search_query=paraview+tutorial)

## openPMD

WarpX’ openPMD files can be visualized with ParaView 5.9+.
ParaView supports ADIOS1, ADIOS2 and HDF5 files, as it implements (like WarpX) against [openPMD-api](https://github.com/openPMD/openPMD-api).

For openPMD output, WarpX automatically creates an `.pmd` file per diagnostics, which can be opened with ParaView.

#### WARNING
[WarpX issue 21162](https://github.com/BLAST-WarpX/warpx/issues/1803):
We currently load WarpX field data with a rotation.
Please apply rotation of `0 -90 0` to mesh data.

#### WARNING
[ParaView issue 21837](https://gitlab.kitware.com/paraview/paraview/-/issues/21837):
In order to visualize particle traces with the `Temporal Particles To Pathlines`, you need to apply the `Merge Blocks` filter first.

If you have multiple species, you may have to extract the species you want with `Extract Block` before applying `Merge Blocks`.

## Plotfiles (AMReX)

ParaView also supports visualizing AMReX plotfiles.
Please see [the AMReX documentation](https://amrex-codes.github.io/amrex/docs_html/Visualization.html#paraview) for more details.

## In Situ Analysis with Catalyst 2

Continue reading [here](catalyst.md#visualization-catalyst).
