# GAP

This is a set of scripts which attempts to perform automated self-calibration of an ASKAP beam, the end result being a full-band Stokes-I image. These scripts are somewhat dusty and feral, I haven't got around to bringing them up-to-date and making them for public consumption, but hopefully they'll be of some use. 

The starting point is assumed to be a flagged and bandpass-calibrated / flux-scaled Measurement Set. From there the scripts perform three cycles of imaging and self-calibration. The first two are traditional (2GC) passes, one with phase corrections derived from an initial image of the data, and second with amplitude and phase corrections, with an intermediate imaging step to refine the model. 

The third pass of self-calibration is direction-dependent (3GC), and includes a step which attempts to identify problem sources in the field: those with high residual artefacts assumed to be present due to beam-to-beam differences and pointing errors. This is based on a measurement of the local pixel RMS around each component, using a thresholded comparison to the median RMS of all components, and weighted by the radial separation of the component from the nominal beam centre (things further out tend to cause more trouble as they sit in regions where the primary beam gradient is higher). The final step of self-cal then uses a hybrid sky model consisting of component models fitted to the identified problem sources, and a re-predicted visibility model based on the clean components with the problem sources excluded. Differential gains are then applied to the problem sources and the residual visibilities are then imaged, with the entire model restored back into the residual map at the final step.

I found this procedure to work quite well, producing maps with artefacts either eliminated or greatly reduced, however like just about everything that tries to automate self-cal it will probably trip up on some edge cases (for example, incomplete deconvolution of an extended source would leave a higher RMS in the residual map at that location, which might be interpreted as residual calibration errors). But then the intention here is really to automatically deal with beams that are edge-cases for regular direction-independent self-cal. There should probably also be more tuning / optimisation of the solution intervals than I've had time to do. Most of the hard-wired parameters should be sensible for ASKAP-12, although as more longer baselines get added the imaging parameters will need adjusting.

**The `run.py` script strings all the steps together, and in principle if you setup all the dependencies and run this in a folder containing a single-beam MS it should run to completion.** Note that it just picks up the MS in the folder that matches `*wtspec.ms`. I use CASA's mstransform tasks to add as WEIGHT_SPECTRUM column to the MS, hence the suffix. This gives much more flexible handling of visibility weights during the subsequent processing, and it's definitely worth doing.

Dependencies are as follows. *By far* the easiest way to deal with this lot is to use Ubuntu and make use of the excellent [KERN](http://kernsuite.info/) project.

* [wsclean](https://sourceforge.net/p/wsclean/wiki/Home/)

* [CubiCal](http://cubical.readthedocs.io/en/latest/): Note that the parsets provided are a bit out of date, and might have some parameters that are deprecated in the latest version. In particular it's worth utilising the madmax flagging algorithm in the latest version which seems to do a great job of removing bad visibilities during calibration. The latest versions of CubiCal can export the gain tables in CASA format, but I don't think that means you can simply utilise them with CASA's `applycal` task, certainly not for the directional solutions. It does make plotting the solutions simpler though, but there's a command line utility [here](https://github.com/IanHeywood/plot_utils) that will plot CubiCal's native parmdb format.

* [PyBDSF](http://www.astron.nl/citt/pybdsf/)

* [Pyrap](https://github.com/casacore/python-casacore): or python-casacore as it's now known.

* [Pyxis](https://github.com/ska-sa/pyxis)

* [Tigger](https://github.com/ska-sa/tigger): The FITS viewer isn't strictly required but the Python API for manipulating sky models is.

* PyFITS: If you're using astropy.io.fits instead (which is probably wise) then `import pyfits` can just be replaced with `from astropy.io import fits as pyfits`. 

* astLib: for astCoords and astWCS.

* [Montage](http://montage.ipac.caltech.edu/): I think this is only used if you want to visualise the model and the assumptions that have been made about the sources that have dE's applied to them (which you should definitely do as a sanity check). That's the `show_models.py` script which presently has to be executed manually.

Here's an example:

![](https://i.imgur.com/fV8pzw2.gif)

This shows the map obtained following the second (a&p) self-cal step, blinking with the final map. The problem source on the flank of the beam is recovered very cleanly, revealing a handful of sources that are not visible in the initial image. The RMS noise in the map centre improves by about 40%.

*"AIPS comes with absolutely no warranty."*
