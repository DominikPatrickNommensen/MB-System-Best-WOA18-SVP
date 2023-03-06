# MB-System-Best-WOA18_SVP
MB-System tool to download World Ocean Atlas 2018 (WOA18) temperature and salinity files, calculate the sound velocity using the UNESCO and the Leroy formula, applying them on a flat seafloor swathfile and find the best World Ocean Atlas 2018 sound velocity profiles for sound speed correction of a survey.

## General idea and preparation

## Programs of this tool

| Program                 | Description
|-------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| mbdownloadwoa18         | Extracts WOA18 temperature and salinity for a specified area, calculates the UNESCO and Leroy sound velocites for each grid cell and standard depth and stores them as netcdf files together with an overview map of the specified area extracted.
| mbbestsvp               | Calculates/extracts several SVP files based on the date and extent of a swathfile (flat seafloor), by either using the mean of 4 surrounding profiles or by using the closest grid point. All SVPs are applied and their capability of flattening the seafloor is evaluated and ranked against one another, including the SVPs stored in the swathfile itself.


### mbdownloadwoa18


### mbdownloadwoa18

## Requirements

### MB-System
Since this tool utilizes MB-System and is meant to work in conjuction with other MB-System programs one should download and install MB-System. Please see also [dwcaress/MB-System](https://github.com/dwcaress/MB-System) from which also the below information was copied from:

The primary source of information about MB-System is the project website at [https://www.mbari.org/products/research-software/mb-system/](https://www.mbari.org/products/research-software/mb-system/), which includes sections on:

- [FAQ](https://www.mbari.org/products/research-software/mb-system/mb-system-faq/)
- [Download and installation](https://www.mbari.org/products/research-software/mb-system/how-to-download-and-install-mb-system/)
- [Documentation](https://www.mbari.org/products/research-software/mb-system/mb-system-documentation/)
- [User and developer email discussion lists](https://www.mbari.org/products/research-software/mb-system/mb-system-discussion-lists/)

### Python modules
