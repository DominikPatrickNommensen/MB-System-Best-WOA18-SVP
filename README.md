# MB-System-Best-WOA18_SVP
MB-System tool to download World Ocean Atlas 2018 (WOA18) temperature and salinity files, calculate the sound velocity using the UNESCO and the Leroy formula, applying them on a flat seafloor swathfile and find the best World Ocean Atlas 2018 sound velocity profiles for sound speed correction of a survey.

## MB-System Programs
--

| Program                 | Description
|-------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| mbdownloadwoa18         | Extracts WOA18 temperature and salinity for a specified area, calculates the UNESCO and Leroy sound velocites for each grid cell and standard depth and stores them as netcdf files together with an overview map of the specified area extracted.
| mbbestsvp               | Calculates/extracts several SVP files based on the date and extent of a swathfile (flat seafloor), by either using the mean of 4 surrounding profiles or by using the closest grid point. All SVPs are applied and their capability of flattening the seafloor is evaluated and ranked against one another, including the SVPs stored in the swathfile itself.
