# Galaxy Builder catalogue column description
| Name | Data type | Description |
|:--|:--|:--|
| `dr7objid` | int64 | The SDSS DR7 object id of this galaxy (not unique) |
| `galaxy_builder_id` | int64 | The Zooniverse id associated with this galaxy builder galaxy (unique) |
| `is_validation` | bool | Whether this galaxy is a repeat in the validation subset|
| `ra` | float64 | Galaxy right-ascension|
| `dec` | float64 | Galaxy Declination |
| `chisq` | float64 | Reduced chi-squared of the model fit |
| `pogson_magnitude` | float64 | The fit r-band Pogson magnitude |
| `flux` | float64 | Total r-band flux of the galaxy (nanomaggies) |
| `disk_Re` | float64 | Disc effective radius (arcseconds) |
| `disk_Re_err` | float64 | Clustering error in disc effective radius (arcseconds) |
| `disk_q` | float64 | Disk ellipticity |
| `disk_q_err` | float64 | Clustering error in disc ellipticity |
| `centre_dx` | float64 | Centre position for the bulge and bar |
| `centre_dy` | float64 | Centre position for the bulge and bar |
| `bulge_Re` | float64 | Bulge effective radius (arcseconds) |
| `bulge_Re_err` | float64 | Clustering error in bulge effective radius (arcseconds) |
| `bulge_q` | float64 | Bulge ellipticity |
| `bulge_q_err` | float64 | Clustering error in bulge ellipticity |
| `bulge_n` | float64 | Bulge Sersic index |
| `bulge_fraction` | float64 | Bulge flux fraction relative to the galaxy as a whole |
| `bar_Re` | float64 | Bar effective radius |
| `bar_Re_err` | float64 | Clustering error in bar effective radius |
| `bar_q` | float64 | Bar ellipticity |
| `bar_q_err` | float64 | Clustering error in bar ellipticity |
| `bar_n` | float64 | Bar Sersic index |
| `bar_c` | float64 | Bar boxyness |
| `bar_fraction` | float64 | Bar flux fraction relative to the galaxy as a whole |
| `pitch_angle_0` | float64 | Pitch angle of spiral arm (degrees) |
| `pitch_angle_1` | float64 | Pitch angle of spiral arm (degrees) |
| `pitch_angle_2` | float64 | Pitch angle of spiral arm (degrees) |
| `pitch_angle_3` | float64 | Pitch angle of spiral arm (degrees) |
| `n_spirals` | float64 | Number of spirals in the fit model |
| `spiral_fraction` | float64 | Spiral flux fraction relative to the galaxy as a whole |
