# Active_Neuromorphic_Space_Imaging_Dataset

### EBVS Liquid-Lens Focus Dataset using GEN 4 HD EVK4

Data and processing scripts for manuscript titled:

**"Active Neuromorphic Space Imaging and Focusing using Liquid Lenses"** 
by Nicholas Owen Ralph, Darren Maybour, Alexandre Marcireau, Nik Dennler, Sami Arja, Nimrod Kruger, and Gregory Cohen.

Affiliations for all authors: International Centre for Neuromorphic Systems (ICNS), MARCS Institute for Brain, Behaviour and Development, Western Sydney University

Manuscript is accessible on Arxiv: 

## Dataset:

The Active Neuromorphic Space Imaging Dataset can currently be accessed from: https://drive.google.com/drive/folders/1p6OVVBlNC_44ZNFSWqZTzlqz-CxaD6mG?usp=sharing

### Format:

FILENAME: *YYYY-MM-DDTHH-MM-SS.SSSSSSZ_outdoor_obs_DD_MM_YYY_field_FIELD_NUM_max_dpt-8.75_min_dpt-7.75_f-FREQUENCYHz*

Date is given by observation date, Field 1 M46 and Field 2 Orion. Max and min diopter is the same for all observations. Oscillation frequency given by 'frequency'. All observations are in seperate folders.

### Each observation folder contains:

- The diopter of the liquid lens timestamped by system clock during observation: *FILENAME_ideal_lens_pos_logs.csv* and *FILENAME_lens_pos_logs.csv* 
    - *ideal_los_pos_logs* contains diopter values of the liquid lens timestamped by the ideal time delta from the start of the recording to when the new amplitude is input to the controller (10 ms)
    - *lens_pos_logs* contain the actual time delta in microseconds from the start of the recording to when the new amplitude is input to the controller
- Event data in .raw and .es
- Accumulated events over the whole recording in ppm
- Video render of the recording in mp4
- Event data and observation log in json, containing liquid lens frequency, bias information, start time, end time.
- Plot comparing the input diopter to the on and off event rate throughout the duration of the recording in *FILENAME_L_ACG_pos_vs_event_rates.png*


## Processing Scripts

 - Main notebook for loading and processing observations from the dataset are found in *main.ipynb*
 - This file also runs the experiments described in the manuscript for both fields, including source detection, astrometric calibration and GAIA catalogue matching for event-based space imaging data
 - Display and plotting in the image and WCS frames in *plot_utils.py* and *plot_utils_astro.py*, respectively
 - Basic source finding and filtering in *source_detection_pixel_space.py* 
 - Astrometric calibration, catalogue matching and display using *astrometry_solver.py*  

#### Dependencies:
```sh
pip install tqdm astropy skyfield scikit-image
```

#### EBVS Specific Dependencies:
```sh
pip install event_stream
```


Contact Nicholas Owen Ralph at N.Ralph@westernsydney.edu.au for support. 
