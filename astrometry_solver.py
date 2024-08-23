

import astropy
import astrometry
import numpy as np
import logging
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.wcs.utils import proj_plane_pixel_area
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.visualization.wcsaxes import WCSAxes
from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.stats import gaussian_kde
 
import plot_utils_astro

class AstrometrySolver:

    def __init__(self, experiment_name, results_filename, field_image, image_detected_sources, frequency, centre_ra, centre_dec, gaia_sources, field_name):
        
        self.centre_ra = centre_ra
        self.centre_dec = centre_dec
        self.gaia_sources = gaia_sources
        self.experiment_name = experiment_name
        self.frequency = frequency
        self.results_filename = results_filename
        self.field_image = field_image
        self.img_detected_sources = image_detected_sources
        self.field_name = field_name

        self.img_detected_source_positions = None

        self.solution = False
        self.wcs_calibration = None
        self.astrometric_solution_found = False

        #unique to astrometry net, just shows the sources it matched from the bright 20
        self.astrometry_detected_sources = {'pos_x':[], 'pos_y':[]} 
        self.sources_calibrated = []

        self.formatted_astrometry_results = None

    def astrometric_calibration(self):

        #number of sources to input for calibration, usually take the first 20 (brightest since sorted)
        number_sources_to_calibrate_from = min(30, len(self.img_detected_sources))
        
        self.img_detected_source_positions = [[float(source['centroid_weighted'][0]), float(source['centroid_weighted'][1])] for source in self.img_detected_sources]
        
        if number_sources_to_calibrate_from > 0:

            solver = astrometry.Solver(
                astrometry.series_5200_heavy.index_files(
                    cache_directory="/Users/30045541/Documents/data/astro/astrometry/astrometry_cache",
                    scales={5},
                )
            )
            self.solution = solver.solve(
                stars=self.img_detected_source_positions[0:number_sources_to_calibrate_from],
                size_hint=astrometry.SizeHint(
                    lower_arcsec_per_pixel=1.5,
                    upper_arcsec_per_pixel=2.0,
                ),
                position_hint=None,
                # position_hint=astrometry.PositionHint(
                #     ra_deg=self.centre_ra,
                #     dec_deg=self.centre_dec,
                #     radius_deg=1, #1 is good, 3 also works, trying 5
                # ),
                solution_parameters=astrometry.SolutionParameters(
                tune_up_logodds_threshold=None, # None disables tune-up (SIP distortion)
                output_logodds_threshold=1.0, #14.0 good, andre uses 1.0
                # logodds_callback=logodds_callback
                logodds_callback=lambda logodds_list: astrometry.Action.CONTINUE #CONTINUE or STOP
                )
            )
            
            if self.solution.has_match():
                self.astrometric_solution_found = True

                print('Solution found')
                print(f'Centre ra, dec: {self.solution.best_match().center_ra_deg=}, {self.solution.best_match().center_dec_deg=}')
                print(f'Pixel scale: {self.solution.best_match().scale_arcsec_per_pixel=}')
                print(f'Number of sources found: {len(self.solution.best_match().stars)}')

                self.wcs_calibration = astropy.wcs.WCS(self.solution.best_match().wcs_fields)
                pixels = self.wcs_calibration.all_world2pix([[star.ra_deg, star.dec_deg] for star in self.solution.best_match().stars],0,)

                for idx, star in enumerate(self.solution.best_match().stars):
                    self.astrometry_detected_sources['pos_x'].append(pixels[idx][0])
                    self.astrometry_detected_sources['pos_y'].append(pixels[idx][1])

            else:
                print(f'\n *** No astrometric solution found ***')

        else:
            print(f'\n *** No astrometric solution found, too few input sources ***')


    def run(self, if_display=False):

        print(f'Astrometrically calibrating field centred at {self.centre_ra}, {self.centre_dec}')

        #run astrometry
        self.astrometric_calibration()
        self.sources_calibrated = self.img_detected_sources.copy()

        #if we have an astrometric solution, format output and associate detected sources to gaia sources
        if self.astrometric_solution_found:

            #make an astro source entry for each source
            for idx in range(len(self.sources_calibrated)):
                centroid_radec_deg = self.wcs_calibration.all_pix2world([self.sources_calibrated[idx]['centroid_weighted']], 0,)
                self.sources_calibrated[idx]['centroid_radec_deg'] = centroid_radec_deg

            #find matching astrophysical sources for each detected source in the pixel space
            print('Searching GAIA for associated astrophyical sources')

            #need to search for every source, not just the ones from the brightest sources sublist or the matching sources sublist
            sources_ra = [source['centroid_radec_deg'][0][0] for source in self.sources_calibrated]
            sources_dec = [source['centroid_radec_deg'][0][1] for source in self.sources_calibrated]

            # Create a list of SkyCoord objects for detected source positions
            source_coords = SkyCoord(ra=sources_ra, dec=sources_dec, unit=(u.degree, u.degree), frame='icrs')
            gaia_source_coords = SkyCoord(ra=self.gaia_sources['ra'], dec=self.gaia_sources['dec'], unit=(u.deg, u.deg), frame='icrs')
            
            # Loop through each source and find the closest match locally
            match_idxs, d2d, _ = match_coordinates_sky(source_coords, gaia_source_coords, nthneighbor=1)
            duplicate_indices = np.where(np.bincount(match_idxs) > 1)[0]

            match_mags = []

            #update all sources with associated astrophyisical characteristics of the matching astrometric source
            for i, match_idx in enumerate(match_idxs):
                    closest_source = self.gaia_sources[match_idx]
                    position_pix = self.wcs_calibration.all_world2pix(closest_source['ra'], closest_source['dec'], 0)

                    self.sources_calibrated[i]['matching_source_ra'] = closest_source['ra']
                    self.sources_calibrated[i]['matching_source_dec'] = closest_source['dec']
                    self.sources_calibrated[i]['matching_source_x'] = (position_pix[0])
                    self.sources_calibrated[i]['matching_source_y'] = (position_pix[1])
                    self.sources_calibrated[i]['match_error_arcsec'] = d2d[i].arcsecond
                    self.sources_calibrated[i]['match_error_pix'] = d2d[i].arcsecond / self.solution.best_match().scale_arcsec_per_pixel
                    self.sources_calibrated[i]['mag'] = closest_source['phot_g_mean_mag']

                    match_mags.append(closest_source['phot_g_mean_mag'])

                    #check to make sure that the associated source is actually within the image, and store if so
                    if 0 < position_pix[0] < 720 and 0 < position_pix[1] < 1280:
                        self.sources_calibrated[i]['match_in_fov'] = True
                    else:
                        self.sources_calibrated[i]['match_in_fov'] = False

                    #flag duplicates 
                    if i in duplicate_indices:
                        self.sources_calibrated[i]['duplicate'] = 1
                    else:
                        self.sources_calibrated[i]['duplicate'] = 0
                        
                    self.sources_calibrated[i]['pixel_scale_arcsec'] = self.solution.best_match().scale_arcsec_per_pixel
                    self.sources_calibrated[i]['field_centre_radec_deg'] = [self.solution.best_match().center_ra_deg, self.solution.best_match().center_dec_deg]
                    self.sources_calibrated[i]['wcs'] = self.wcs_calibration
                    self.sources_calibrated[i]['detected_sources_pix'] = len(self.sources_calibrated)

            #get the most important values and save to disk
            matching_errors = [float(source['match_error_pix']) for source in self.sources_calibrated]
            matching_errors_arcsec = [float(source['match_error_arcsec']) for source in self.sources_calibrated]
            mags = [float(source['mag']) for source in self.sources_calibrated]
            areas = [float(source['equivalent_diameter_area']) for source in self.sources_calibrated]
            duplicate_flag = [source['duplicate'] for source in self.sources_calibrated]

            self.formatted_astrometry_results = {
                'experiment_name':self.experiment_name, 
                'frequency': self.frequency,
                'matching_errors': matching_errors,
                'mean_matching_errors': np.mean(matching_errors),
                'std_matching_errors':np.std(matching_errors), 
                'matching_errors_arcsec': matching_errors_arcsec,
                'mean_matching_errors_arcsec': np.mean(matching_errors_arcsec),
                'std_matching_errors_arcsec':np.std(matching_errors_arcsec), 
                'event_counts': [float(source['event_count']) for source in self.sources_calibrated],
                'event_rates': [float(source['event_rate']) for source in self.sources_calibrated],
                'mags': mags,
                'areas': areas,
                'astrometric_solution': True, 
                'detected_sources': len(self.sources_calibrated),
                'duplicate_flags': duplicate_flag,
                'matched_sources': len(mags),
                'matched_sources_within_fov':[source['match_in_fov'] for source in self.sources_calibrated],
                'wcs_calibration':self.wcs_calibration
            }

        else:
            self.formatted_astrometry_results = {'experiment_name':self.experiment_name,
                            'frequency':self.frequency,
                            'astrometric_solution':False, 
                            'detected_sources': len(self.img_detected_sources),
                            'matched_sources':0}

        if if_display and self.astrometric_solution_found:
            plot_utils_astro.display_results(self.field_image, 
                                             self.sources_calibrated, 
                                             self.formatted_astrometry_results, 
                                             self.wcs_calibration, 
                                             self.results_filename, 
                                             self.frequency,
                                             self.field_name,
                                             save_fig=True)

        print('***                  Done                ***')

        return self.sources_calibrated, self.formatted_astrometry_results, self.wcs_calibration
    

    