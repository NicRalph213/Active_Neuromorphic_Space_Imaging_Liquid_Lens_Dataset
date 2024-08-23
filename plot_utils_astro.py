

import numpy as np
from astropy import units as u
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.wcs.utils import proj_plane_pixel_area
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.visualization.wcsaxes import WCSAxes
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.stats import gaussian_kde
 

final_results_path = '../final_results/'


def display_results(field_image, sources_calibrated, output_results, wcs_calibration, results_filename, frequency, field, save_fig=False):

    #display calibration results
    #*******************************************************************************************
    results_filename = f'{results_filename}_{frequency}-hz'
    filtered_frequencies_flag = True
    filtered_frequencies = [0,0.1, 0.75, 1.25, 1.5]


    if frequency in filtered_frequencies or filtered_frequencies_flag is False:
        # Plotting the image with WCS axis
        fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': wcs_calibration})
        ax.set_title('Comparison of all detected sources (magenta) with associated astrophysical sources (cyan)')

        # Make the x-axis (RA) tick intervals and grid dense
        ax.coords[0].set_ticks(spacing=0.1*u.deg, color='white', size=6)
        ax.coords[0].grid(color='white', linestyle='solid', alpha=0.7)
        # Customize the y-axis (Dec) tick intervals and labels
        ax.coords[1].set_ticks(spacing=0.1*u.deg, color='white', size=6)
        ax.coords[1].grid(color='white', linestyle='solid', alpha=0.7)

        # Set the number of ticks and labels for the y-axis
        ax.coords[1].set_ticks(number=10)
        ax.set_xlabel('Right Ascention (ICRS deg)')
        ax.set_ylabel('Declination (ICRS deg)')

        # Display the image
        im = ax.imshow(field_image, origin='lower', cmap='viridis', vmin=0, vmax=3)

        # for source in sources_calibrated:
        #     circle = Circle((source['centroid'][1], source['centroid'][0]), color='magenta', radius=3, fill=False,alpha=0.25)
        #     ax.add_patch(circle)

        # for source in sources_calibrated:
        #     circle = Circle((source['matching_source_y'], source['matching_source_x']), color='cyan', radius=7.5, fill=False,alpha=0.25)
        #     ax.add_patch(circle)

        plt.tight_layout()
        if save_fig: plt.savefig(f'{results_filename}_astro_detected_vs_matched_from_gaia.png', bbox_inches='tight', pad_inches=0)
        plt.show()


        #*******************************************************************************************

        #scatter all sources detected and catalogued in image
        plt.figure(figsize=(12, 7))
        plt.title('Comparison of input brightest detected sources (magenta) with astometric calibrated and associated sources (cyan)')
        plt.xlabel('X (pix)')
        plt.ylabel('Y (pix)')
        plt.axvline(x=0, color='magenta', linestyle='--')
        plt.axvline(x=1280, color='magenta', linestyle='--')
        plt.axhline(y=0, color='magenta', linestyle='--')
        plt.axhline(y=720, color='magenta', linestyle='--')

        catalog_sources_x = [source['centroid'][0] for source in sources_calibrated]
        catalog_sources_y = [source['centroid'][1] for source in sources_calibrated]

        detected_sources_x = [source['matching_source_x'] for source in sources_calibrated]
        detected_sources_y = [source['matching_source_y'] for source in sources_calibrated]

        plt.scatter(catalog_sources_y, catalog_sources_x)
        plt.scatter(detected_sources_y, detected_sources_x)
        plt.grid(color='white', linestyle='solid', alpha=0.7)
        plt.tight_layout()
        if save_fig: plt.savefig(f'{results_filename}_astro_input_sources_vs_matched_astrometry.png', bbox_inches='tight', pad_inches=0)
        plt.show()

        #*******************************************************************************************

        plt.figure(figsize=(12, 7))
        plt.imshow(field_image, vmin=0, vmax=1)

        fig, ax = plt.subplots()

        # Display the image
        img = ax.imshow(field_image, vmin=0, vmax=1)

        # Add the colorbar, correctly linking it to the image and the axis
        cbar1 = fig.colorbar(img, ax=ax, fraction=0.046, pad=0.04)

        # Get the maximum value in the field_image for colorbar labels
        max_val_frame = np.max(field_image)

        # Set colorbar ticks and labels
        cbar1.set_ticks([0, 0.5, 1])
        cbar1.set_ticklabels([0, int(max_val_frame/2), int(max_val_frame)])
        cbar1.set_label('Event Count')

        # Set axis limits and labels
        ax.set_xlim([550, 750])
        ax.set_ylim([250, 450])
        ax.set_xlabel('X (pix)')
        ax.set_ylabel('Y (pix)')
        ax.set_title(f'{field} - {frequency} Hz')

        # Add circles to the plot
        for source in sources_calibrated:
            circle = Circle((source['centroid_weighted'][1], source['centroid_weighted'][0]), 
                            color='red', radius=3.5, fill=False, alpha=0.75, linewidth=1)
            ax.add_patch(circle)

        for source in sources_calibrated:
            circle = Circle((source['matching_source_y'], source['matching_source_x']), 
                            color='cyan', radius=5, fill=False, alpha=0.75, linewidth=1, linestyle='--')
            ax.add_patch(circle)

        # Add grid
        ax.grid(color='white', linestyle='dotted', alpha=0.25)

        # Adjust layout and save figure
        plt.tight_layout()

        if save_fig: 
            plt.savefig(f'{results_filename}_centre_crop-astro_detected_vs_matched_from_gaia.png', bbox_inches='tight', pad_inches=0)
        if save_fig and frequency in filtered_frequencies: 
            plt.savefig(f'{final_results_path}{field}_{frequency}-Hz_centre_crop-astro_detected_vs_matched_from_gaia.png', bbox_inches='tight', pad_inches=0)

        plt.show()

        #*******************************************************************************************

        plt.figure(figsize=(10, 6))
        plt.hist(output_results['matching_errors'], bins=50)
        plt.xlim([0,3])
        plt.ylabel('Occurances')
        plt.xlabel('Astrometric association error magnitude (pix)')
        plt.axvline(x=3, color='magenta', linestyle='--', label='Matching error cutt-off (3 pix)')
        plt.legend()
        plt.grid('on')
        plt.tight_layout()
        if save_fig: plt.savefig(f'{results_filename}_match_error_hist.png', bbox_inches='tight', pad_inches=0)
        plt.show()

        #*******************************************************************************************

        plt.figure(figsize=(10, 6))
        plt.title(f'{field} - {frequency} Hz')
        sc = plt.scatter(output_results['mags'], output_results['event_counts'], c=output_results['matching_errors'])
        plt.ylabel('Event Rate (Eps)')
        plt.xlabel('G-Band Magnitude')
        plt.axvline(x=14.45, color='magenta', linestyle='--', label='Previusly reported\nmag limit\n(Ralph et al. 2022)')
        plt.legend()
        plt.yscale('log')
        plt.grid('on')
        cbar = plt.colorbar(sc)
        cbar.set_label('Astrometric association error magnitude (pix)')
        plt.tight_layout()
        if save_fig: plt.savefig(f'{results_filename}_scatter_eps_vs_mag_cmap_match_error.png', bbox_inches='tight', pad_inches=0)
        if save_fig and frequency in filtered_frequencies: plt.savefig(f'{final_results_path}{field}_{frequency}-hz_scatter_eps_vs_mag_cmap_match_error.png', bbox_inches='tight', pad_inches=0)
        plt.show()

        #*******************************************************************************************

        plt.figure(figsize=(10, 6))
        # #calculate point density for display
        # kde = gaussian_kde([output_results['mags'], output_results['event_rates']])
        # density = kde([output_results['mags'], output_results['event_rates']])
        # normalized_density = (density - density.min()) / (density.max() - density.min())

        sc = plt.scatter(output_results['mags'], output_results['event_rates'])#, c=normalized_density)
        plt.ylabel('Event Rate (Eps)')
        plt.xlabel('G-Band Magnitude')
        plt.axvline(x=14.45, color='magenta', linestyle='--', label='Previusly reported\nmag limit\n(Ralph et al. 2022)')
        plt.legend()
        plt.grid('on')
        cbar = plt.colorbar(sc)
        cbar.set_label('Normalised Sample Density')
        plt.tight_layout()
        if save_fig: plt.savefig(f'{results_filename}_scatter_eps_vs_mag_cmap_density.png', bbox_inches='tight', pad_inches=0)
        plt.show()

        #*******************************************************************************************

        sc = plt.scatter(output_results['mags'], output_results['matching_errors'])#, c=normalized_density)
        plt.ylabel('Event Rate (Eps)')
        plt.xlabel('G-Band Magnitude')
        plt.axvline(x=14.45, color='magenta', linestyle='--', label='Previusly reported\nmag limit\n(Ralph et al. 2022)')
        plt.legend()
        plt.grid('on')
        cbar = plt.colorbar(sc)
        cbar.set_label('Normalised Sample Density')
        plt.tight_layout()
        if save_fig: plt.savefig(f'{results_filename}_scatter_eps_vs_mag_cmap_density.png', bbox_inches='tight', pad_inches=0)
        plt.show()





def display_source_counts_per_frequency(all_formatted_astro_results, results_folder, results_filename, field):
    #plot source count vs frequency, detected and associated

    detected_sources_pix = [observation['detected_sources'] for observation in all_formatted_astro_results]
    detected_sources_astro = [observation['matched_sources'] for observation in all_formatted_astro_results]
    frequencies = [observation['frequency'] for observation in all_formatted_astro_results]

    # 
    plt.figure(figsize=(6,4))
    plt.plot(frequencies, detected_sources_pix, marker='+')
    # plt.plot(frequencies, detected_sources_astro,  marker='+', label='Astrometrically Associated Sources')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Source Counts')
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.9)
    # plt.legend()
    plt.title(f'{field}')
    plt.savefig(f'{final_results_path}{results_filename}_detected_and_matches_sources_vs_frequency.png', bbox_inches='tight')
    # plt.xscale('log')
    plt.show()

    plt.figure(figsize=(6,4))

    plt.plot(frequencies, detected_sources_pix, marker='+')
    # plt.plot(frequencies, detected_sources_astro,  marker='+', label='Astrometrically Associated Sources')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Source Counts')
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.9)
    # plt.legend()
    plt.title(f'{field}')
    # plt.xlim([0,100])
    # plt.axis([0,16,0,500])
    plt.savefig(f'{final_results_path}{results_filename}_detected_and_matches_sources_vs_frequency_log.png', bbox_inches='tight')
    plt.xscale('log')
    plt.show()


def display_astrometric_error_per_freq(all_formatted_astro_results, results_folder, results_filename, field):
    #plot mean and std association error for each frequency


    match_errors= [observation['mean_matching_errors'] for observation in all_formatted_astro_results if observation['astrometric_solution']]
    match_errors_std= [observation['std_matching_errors']/2 for observation in all_formatted_astro_results if observation['astrometric_solution']]

    frequencies = [observation['frequency'] for observation in all_formatted_astro_results if observation['astrometric_solution']]

    # 
    plt.plot(frequencies, match_errors, marker='+', label='Mean Error')
    plt.plot(frequencies, match_errors_std, marker='o', label='STD Error')

    # plt.plot(frequencies, detected_sources_astro,  marker='+', label='Astrometrically Associated Sources')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Mean Astrometric Association Error (pix)')
    plt.title(f'{field}')
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.9)
    plt.legend()
    plt.xscale('log')
    plt.savefig(f'{final_results_path}{results_filename}_match_error_vs_frequency.png', bbox_inches='tight')

    plt.show()

    plt.errorbar(frequencies, match_errors, yerr=match_errors_std, marker='o')
    plt.xscale('log')
    plt.xlabel('L-ACG Frequency (Hz)')
    plt.ylabel('Mean Astrometric Association Error (pix)')
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.9)
    plt.savefig(f'{results_folder}{results_filename}_match_error_vs_frequency_eb.png', bbox_inches='tight')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

    plt.errorbar(frequencies, match_errors, yerr=match_errors_std, marker='o')
    plt.xlabel('L-ACG Frequency (Hz)')
    plt.ylabel('Mean Astrometric Association Error (pix)')
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.9)
    plt.savefig(f'{results_folder}{results_filename}_match_error_vs_frequency_eb_log.png', bbox_inches='tight')
    plt.show()

