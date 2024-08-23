
import matplotlib.pyplot as plt
import numpy as np

final_results_path = '../final_results/'


def plot_rates_vs_dpt(obs_data, frequencies, start_time, end_time, results_folder, results_filename):
    """
    Plots the on-event and off-event rates against time for each frequency, 
    along with L-ACG optical power (DPT). Optionally filters frequencies for 
    generating specific plots.
    """
    # Filtered frequencies for specific plots
    filtered_frequencies = [0, 0.03, 0.1, 1.0, 15.0]
    
    def plot_single_frequency(idx, singular_plot=False):
        """
        Helper function to plot the data for a single frequency.
        """
        fig, ax = plt.subplots(figsize=(9, 3))
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Event Rate (Eps)')
        ax.minorticks_on()
        ax.grid(True)
        
        # Plot on-event and off-event rates
        ax.plot(np.array(obs_data[idx]['t_on']) / 1e6, np.array(obs_data[idx]['event_rates_on']) / 1e6, color='orange', label='On-event rate')
        ax.plot(np.array(obs_data[idx]['t_off']) / 1e6, np.array(obs_data[idx]['event_rates_off']) / 1e6, color='purple', label='Off-event rate')
        ax.set_yscale('log')
        
        # X-axis limit based on frequency
        time_to_10_oscillations = 1 / (frequencies[idx] + 1e-6) * 60 * 1e6  # Time to 30 oscillations (in microseconds)
        obs_duration = obs_data[idx]['duration']
        x_limit_based_on_f = min(time_to_10_oscillations, obs_duration - (start_time + (obs_duration - end_time)))
        ax.set_xlim([start_time / 1e6, x_limit_based_on_f / 1e6])
        
        # Plot L-ACG optical power (DPT) on secondary y-axis
        ax2 = ax.twinx()
        ax2.set_ylabel('Optical Power (Diopter)')
        if frequencies[idx] == 0:
            ax2.hlines(y=8, xmin=0, xmax=max(obs_data[idx]['t_on']) / 1e6, color='grey', linestyle='solid', linewidth=2, alpha=0.5, label='L-ACG DPT')
        else:
            ax2.plot(np.array(obs_data[idx]['l_acg_t']) / 1e6, obs_data[idx]['l_acg_dpt'], color='grey', linestyle='solid', linewidth=2, alpha=0.5, label='L-ACG DPT')
        
        ax2.legend(loc='upper right' if singular_plot else 'lower right')
        ax.set_title(f'Frequency: {frequencies[idx]} Hz')
        ax.legend(loc='upper left')
        
        # Save the plot
        save_filename = f'{results_filename}_events_vs_L-ACG_pos' + (f'_f-{frequencies[idx]}.png' if singular_plot else '.png')
        plt.savefig(f'{results_folder}{save_filename}', bbox_inches='tight')
        plt.show()

    # Plot all frequencies or only filtered frequencies
    for idx, freq in enumerate(frequencies):
        if freq in filtered_frequencies:
            plot_single_frequency(idx, singular_plot=True)


def plot_event_rates(obs_data, frequencies, event_rate_bin_width, results_folder, results_filename, field_name):
    """
    Plots the event rates against frequencies for different observation data.
    Generates both subplots and individual plots based on the `run_as_subplot` flag.
    """

    run_as_subplot = False  # Set to True for subplots, False for individual plots

    def plot_mean_event_rates(ax, frequencies, obs_data, label_suffix=""):
        """
        Helper function to plot mean event rates.
        """
        ax.plot(frequencies, [np.array(inner_dict[f'mean_event_rates{label_suffix}_global'])/1e6 for inner_dict in obs_data.values()],
                '-+', label=f'Mean {label_suffix.strip()}event rate')
        ax.plot(frequencies, [np.array(inner_dict[f'mean_event_rates{label_suffix}_on'])/1e6 for inner_dict in obs_data.values()],
                '-+', label=f'Mean On{label_suffix} event rate')
        ax.plot(frequencies, [np.array(inner_dict[f'mean_event_rates{label_suffix}_off'])/1e6 for inner_dict in obs_data.values()],
                '-+', label=f'Mean Off{label_suffix} event rate')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(f'Event Rate (Eps)')
        ax.minorticks_on()
        ax.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.9)
        ax.legend()

    def plot_std_event_rates(ax, frequencies, obs_data, label_suffix=""):
        """
        Helper function to plot standard deviation of event rates.
        """
        ax.plot(frequencies, [np.array(inner_dict[f'std_event_rates{label_suffix}_on'])/1e6 for inner_dict in obs_data.values()],
                '-+', label=f'STD On{label_suffix} event rate')
        ax.set_xscale('log')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel(f'Event Rate (Eps)')
        ax.minorticks_on()
        ax.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.9)
        ax.legend()

    if run_as_subplot:
        fig, axs = plt.subplots(1, 2, figsize=(15, 6))
        plt.subplots_adjust(hspace=1)

        plot_mean_event_rates(axs[0], frequencies, obs_data)
        plot_std_event_rates(axs[1], frequencies, obs_data)

        plt.savefig(f'{results_folder}{results_filename}_event_stats_all_frequencies.png', bbox_inches='tight')
    else:
        fig, ax = plt.subplots(figsize=(9, 6))
        plot_mean_event_rates(ax, frequencies, obs_data)
        ax.set_title(f'Observation: {field_name}')
        plt.savefig(f'{results_folder}{results_filename}_event_stats_all_frequencies_means.png', bbox_inches='tight')
        plt.savefig(f'{final_results_path}{results_filename}_event_stats_all_frequencies_means.png', bbox_inches='tight')

        fig, ax = plt.subplots(figsize=(9, 6))
        plot_std_event_rates(ax, frequencies, obs_data)
        ax.set_title(f'Observation: {field_name}')
        plt.savefig(f'{results_folder}{results_filename}_event_stats_focussed.png', bbox_inches='tight')
        plt.savefig(f'{final_results_path}{results_filename}_event_stats_focussed.png', bbox_inches='tight')

    plt.show()




def plot_event_peaks(obs_data, frequencies, processed_peaks):

    #number of cols is 1, rows is given by the number of frequencies
    fig, axs = plt.subplots(len(frequencies), 1, figsize=(9, 3 * len(frequencies)))
    plt.subplots_adjust(hspace=0.5)

    for idx in range(len(frequencies)):
        ax = axs[idx]
        ax.set_title(f'Frequency: {frequencies[idx]}Hz')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel(f'EPus')
        ax.minorticks_on()
        ax.grid('on')
        ax.plot(np.array(obs_data[idx]['t_on'])/1e6, obs_data[idx]['event_rates_on'],color='magenta', label='On-event rate')
        diff_event_rate = np.diff(np.array(obs_data[idx]['event_rates_on'])/1e6)
        ax.plot(np.array(obs_data[idx]['t_on'][:-1])/1e6, diff_event_rate, color='cyan', label='Diff On-event rate')

        #find peaks in the event rate plots
        peaks = processed_peaks[idx]
        event_rate_peaks = [obs_data[idx]['event_rates_on'][peak] for peak in peaks]
        event_rate_peak_times = [obs_data[idx]['t_on'][peak] for peak in peaks]
        if peaks.size > 0:
            # Convert peaks to integers
            ax.plot(np.array(event_rate_peak_times)/1e6, event_rate_peaks, "x", color='k', label='Detected Peak')

        ax.legend(loc = 'upper left')

        #if frequency too high, cant show without aliasing and crowding 
        if frequencies[idx] < 2:
            ax2 = ax.twinx()

        if peaks.size > 0:
            for peak_time in event_rate_peak_times:
                dpt_at_peak_time = obs_data[idx]['l_acg_dpt'][np.abs(obs_data[idx]['l_acg_t'] - peak_time).argmin()]
                # ax2.axhline(y=dpt_at_peak_time, color='k', alpha=0.5)

def focus_frame_comparisons(obs_data, focus_frames, height, width, results_folder, results_filename):
    separated_frames = True
    filtered_frequencies = [0, 0.1, 5, 10]

    def plot_frame(frame, title, filename, vmin=0, vmax=1):
        fig, ax = plt.subplots(figsize=(10, 10))
        im = ax.imshow(frame, vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.set_xlabel('x (pix)')
        ax.set_ylabel('y (pix)')

        # Determine the actual range of the image data
        actual_vmax = np.max(frame)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink=0.8, aspect=10)
        cbar.set_ticks([0, vmax/2, vmax])
        cbar.set_ticklabels([0, actual_vmax/2, actual_vmax])
        cbar.set_label('Event Count')
        
        plt.savefig(f'{results_folder}{filename}.png', bbox_inches='tight')
        plt.show()

    for idx, (key, data) in enumerate(obs_data.items()):
        curr_frequency = data["frequency"]
        
        if not separated_frames or curr_frequency in filtered_frequencies:
            frame = np.zeros((height, width))
            on_frame = np.zeros((height, width))

            # Accumulate events
            np.add.at(frame, (data['e_x'], data['e_y']), 1)
            np.add.at(on_frame, (data['e_x_on'], data['e_y_on']), 1)

            # Focus frame handling
            focus_frame = focus_frames[idx] if curr_frequency != 0 else on_frame
            focus_frame_title = f'On-events in focus period at {curr_frequency}Hz' if curr_frequency != 0 else 'On-events in focus period (L-ACG off at 0 Hz)'

            # Plotting
            plot_frame(frame, f'Accumulated events {curr_frequency}Hz', f'{results_filename}_event_frame_all_events_{curr_frequency}Hz', vmax=3)
            plot_frame(on_frame, f'Accumulated on-events {curr_frequency}Hz', f'{results_filename}_event_frame_on_only_{curr_frequency}Hz', vmax=3)
            plot_frame(focus_frame, focus_frame_title, f'{results_filename}_event_frame_on_during_focus_{curr_frequency}Hz', vmax=3)


