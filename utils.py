
import numpy as np
import csv
import json
import bisect
import os

import matplotlib.pyplot as plt

import event_stream



def read_file(filename, start_time, end_time, start_idx=0, keep_only_off_events=False):

    decoder = event_stream.Decoder(filename)
    packets = []
    
    for packet in decoder:
        packets.append(packet)

    events = np.concatenate(packets)
    # NOTE FLIPPED x, y
    height = decoder.width
    width = decoder.height

    # Extract the details of each event, NOTE FLIPPED x, y
    x_c = events[start_idx:]['y']
    y_c = events[start_idx:]['x']
    p = events[start_idx:]['on']
    t = events[start_idx:]['t']
    
    # Apply start and end time filters
    mask = (t >= start_time) & (t <= end_time)
    x_c, y_c, p, t = x_c[mask], y_c[mask], p[mask], t[mask]

    #remove events outside of range (in event of bad data)
    p[(x_c >= width) | (x_c < 0) | (y_c >= height) | (y_c < 0)] = 0

    # if not keep_only_off_events:
    x_on, y_on, t_on = x_c[p==1], y_c[p==1], t[p==1]
    x_off, y_off, t_off = x_c[p==0], y_c[p==0], t[p==0]

    events_dict = {'e_x': x_c, 'e_y': y_c, 'e_p': p, 'e_t': t, 
                   'e_x_on': x_on, 'e_y_on': y_on, 'e_t_on': t_on,
                   'e_x_off': x_off, 'e_y_off': y_off, 'e_t_off': t_off}
    
    return events_dict


# read the lens log csv
def load_csv_to_dict(csv_file_path):
    dict = {'l_acg_t':[], 'l_acg_dpt':[]}

    with open(csv_file_path, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            dict['l_acg_t'].append(float(row['time'])*1e6) #convert from s to us
            dict['l_acg_dpt'].append(float(row['amplitude']))

    return dict

def load_json_to_dict(json_file_path):
    with open(json_file_path, 'r') as json_file:
        dict = json.load(json_file)

    return dict

# read the experiment data json

def get_files_from_experiment_folder(experiment_main_folder, experiment_name):

    experiment_folders = []
    individual_experiment_data_name = []

    # Walk through the directory and find folders with the specified prefix
    for root, dirs, files in os.walk(experiment_main_folder):
        for folder in dirs:
            if folder.startswith(experiment_name):
                experiment_folder = f'{os.path.join(root, folder)}/'
                experiment_folders.append(experiment_folder)

                # Search for files within the folder with the specified prefix
                for file in os.listdir(experiment_folder):
                    if '.es' in file:
                        individual_experiment_data_name.append(os.path.splitext(file)[0])  # Remove file format extension
                        break

        return experiment_folders, individual_experiment_data_name

def calc_event_rate(t, bin_width=1e3):
    delta_t = 0
    previous_interval_t = 0

    event_rates = {"rate":[], "t":[]}
    n_events = 0

    for t_curr in t:
        n_events += 1
        delta_t = (t_curr - previous_interval_t)

        if delta_t > bin_width:
            previous_interval_t = t_curr
            event_rates['rate'].append(n_events)
            event_rates['t'].append(t_curr)
            n_events = 0

    return event_rates


def calc_event_rate_stats(events, bin_width):
    
    event_stats = {}
    event_rates_per_bin = calc_event_rate(events['e_t'], bin_width)
    event_stats['t_global'] = event_rates_per_bin['t']
    event_stats['event_rates_global'] = event_rates_per_bin['rate']
    event_stats['mean_event_rates_global'] = np.mean(event_rates_per_bin['rate'])
    event_stats['std_event_rates_global'] = np.std(event_rates_per_bin['rate'])

    event_rates_per_bin = calc_event_rate(events['e_t_on'], bin_width)
    event_stats['t_on'] = event_rates_per_bin['t']
    event_stats['event_rates_on'] = event_rates_per_bin['rate']
    event_stats['mean_event_rates_on'] = np.mean(event_rates_per_bin['rate'])
    event_stats['std_event_rates_on'] = np.std(event_rates_per_bin['rate'])

    event_rates_per_bin = calc_event_rate(events['e_t_off'], bin_width)
    event_stats['t_off'] = event_rates_per_bin['t']
    event_stats['event_rates_off'] = event_rates_per_bin['rate']
    event_stats['mean_event_rates_off'] = np.mean(event_rates_per_bin['rate'])
    event_stats['std_event_rates_off'] = np.std(event_rates_per_bin['rate'])

    return event_stats

def get_dpt_per_event_t(t1, t2, L_ACG_positions):
    closest_indexes = [bisect.bisect_left(t2, val) for val in t1]
    dpt_per_event = [L_ACG_positions[index-1] for index in closest_indexes]

    return dpt_per_event


def plot_event_rate_vs_L_ACG(event_stats, event_rate_bin_width, L_ACG_positions, plot_both_polarities=False, plot_L_ACG=False):

    plt.figure()
    plt.title('Event Rate vs L_ACG Optical Power')
    if plot_both_polarities:
        plt.plot(event_stats['t_on'], event_stats['event_rates_on'], color = 'magenta', label='On Polarity Event Rate')
        plt.plot(event_stats['t_off'], event_stats['event_rates_off'], color = 'cyan', label='Off Polarity Event Rate')
    else:
        plt.plot(event_stats['t_global'], event_stats['event_rates_global'], color = 'magenta', label='Event Rate')
    plt.xlabel('Time (us)')
    plt.ylabel(f'Event rate (bin width: {event_rate_bin_width}us)')
    plt.legend(loc='upper right')

    if plot_L_ACG:
        ax2 = plt.gca().twinx()
        ax2.plot(L_ACG_positions['t'], L_ACG_positions['dpt'], color = 'black', label='L_ACG Optical Power')
        ax2.set_ylabel('L_ACG optical power (diopter, 1/f)')
        ax2.legend(loc='upper left')

    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.5)
    plt.show()
    
def ensure_folder_exists(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Folder created at: {path}")
    else:
        print(f"Folder already exists at: {path}")
