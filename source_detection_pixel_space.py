
import matplotlib.pyplot as plt

from skimage import measure
from matplotlib.patches import Ellipse
from matplotlib.colors import to_rgba

import numpy as np



def display_img_sources(regions, intensity_img, img_mask, title='Detected Sources', scaling_factor=3, colour_source_overlay=False, save_fig=False):
    
    fig, ax = plt.subplots(figsize=([20, 10]))
    ax.set_title(title)

    for region in regions:
        #filter single pixel regions
        if region['intensity_max'] > 1:
            # Draw ellipse around the centroid with scaling factor
            centroid = region.centroid
            minr, minc, maxr, maxc = region.bbox
            height = (maxr - minr) * scaling_factor
            width = (maxc - minc) * scaling_factor

            # Overlay pixels with unique alpha color
            alpha_color = to_rgba(f"C{region.label}", alpha=0.8)
            
            if colour_source_overlay:
                # Create a mask for the region and Plot only the pixels within the region
                mask = img_mask[minr:maxr, minc:maxc] == region.label
                ax.plot(np.where(mask)[1] + minc, np.where(mask)[0] + minr, 's', color=alpha_color, markersize=3)

            #label each source with an ellipse scaled by the scale factor
            ellipse = Ellipse((centroid[1], centroid[0]), width, height, edgecolor=alpha_color, linewidth=1, fill=False, alpha=0.75)
            ax.add_patch(ellipse)


    # display the original image with labeled regions outlined
    ax.imshow(intensity_img, cmap='viridis', interpolation='nearest', vmin=0, vmax=10)


#detect sources in target frequency:
def detect_sources(input_img):

    print('Detecting sources in image frame using masking and region props')


    #do we try DBSCAN

    # detect sources using region props
    input_img_mask = input_img
    input_img_mask = (input_img_mask>0).astype(int)
    input_img_mask[input_img_mask >0] == 1
    input_img_mask_labelled = measure.label(input_img_mask)

    #detect sources as pixel regions covered by event mask
    regions = measure.regionprops(label_image=input_img_mask_labelled, intensity_image=input_img)

    #filter regions with an area of only 1 pixel or number of events at 1 pixel
    filtered_regions = [region for region in regions if (region['intensity_max'] > 1 and region['num_pixels'] > 1)]

    print(f'Removed {len(regions)-len(filtered_regions)} of {len(regions)} as single event sources with {len(filtered_regions)} sources remaining')

    #sort sources by the area in descending order
    sources = sorted(filtered_regions, key=lambda x: x['area'], reverse=True)

    #extract the region properties from the skimage.measure region props properties as list of dicts

    img_detected_sources = []

    for source in sources:
        properties_dict = {}
        for prop in source:
            properties_dict[prop] = source[prop]
        img_detected_sources.append(properties_dict)

    return img_detected_sources

