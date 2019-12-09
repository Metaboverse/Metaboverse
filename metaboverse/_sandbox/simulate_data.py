import json
import matplotlib
import matplotlib.pyplot as plt
cmap = matplotlib.cm.get_cmap('seismic')
from random import randrange
import random
import numpy as np

dir = '/Users/jordan/Desktop/'
mu, sigma = 0, 1

def extract_value(
        expression_value,
        max_value):

    position = (expression_value + max_value) / (2 * max_value)
    rgba = cmap(position)

    return rgba

def convert_rgba(
        rgba_tuple):

    rgba_list = list(rgba_tuple)

    rgba_new = []
    for x in rgba_list[:3]:
        rgba_new.append(int(x * 255))

    rgba_new.append(rgba_list[3])
    return tuple(rgba_new)


with open(dir + 'HSA_global_reactions.json') as json_file:
    data = json.load(json_file)

len(data[0]['nodes'])


for x in range(len(data[0]['nodes'])):

    # if no value already there
    if data[0]['nodes'][x]['expression']['0'] == 'None':

        # Make random number
        # If odd, add value to graph
        #if randrange(10) % 2 == 0:

        # Make random number and get rgba
        # add to graph
        value = np.random.normal(mu, sigma, 1)[0]
        rgba = extract_value(
            expression_value=value,
            max_value=6)
        rgba_js = convert_rgba(
            rgba_tuple=rgba)

        data[0]['nodes'][x]['expression']['0'] = value
        data[0]['nodes'][x]['rgba']['0'] = rgba
        data[0]['nodes'][x]['rgba_js']['0'] = rgba_js



with open(dir + 'HSA_global_reactions_all_filled.json', 'w') as f:
    json.dump(data, f, indent=4) # Parse out as array for javascript
