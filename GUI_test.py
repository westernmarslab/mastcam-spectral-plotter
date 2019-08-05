# Begin imports
import tkinter as tk  # Tkinter is used for making GUIs
import tkinter.messagebox
from tkinter import filedialog
from PIL import ImageTk, Image  # PIL, or Pillow, is used to display images in Tkinter
import os  # Used to interact with the computer
import numpy as np  # Standard data analysis package
import pandas as pd  # Standard data analysis package
import matplotlib as mpl  # Standard plotting package
mpl.use('TkAgg')
from matplotlib import pyplot as pyplot  # Standard plotting package
from uncertainties import ufloat  # Used to propogate error
import copy  # Used to copy python objects
from sklearn.decomposition import PCA  # Used for PCA, tends to cause errors with Pyinstaller
import pickle  # Used to save python objcts as files, this package is awesome
import csv  # Used to read csvs
import sys  # Used to get system info, eg. mac or windows
import matplotlib.backends.backend_tkagg as tkagg  # Needed for freezing, has to be imported last

root = tk.Tk() # starts the root of the app


class datapoint:
    """
    The datapoint class represents the reflectance at a given filter in a Mastcam spectrum.
    It has three attributes: filter, reflectance, error, and eye.
    """
    def __init__(self, filter, reflectance, error, eye):
        self.filter = filter  # Filter, stored as wavelength eye: "527 R"
        self.reflectance = reflectance  # Reflectance value
        self.error = error  # Error in reflectance
        self.eye = eye
        if isinstance(self.reflectance, str):
            self.reflectance = self.reflectance.strip(' ')
        try: self.reflectance = float(reflectance) # If it's a float
        except ValueError:
            self.reflectance = float('nan')
        # If the datapoint has a bad reflectance, it is a nan, and will be picked up by the below function

    def has_data(self):
        '''Used to determine if the datapoint has an actual reflectance, or if the reflectance is not valid'''
        if pd.isnull(self.reflectance): return False
        if isinstance(self.reflectance, float): return True
        if self.reflectance.strip(' ') == '': return False
        return True


class Spectra:
    '''
    The Spectra object contains spectral data. When loading data from the databas, one spectra object is created for
    each ROI in the database. Spectra objects have two dictionary attributes: data and metadata. Data is a dictionary
    containing datapoint objects, one for each spectrum. The metadata attribute is a dictionary with strings as keys,
    that contain various metadata on the observation, and other information used for plotting.
    '''
    def __init__(self):
        self = self

    def __init__(self, data, metadata):
        self.data = data  # Dictionary holding spectral data
        self.metadata = metadata  # Dictionary holding metadata
        self.assign_metadata()  # Assigns additional metadata

    def assign_metadata(self):
        ''' This section assigns member formation and group based off of elevation, and sets a simple color scene. '''
        elev2member = {
            'Jura': [0,-4170],
            'Pettegrove Point': [-4170,-4210],
            'Blunts Point': [-4210,-4280],
            'Sutton Island': [-4280,-4370],
            'Karasburg': [-4370,-4420],
            'Hartmann\'s Valley': [-4420,-4435],
            'Pahrump Hills': [-4435,-4460]
        }
        elev2formation = {
            'Murray':[-4050,-4435],
            'Kimberly':[-4473,-4493],
            'Yellowknife Bay':[-4510,-4520]
        }
        self.metadata['Member'] = ''
        self.metadata['Formation'] = ''
        self.metadata['Group'] = ''
        if self.has_rover_elevation() or self.has_target_elevation():
            elevation = ''
            if self.has_target_elevation():
                if self.metadata['Target Elevation'] == float('inf') or \
                    self.metadata['Target Elevation'] == -float('inf'):
                    elevation = self.metadata['Target Elevation']
                else:
                    elevation = self.metadata['Rover Elevation']
            else:
                elevation = self.metadata['Rover Elevation']
            for name, elev in elev2member.items():
                high = elev[0]
                low = elev[1]
                if elevation >= low and elevation < high:
                    self.metadata['Member'] = name
            for name, elev in elev2formation.items():
                high = elev[0]
                low = elev[1]
                if elevation >= low and elevation < high:
                    self.metadata['Formation'] = name
            if elevation >= -4460:
                self.metadata['Group'] = 'Mount Sharp'
            else: self.metadata['Group'] = 'Bradbury'

        # Tina lithology overrides lithology from elevation
        if not self.metadata['Tina Member'] == '':
            self.metadata['Member'] = self.metadata['Tina Member']
        if not self.metadata['Tina Formation'] == '':
            self.metadata['Formation'] = self.metadata['Tina Formation']

        # Set simplified color scheme
        if self.metadata['Feature'] == 'Dusty Rock' or self.metadata['Feature'] == 'Nodule-rich Rock':
            self.metadata['SimpleRGB'] = (0, 0, 1)
            self.metadata['SimpleFeature'] = 'Dusty Rock'
        if self.metadata['Feature'] == 'Dust-Cleared Rock' or self.metadata['Feature'] == 'Broken Rock Face':
            self.metadata['SimpleRGB'] = (109 / 255, 225 / 255, 0)
            self.metadata['SimpleFeature'] = 'Dust-Cleared Rock'
        if self.metadata['Feature'] == 'Drill Tailings' or self.metadata['Feature'] == 'Dump Piles':
            self.metadata['SimpleRGB'] = (254 / 255, 0, 255 / 255)
            self.metadata['SimpleFeature'] = 'Drill Fines'
        if 'Soil' in self.metadata['Feature']:
            self.metadata['SimpleRGB'] = (212 / 255, 0, 42 / 255)
            self.metadata['SimpleFeature'] = 'Soil'
        if self.metadata['Feature'] == 'Other' or self.metadata['Feature'] == 'Veins':
            self.metadata['SimpleRGB'] = (255 / 255, 255 / 255, 0)
            self.metadata['SimpleFeature'] = 'Other'


    # These functions are used to get various band parameters
    def get_band_depth(self, left, middle, right):
        # Returns band depth
        if self.has_data(left) and self.has_data(middle) and self.has_data(right):  # Has to have all data
            lx = float(left.split(' ')[0])
            mx = float(middle.split(' ')[0])
            rx = float(right.split(' ')[0])
            ly = float(self.data[left].reflectance)
            my = float(self.data[middle].reflectance)
            ry = float(self.data[right].reflectance)
            return 1 - my/(((ly - ry) / (lx - rx) * (mx - lx)) + ly)  # Equation is from Melissa's Veins paper
        else:
            return 'invlaid'

    def get_band_depth_error(self, left, middle, right):
        # Gets band depth error
        if self.has_data(left) and self.has_data(middle) and self.has_data(right):
            lx = float(left.split(' ')[0])
            mx = float(middle.split(' ')[0])
            rx = float(right.split(' ')[0])
            ly = ufloat(self.data[left].reflectance, self.data[left].error)
            my = ufloat(self.data[middle].reflectance, self.data[middle].error)
            ry = ufloat(self.data[right].reflectance, self.data[right].error)
            bd = 1 - my/(((ly - ry) / (lx - rx) * (mx - lx)) + ly)
            return bd.s # return standard deviation of band depth
        else:
            return 'invlaid'

    def get_slope(self, left, right):
        # Get slope
        if self.has_data(left) and self.has_data(right):  # Has to have all data
            lx = float(left.split(' ')[0])
            rx = float(right.split(' ')[0])
            ly = float(self.data[left].reflectance)
            ry = float(self.data[right].reflectance)
            return ((ly - ry) / (lx - rx))
        else:
            return 'invalid'

    def get_slope_error(self, left, right):
        # Get error in slope
        if self.has_data(left) and self.has_data(right):
            lx = float(left.split(' ')[0])
            rx = float(right.split(' ')[0])
            ly = ufloat(self.data[left].reflectance, self.data[left].error)
            ry = ufloat(self.data[right].reflectance, self.data[right].error)
            return (((ly - ry) / (lx - rx))).s
        else:
            return 'invalid'

    def get_ratio(self, top, bottom):
        if self.has_data(top) and self.has_data(bottom):
            ty = float(self.data[top].reflectance)
            by = float(self.data[bottom].reflectance)
            return (ty / by)
        else:
            return 'invalid'

    def get_ratio_error(self, top, bottom):
        if self.has_data(top) and self.has_data(bottom):
            ty = ufloat(self.data[top].reflectance, self.data[top].error)
            by = ufloat(self.data[bottom].reflectance, self.data[bottom].error)
            return (ty / by).s
        else:
            return 'invalid'

    def has_data(self, nm):
        # Used to determine if a spectrum contains a certain filter
        if nm in self.data:
            if not self.data[nm].reflectance == '' and not self.data[nm].reflectance == ' ' and not pd.isnull(self.data[nm].reflectance):
                return True
        return False

    def get_reflectance(self, filter):
        # Returns the reflectance at a given filter
        if self.has_data(filter):
            return float(self.data[filter].reflectance)
        else:
            return 'invalid'

    def get_reflectance_error(self, filter):
        # Gets error in reflectance for a given filter
        if self.has_data(filter):
            return float(self.data[filter].error)
        else:
            return 'invalid'


    # The following funtions are used to determine if a spectrum has certain metadata
    def has_LTST(self):
        if pd.isnull(self.metadata['LTST']) or self.metadata['LTST'] == '':
            return False
        return True

    def has_target_elevation(self):
        if pd.isnull(self.metadata['Target Elevation']) or self.metadata['Target Elevation'] == '':
            return False
        return True

    def has_rover_elevation(self):
        if pd.isnull(self.metadata['Rover Elevation']) or self.metadata['Rover Elevation'] == '':
            return False
        return True

    def has_tau(self):
        if pd.isnull(self.metadata['Tau']) or self.metadata['Tau'] == '':
            return False
        return True

    def has_member(self):
        if self.metadata['Member'] == '':
            return False
        return True

    def has_formation(self):
        if self.metadata['Formation'] == '':
            return False
        return True

    def has_tina_member(self):
        if self.metadata['Tina Member'] == '':
            return False
        return True

    def has_tina_formation(self):
        if self.metadata['Tina Formation'] == '':
            return False
        return True

    def has_group(self):
        if self.metadata['Group'] == '':
            return False
        return True
    #All spectra have a float designations

    def has_ls(self):
        if pd.isnull(self.metadata['Ls']) or self.metadata['Ls'] == '':
            return False
        return True

    def has_lat(self):
        if pd.isnull(self.metadata['Lat']) or self.metadata['Lat'] == '':
            return False
        return True

    def has_long(self):
        if pd.isnull(self.metadata['Long']) or self.metadata['Long'] == '':
            return False
        return True

    def has_traverse(self):
        if pd.isnull(self.metadata['Traverse']) or self.metadata['Traverse'] == '':
            return False
        return True


class Extraction:
    '''Extraction objects are created when the user selects options and them presses "plot" or "make csv".
    Extraction objects simply hold all the user decisions as attributes. For example, it stores:
    X and Y axis choices, any relevant filters, whether data is R* corrected or not, data constraints,
    normalization point, etc.'''
    def __init__(self, App):
        self.x_choice = App.x_options.get(tk.ACTIVE)  # Slope, reflectance, band depth, tau, PC, etc.
        self.x_1 = App.x_expand1.get(tk.ACTIVE)  # Get filters
        self.x_2 = App.x_expand2.get(tk.ACTIVE)
        self.x_3 = App.x_expand3.get(tk.ACTIVE)
        self.x_pca = App.x_pca.get(tk.ACTIVE)  # Get PC selection
        self.y_choice = App.y_options.get(tk.ACTIVE)  # Slope, reflectance, band depth, tau, PC, etc.
        self.y_1 = App.y_expand1.get(tk.ACTIVE)  # Get filters
        self.y_2 = App.y_expand2.get(tk.ACTIVE)
        self.y_3 = App.y_expand3.get(tk.ACTIVE)
        self.y_pca = App.y_pca.get(tk.ACTIVE)  # Get PC selection
        self.conditions_list = App.condition_list.get(0,tk.END)
        self.conditions = {'Sols': [], 'LTST': '', 'Tau': [], 'Features': [], 'Elevation': [],
                           'Mcam': [], 'Units':[], 'Float':[], 'Ls':[], 'Latitude':[], 'Longitude':[],'Traverse':[]}
        self.highlight = {'Sols': [], 'LTST': '', 'Tau': [], 'Features': [], 'Elevation': [],
                          'Mcam': [], 'Units':[], 'Float':[], 'Ls':[], 'Latitude':[], 'Longitude':[],'Traverse':[]}
        self.analyze_conditions(self.conditions_list)
        self.r_star = App.r_star.get()
        self.norm_point = App.norm_point.get()
        self.propogate_error = App.propogate_error.get()
        self.color_type = App.color_type.get()
        self.geologic_members = App.geologic_members.get()
        self.user_spectra = App.user_spectra.get()
        self.our_spectra = App.our_spectra.get()

        self.set_labels()

        self.userhighlight = False
        for key, item in self.highlight.items():
            if key == 'LTST' and not item == '':
                self.userhighlight = True
            if key != 'LTST' and not item == []:
                self.userhighlight = True

    def set_labels(self):
        '''This function sets the plotting labels.'''
        if self.x_choice == 'Reflectance':
            self.x_label = 'Reflectance: ' + str(self.x_1)
        elif self.x_choice == 'Slope':
            self.x_label = 'Slope: ' + str(self.x_1) + ' to ' + str(self.x_2)
        elif self.x_choice == 'Ratio':
            self.x_label = 'Ratio: ' + str(self.x_1) + ' over ' + str(self.x_2)
        elif self.x_choice == 'Band Depth':
            self.x_label = 'Band Depth: ' + str(self.x_1) + '_' + str(self.x_2) + '_' + str(self.x_3)
        elif self.x_choice == 'Elevation':
            self.x_label = 'Elevation (m)'
        elif self.x_choice == 'PCA':
            self.x_label = self.x_pca
        else: self.x_label = self.x_choice

        if self.y_choice == 'Reflectance':
            self.y_label = 'Reflectance: ' + str(self.y_1)
        elif self.y_choice == 'Slope':
            self.y_label = 'Slope: ' + str(self.y_1) + ' to ' + str(self.y_2)
        elif self.y_choice == 'Ratio':
            self.y_label = 'Ratio: ' + str(self.y_1) + ' over ' + str(self.y_2)
        elif self.y_choice == 'Band Depth':
            self.y_label = 'Band Depth: ' + str(self.y_1) + '_' + str(self.y_2) + '_' + str(self.y_3)
        elif self.y_choice == 'Elevation':
            self.y_label = 'Elevation (m)'
        elif self.y_choice == 'PCA':
            self.y_label = self.y_pca
        else: self.y_label = self.y_choice

    def get_relevant_wavelengths(self):
        '''This method returns all the relevant wavelengths. For example, if the user wants to plot some band
        parameters, this function returns those bands. This is important for plotting purposes, because you can
        only plot observations that have all the bands needed for a given band parameter.'''
        ret = []
        if self.x_choice == 'Reflectance':
            ret.append(self.x_1)
        if self.x_choice == 'Slope' or self.x_choice == 'Ratio':
            ret.append(self.x_1)
            ret.append(self.x_2)
        if self.x_choice == 'Band Depth':
            ret.append(self.x_1)
            ret.append(self.x_2)
            ret.append(self.x_3)
        if self.y_choice == 'Reflectance':
            ret.append(self.y_1)
        if self.y_choice == 'Slope' or self.x_choice == 'Ratio':
            ret.append(self.y_1)
            ret.append(self.y_2)
        if self.y_choice == 'Band Depth':
            ret.append(self.y_1)
            ret.append(self.y_2)
            ret.append(self.y_3)
        if self.y_choice == 'PCA' or self.x_choice == 'PCA':
            for wavelength in ['445 L', '447 R', '493 R', '495 L', '527 L', '527 R','551 R', '554 L', '638 R',
                               '640 L', '676 L', '751 L', '805 R', '867 L', '908 R', '937 R', '1012 L', '1013 R']:
                ret.append(wavelength)
        return ret

    # The following functions are used to check whether a given spectrum is in the constraints of the constraints dict.
    # In retrospect, all of these could have been two functions, called check_conditions() and check_highlight().
    def analyze_conditions(self, conditions):
        '''This function turn the user constraints, which are read in as strings, into a conditions dictionary. The
        conditions dictionary is an attribute of the extraction object. This dictionary is used to filter what data is
        plotted when the user plots things. For example, the plotter will not plot observations that don't match the
        sol restriction in the constraints dictionary.'''
        for item in conditions:
            if item.startswith('Sols:'):
                self.conditions['Sols'] = self.analyze_sols(item[5:])
            if item.startswith('LTST:'):
                self.conditions['LTST'] = self.conditions['LTST']+item[5:]
            if item.startswith('Tau:'):
                self.conditions['Tau'] = self.conditions['Tau'] + item[4:].split('_')
            if item.startswith('Feature:'):
                if item.startswith('Feature:All Soils'):
                    self.conditions['Features'] = self.conditions['Features'] + ['Undisturbed Soil','Disturbed Soil']
                elif item.startswith('Feature:All Drill Fines'):
                    self.conditions['Features'] = self.conditions['Features'] + ['Drill Tailings','Dump Piles']
                elif item.startswith('Feature:All Rocks'): self.conditions['Features'] = self.conditions['Features'] + \
                                        ['Dust-Cleared Rock','Dusty Rock','Nodule-Rich Rock','Broken Rock Face']
                else: self.conditions['Features'] = self.conditions['Features']+[item[8:]]
            if item.startswith('Elevation:'):
                self.conditions['Elevation'] = self.conditions['Elevation']+[item[10:]]
            if item.startswith('Mcam:'):
                self.conditions['Mcam'] = self.conditions['Mcam'] + item[5:].split(',')
            if item.startswith('Float:'):
                self.conditions['Float'] = self.conditions['Float'] + [item[6:]]
            if item.startswith('Ls:'):
                self.conditions['Ls'] = self.conditions['Ls'] + [item[3:]]
            if item.startswith('Formation:'):
                self.conditions['Units'] = self.conditions['Units'] + [item[10:]]
            if item.startswith('Member:'):
                self.conditions['Units'] = self.conditions['Units'] + [item[7:]]
            if item.startswith('Group:'):
                self.conditions['Units'] = self.conditions['Units'] + [item[6:]]
            if item.startswith('Latitude:'):
                self.conditions['Latitude'] = self.conditions['Latitude'] + [item[9:]]
            if item.startswith('Longitude:'):
                self.conditions['Longitude'] = self.conditions['Longitude'] + [item[10:]]
            if item.startswith('Traverse:'):
                self.conditions['Traverse'] = self.conditions['Traverse'] + [item[9:]]

            if item.startswith('HighlightedSols:'):
                self.highlight['Sols'] = self.analyze_sols(item[16:])
            if item.startswith('HighlightedLTST:'):
                self.highlight['LTST'] = self.highlight['LTST']+item[16:]
            if item.startswith('HighlightedTau:'):
                self.highlight['Tau'] = self.highlight['Tau'] + item[15:].split('_')
            if item.startswith('HighlightedFeature:'):
                if item.startswith('HighlightedFeature:All Soils'):
                    self.highlight['Features'] = self.highlight['Features'] + ['Undisturbed Soil','Disturbed Soil']
                elif item.startswith('HighlightedFeature:All Drill Fines'):
                    self.highlight['Features'] = self.highlight['Features'] + ['Drill Tailings','Dump Piles']
                elif item.startswith('HighlightedFeature:All Rocks'): self.highlight['Features'] = self.highlight['Features'] + \
                                        ['Dust-Cleared Rock','Dusty Rock','Nodule-Rich Rock','Broken Rock Face']
                else:self.highlight['Features'] = self.highlight['Features']+[item[19:]]
            if item.startswith('HighlightedElevation:'):
                self.highlight['Elevation'] = self.highlight['Elevation']+[item[21:]]
            if item.startswith('HighlightedMcam:'):
                self.highlight['Mcam'] = self.highlight['Mcam'] + item[16:].split(',')
            if item.startswith('HighlightedFloat:'):
                self.highlight['Float'] = self.highlight['Float'] + [item[17:]]
            if item.startswith('HighlightedLs:'):
                self.highlight['Ls'] = self.highlight['Ls'] + [item[14:]]
            if item.startswith('HighlightedFormation:'):
                self.highlight['Units'] = self.highlight['Units'] + [item[21:]]
            if item.startswith('HighlightedMember:'):
                self.highlight['Units'] = self.highlight['Units'] + [item[18:]]
            if item.startswith('HighlightedGroup:'):
                self.highlight['Units'] = self.highlight['Units'] + [item[17:]]
            if item.startswith('HighlightedLatitude:'):
                self.highlight['Latitude'] = self.highlight['Latitude'] + [item[20:]]
            if item.startswith('HighlightedLongitude:'):
                self.highlight['Longitude'] = self.highlight['Longitude'] + [item[21:]]
            if item.startswith('HighlightedTraverse:'):
                self.highlight['Traverse'] = self.highlight['Traverse'] + [item[20:]]

    def analyze_sols(self, sol_string):
        '''This method turns a sol input string into a list, which is returned and added to the constraints dict.
        Example of viable input: 1,2,3_7,10_100,106,1010. Commas seperate individual days, underscores represent
        ranges.'''
        parts = sol_string.split(',')
        temp = []
        for sol in parts:
            if not '_' in sol:  # appends a single sol to the list
                temp.append(int(sol))
            if '_' in sol:  # appends a range of sols to the list
                first = int(sol.split('_')[0])
                last = int(sol.split('_')[1])
                for item in np.arange(first, last + 1, 1):
                    temp.append(item)
        return temp

    def in_LTST(self, LTST):

        if self.conditions['LTST'] == '': return True # If user did not select give an LTST limitation
        if LTST == '' or LTST == []: return False
        if isinstance(LTST, float):
            if pd.isnull(LTST): return False
        #LTST Example : 12:23:14_13:14:56
        hour1 = self.conditions['LTST'].split('_')[0].split(':')[0]
        min1 = self.conditions['LTST'].split('_')[0].split(':')[1]
        sec1 = self.conditions['LTST'].split('_')[0].split(':')[2]
        hour2 = self.conditions['LTST'].split('_')[1].split(':')[0]
        min2 = self.conditions['LTST'].split('_')[1].split(':')[1]
        sec2 = self.conditions['LTST'].split('_')[1].split(':')[2]
        hour = str(LTST).split(':')[0]
        min = str(LTST).split(':')[1]
        sec = str(LTST).split(':')[2]
        if hour < hour1 or hour > hour2: return False #The given LTST is out of hour range
        if hour == hour1 and min < min1: return False # hour is equal, but minute value is too little
        if hour == hour2 and min > min2: return False # Hour is equal, but minute value too high
        if hour == hour1 and min == min1 and sec < sec1 : return False # very unlikely, but you get the idea
        if hour == hour2 and min == min2 and sec > sec2 : return False #second value is out of range
        return True # if this point is reached, the LTST is in range!

    def in_tau(self,tau): #returns true or false
        if self.conditions['Tau'] == [] : return True # If user did not select Tau limitation
        elif tau >= float(self.conditions['Tau'][0]) and tau <= float(self.conditions['Tau'][1]): return True
        else: return False

    def in_features(self,feature):
        if self.conditions['Features'] == [] : return True # If user did not select Feature limitation
        elif feature in self.conditions['Features'] : return True
        else : return False

    def in_elevation(self,elevation):
        if self.conditions['Elevation'] == []: return True  # If user did not select Feature limitation
        low = float(self.conditions['Elevation'][0].split('_')[0])
        high = float(self.conditions['Elevation'][0].split('_')[1])
        if elevation >= low and elevation <= high: return True
        return False

    def in_sol_range(self, sol):
        if self.conditions['Sols'] == [] : return True # If user does not input sol limitation
        elif int(sol) in self.conditions['Sols'] : return True
        else : return False

    def in_mcam(self, seq):
        if self.conditions['Mcam'] == []: return True
        if seq in self.conditions['Mcam']: return True
        else: return False

    def in_unit(self, unit):
        if self.conditions['Units'] == []: return True
        if unit in [unit.lower() for unit in self.conditions['Units']]: return True
        else: return False

    def in_highlighted_LTST(self, LTST): #returns true or false
        if self.highlight['LTST'] == '': return False
        if LTST == '' or LTST == []: return False
        if isinstance(LTST, float):
            if math.isnan(LTST): return False
        #LTST Example : 12:23:14_13:14:56
        hour1 = self.highlight['LTST'].split('_')[0].split(':')[0]
        min1 = self.highlight['LTST'].split('_')[0].split(':')[1]
        sec1 = self.highlight['LTST'].split('_')[0].split(':')[2]
        hour2 = self.highlight['LTST'].split('_')[1].split(':')[0]
        min2 = self.highlight['LTST'].split('_')[1].split(':')[1]
        sec2 = self.highlight['LTST'].split('_')[1].split(':')[2]
        hour = str(LTST).split(':')[0]
        min = str(LTST).split(':')[1]
        sec = str(LTST).split(':')[2]
        if hour < hour1 or hour > hour2: return False #The given LTST is out of hour range
        if hour == hour1 and min < min1: return False # hour is equal, but minute value is too little
        if hour == hour2 and min > min2: return False # Hour is equal, but minute value too high
        if hour == hour1 and min == min1 and sec < sec1 : return False # very unlikely, but you get the idea
        if hour == hour2 and min == min2 and sec > sec2 : return False #second value is out of range
        return True # if this point is reached, the LTST is in range!

    def in_highlighted_tau(self,tau): #returns true or false
        if self.highlight['Tau'] == []: return False
        if tau >= float(self.highlight['Tau'][0]) and tau <= float(self.highlight['Tau'][1]): return True
        return False

    def in_highlighted_features(self,feature):
        if self.highlight['Features'] == []: return False
        if feature in self.highlight['Features'] : return True
        return False

    def in_highlighted_elevation(self,elevation):
        if self.highlight['Elevation'] == []: return False
        low = float(self.highlight['Elevation'][0].split('_')[0])
        high = float(self.highlight['Elevation'][0].split('_')[1])
        if elevation >= low and elevation <= high: return True
        return False

    def in_highlighted_sol_range(self, sol):
        if self.highlight['Sols'] == []: return False
        if int(sol) in self.highlight['Sols'] : return True
        return False

    def in_highlighted_mcam(self, seq):
        while seq[0] == '0':
            seq = seq[1:]
        if self.highlight['Mcam'] == []: return False
        if str(seq) in self.highlight['Mcam']: return True
        return False

    def in_highlighted_unit(self, unit):
        if self.highlight['Units'] == []: return False
        if unit in self.highlight['Units']: return True
        return False

    def in_highlighted_tina(self, unit):
        if self.highlight['Tina'] == []: return False
        if unit in self.highlight['Tina']: return True
        return False

    def in_float(self, designation):
        if self.conditions['Float'] == []: return True
        if designation in self.conditions['Float']: return True
        else: return False

    def in_highlighted_float(self, designation):
        if self.highlight['Float'] == []: return False
        if designation in self.highlight['Float']: return True
        return False

    def in_ls(self, ls):
        if self.conditions['Ls'] == [] : return True # If user does not input Ls limitation
        # ['1_30,40_50','70_80'] Example input

        if self.conditions['Ls'] != []:
            for item in self.conditions['Ls']:
                for pair in item.split(','):
                    lower = float(pair.split('_')[0])
                    upper = float(pair.split('_')[1])
                    if ls <= upper and ls >= lower:
                        return True
        else : return False

    def in_highlighted_ls(self, ls):
        if self.highlight['Ls'] == [] : return False # If user does not input Ls limitation
        # [1_30, 40_50, 70_80] Example input
        if self.highlight['Ls'] != []:
            for item in self.highlight['Ls']:
                for pair in item.split(','):
                    lower = float(pair.split('_')[0])
                    upper = float(pair.split('_')[1])
                    if ls <= upper and ls >= lower:
                        return True
        else : return False

    def in_latitude(self, lat):
        if self.conditions['Latitude'] == [] : return True # If user does not input Ls limitation
        # ['1_30,40_50','70_80'] Example input

        if self.conditions['Latitude'] != []:
            for item in self.conditions['Latitude']:
                for pair in item.split(','):
                    lower = float(pair.split('_')[0])
                    upper = float(pair.split('_')[1])
                    if lower > upper:
                        lower, upper = upper, lower
                    if lat <= upper and lat >= lower:
                        return True
        return False

    def in_highlighted_latitude(self, lat):
        if self.highlight['Latitude'] == [] : return False # If user does not input Ls limitation
        # ['1_30,40_50','70_80'] Example input

        if self.highlight['Latitude'] != []:
            for item in self.highlight['Latitude']:
                for pair in item.split(','):
                    lower = float(pair.split('_')[0])
                    upper = float(pair.split('_')[1])
                    if lat <= upper and lat >= lower:
                        return True
        else : return False

    def in_longitude(self, long):
        if self.conditions['Longitude'] == []: return True  # If user does not input Ls limitation
        # ['1_30,40_50','70_80'] Example input

        if self.conditions['Longitude'] != []:
            for item in self.conditions['Longitude']:
                for pair in item.split(','):
                    lower = float(pair.split('_')[0])
                    upper = float(pair.split('_')[1])
                    if long <= upper and long >= lower:
                        return True
        return False

    def in_highlighted_longitude(self, long):
        if self.highlight['Longitude'] == []: return False  # If user does not input Ls limitation
        # ['1_30,40_50','70_80'] Example input

        if self.highlight['Longitude'] != []:
            for item in self.highlight['Longitude']:
                for pair in item.split(','):
                    lower = float(pair.split('_')[0])
                    upper = float(pair.split('_')[1])
                    if long <= upper and long >= lower:
                        return True
        else:
            return False

    def in_traverse(self, trav):
        if self.conditions['Traverse'] == []: return True  # If user does not input Ls limitation
        # ['1_30,40_50','70_80'] Example input

        if self.conditions['Traverse'] != []:
            for item in self.conditions['Traverse']:
                for pair in item.split(','):
                    lower = float(pair.split('_')[0])
                    upper = float(pair.split('_')[1])
                    if trav <= upper and trav>= lower:
                        return True
        return False

    def in_highlighted_traverse(self, trav):
        if self.highlight['Traverse'] == []: return False  # If user does not input Ls limitation
        # ['1_30,40_50','70_80'] Example input

        if self.highlight['Traverse'] != []:
            for item in self.highlight['Traverse']:
                for pair in item.split(','):
                    lower = float(pair.split('_')[0])
                    upper = float(pair.split('_')[1])
                    if trav <= upper and trav >= lower:
                        return True
        else:
            return False


class Application(tk.Frame):
    '''This class is a bit of a mess, but is responsible for the GUI that the user interacts with. There's a lot of
    methods and things that make this class function.'''

    #These dictionaries are used by a few functions, they correlate color between strings, features, and RGB values.
    color2feature = {  # Converts color to associated feature
        'red': 'Undisturbed Soil',
        'light green': 'Dust-Cleared Rock',
        'light blue': 'Dusty Rock',
        'dark blue': 'Dusty Rock',
        'teal': 'Dusty Rock',
        'light cyan': 'Nodule-rich Rock',
        'dark green': 'Broken Rock Face',
        'light purple': 'Drill Tailings',
        'pink': 'Other',
        'sienna': 'Other',
        'bright red': 'Undisturbed Soil',
        'yellow': 'Other',
        'goldenrod': 'Veins',
        'dark red': 'Disturbed Soil',
        'dark purple': 'Dump Piles'}
    feature2color = {
        'Undisturbed Soil': 'red',
        'Dust-Cleared Rock': 'light green',
        'Dusty Rock': 'light blue',
        'Nodule-rich Rock': 'light cyan',
        'Broken Rock Face': 'dark green',
        'Drill Tailings': 'light purple',
        'Other': 'yellow',
        'Veins': 'goldenrod',
        'Disturbed Soil': 'dark red',
        'Dump Piles': 'dark purple'
    }
    color2RGB = {  # Converts color to associated feature
        'red': (212 / 255, 0, 42 / 255),
        'light green': (109 / 255, 225 / 255, 0),
        'light blue': (0, 0, 225 / 255),
        'dark blue': (0, 0, 111 / 255),
        'teal': (0, 110 / 255, 109 / 255),
        'light cyan': (0, 255 / 255, 255 / 255),
        'dark green': (16 / 255, 113 / 255, 1 / 255),
        'light purple': (254 / 255, 0, 255 / 255),
        'pink': (249 / 255, 104 / 255, 92 / 255),
        'sienna': (143 / 255, 63 / 255, 30 / 255),
        'bright red': (255 / 255, 0, 0),
        'yellow': (255 / 255, 255 / 255, 0),
        'goldenrod': (170 / 255, 117 / 255, 0),
        'dark red': (109 / 255, 0, 0),
        'dark purple': (108 / 255, 0, 110 / 255)}
    eye2letter = {
        ' LEFT': ' L',
        ' RIGHT': ' R'
    }

    def __init__(self, master=None):
        super().__init__(master)
        self.pack()
        self.winfo_toplevel().title('WWU Mastcam Spectral Plotter')
        root.wm_iconbitmap('MSP.ico')
        self.create_widgets()
        self.create_menu()
        self.spectra_list = pickle.load(open('spectra_list.p','rb')) # Loads in spectra
        # self.user_spectra_list = pickle.load(open('user_spectra.p','rb'))
        root.bind('<space>', lambda x:self.plot(Extraction(self))) # Can plot with space bar
        self.home_path = os.getcwd()

    def create_menu(self):
        '''This function creates the menu that holds various plotter options.'''
        self.menubar = tk.Menu(self)
        menu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="File", menu=menu)  # File Menu
        menu.add_command(label='Import from Server', command = lambda: self.import_spectra())
        menu.add_separator()
        self.our_spectra = tk.IntVar(value=1)
        self.user_spectra = tk.IntVar(value=1)
        menu.add_checkbutton(label='Include WWU Spectra', variable=self.our_spectra)
        menu.add_checkbutton(label='Include User Spectra', variable=self.user_spectra)
        menu.add_separator()
        menu.add_command(label="Plot", command= lambda: self.plot(Extraction(self)))
        menu.add_command(label='Export CSV', command = lambda: self.export_csv(Extraction(self)))
        menu.add_separator()
        menu.add_command(label='Save Template', command = lambda: self.save_template(Extraction(self)))
        menu.add_command(label='Load Template', command=lambda: self.load_template())
        menu.add_separator()
        menu.add_command(label='Quit', command=root.destroy)

        menu = tk.Menu(self.menubar, tearoff=0)
        submenu = tk.Menu(self.menubar, tearoff=0)
        self.norm_point = tk.IntVar(value = 6)
        submenu.add_radiobutton(label='None', variable=self.norm_point, value=0)
        submenu.add_radiobutton(label = '446',variable = self.norm_point, value = 1)
        submenu.add_radiobutton(label='494',variable = self.norm_point, value = 2)
        submenu.add_radiobutton(label='527',variable = self.norm_point, value = 3)
        submenu.add_radiobutton(label='552.5',variable = self.norm_point, value = 4)
        submenu.add_radiobutton(label='639',variable = self.norm_point, value = 5)
        submenu.add_radiobutton(label='1012.5',variable = self.norm_point, value = 6)
        self.menubar.add_cascade(label="Options", menu=menu)  # Options Menu
        self.r_star = tk.IntVar(value=1)
        menu.add_checkbutton(label='R* correction',variable = self.r_star)
        self.geologic_members = tk.IntVar(value = 1)
        menu.add_checkbutton(label='Display geologic boundaries in elevation plot', variable = self.geologic_members)
        menu.add_separator()
        menu.add_cascade(label='Normalization Point', menu = submenu)
        menu.add_separator()
        self.propogate_error = tk.IntVar(value = 0)
        menu.add_radiobutton(label='No Error Bars', variable = self.propogate_error, value = 0)
        menu.add_radiobutton(label='All Error Bars', variable=self.propogate_error, value=1)
        menu.add_radiobutton(label='Representative Error Bars', variable=self.propogate_error, value=2)
        menu.add_separator()
        colormenu = tk.Menu(self.menubar, tearoff=0)
        self.color_type = tk.IntVar(value = 3)
        colormenu.add_radiobutton(label='Black', variable=self.color_type, value=0)
        colormenu.add_radiobutton(label='Simplified Feature Color Scheme', variable=self.color_type, value=3)
        colormenu.add_radiobutton(label='Standard Feature Color Scheme', variable=self.color_type, value=1)
        colormenu.add_radiobutton(label='Color from Target Elevation', variable=self.color_type, value=2)
        colormenu.add_radiobutton(label='Color from Sol', variable=self.color_type, value=4)
        colormenu.add_radiobutton(label='Color from Tau', variable=self.color_type, value=5)
        colormenu.add_radiobutton(label='Color from Ls', variable=self.color_type, value=6)
        menu.add_cascade(label='Color Schemes', menu=colormenu)
        self.master.config(menu=self.menubar)

    def create_widgets(self):
        '''This function creates all widgets. The widgets are things the user can click and select, like the filters,
        constraints, the scrollbars, and X and Y axis options.'''
        # Basic Labels
        tk.Label(self, text="X axis:").grid(row=1, column=0)
        tk.Label(self, text=" ").grid(row=2, column=0)
        tk.Label(self, text="Y axis:").grid(row=3, column=0)
        tk.Label(self, text=" ").grid(row=4, column=0)
        tk.Label(self, text = "Conditions").grid(row = 5, column = 0)
        tk.Label(self, text="Selection").grid(row=0, column=1)

        # Dynamic X axis labels
        self.x_label1 = tk.Label(self, text="1")
        self.x_label1.grid(row=0, column=2)
        self.x_label2 = tk.Label(self, text="2")
        self.x_label2.grid(row=0, column=3)
        self.x_label3 = tk.Label(self, text="3")
        self.x_label3.grid(row=0, column=4)

        # Initialize X listbox widgets
        self.x_expand1 = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        self.x_expand1.grid(row=1, column=2)
        self.x_expand2 = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        self.x_expand2.grid(row=1, column=3)
        self.x_expand3 = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        self.x_expand3.grid(row=1, column=4, sticky = 'W')
        self.populate_list(self.x_expand1)
        self.populate_list(self.x_expand2)
        self.populate_list(self.x_expand3)
        self.x_pca = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        for n in np.arange(1, 10, 1):
            self.x_pca.insert(tk.END, 'PC '+str(n))

        # X widget scrollbar
        self.scrollbar = tk.Scrollbar(self)
        self.x_expand1.config(yscrollcommand=self.scrollbar.set)
        self.x_expand2.config(yscrollcommand=self.scrollbar.set)
        self.x_expand3.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=self.yview)
        self.scrollbar.grid(row=1, column=6, sticky=tk.NS)

        # Initialize Y axis widgets
        self.y_label1 = tk.Label(self, text="1")
        self.y_label1.grid(row=2, column=2)
        self.y_label2 = tk.Label(self, text="2")
        self.y_label2.grid(row=2, column=3)
        self.y_label3 = tk.Label(self, text="3")
        self.y_label3.grid(row=2, column=4)
        self.y_expand1 = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        self.y_expand1.grid(row=3, column=2)
        self.y_expand2 = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        self.y_expand2.grid(row=3, column=3)
        self.y_expand3 = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        self.y_expand3.grid(row=3, column=4, sticky = 'W')
        self.populate_list(self.y_expand1)
        self.populate_list(self.y_expand2)
        self.populate_list(self.y_expand3)
        self.y_pca = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        for n in np.arange(1, 10, 1):
            self.y_pca.insert(tk.END, 'PC '+str(n))

        # Y widget scrollbar
        self.scrollbar2 = tk.Scrollbar(self)
        self.y_expand1.config(yscrollcommand=self.scrollbar2.set)
        self.y_expand2.config(yscrollcommand=self.scrollbar2.set)
        self.y_expand3.config(yscrollcommand=self.scrollbar2.set)
        self.scrollbar2.config(command=self.yview2)
        self.scrollbar2.grid(row=3, column=6, sticky=tk.NS)

        # Populate X axis selection listbox
        self.x_options = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        for item in ['Reflectance', 'Band Depth', 'Slope', 'Ratio', 'Sol', 'Rover Elevation', 'Target Elevation',
                     'LTST','Tau','PCA', 'Ls', 'Latitude', 'Longitude', 'Traverse Distance']:
            self.x_options.insert(tk.END, item)
        self.x_options.grid(row=1, column=1)
        self.x_options.bind('<<ListboxSelect>>', self.x_show_options)

        # Y axis selection listbox
        self.y_options = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        for item in ['Reflectance', 'Band Depth', 'Slope', 'Ratio', 'Sol', 'Rover Elevation', 'Target Elevation',
                     'LTST','Tau','PCA', 'Ls', 'Latitude', 'Longitude', 'Traverse Distance']:
            self.y_options.insert(tk.END, item)
        self.y_options.grid(row=3, column=1)
        self.y_options.bind('<<ListboxSelect>>', self.y_show_options)

        # Condition selection
        self.conditions = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        self.conditions.grid(row=5, column=1, rowspan=5)
        for item in ['Feature Types','Sol Range','LTST Range','Tau Range','Elevation Range','Mcam IDs',
                     'Group','Member','Formation','Float','Ls','Latitude','Longitude','Traverse Distance']:
            self.conditions.insert(tk.END, item)
        self.conditions.bind('<<ListboxSelect>>', self.condition_options)

        #Condition Option Widgets
        self.feature_types_input = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        for item in ['All Soils','All Rocks','All Drill Fines','Undisturbed Soil','Disturbed Soil','Dust-Cleared Rock',
                     'Dusty Rock','Nodule-Rich Rock',
                     'Broken Rock Face','Drill Tailings','Veins','Dump Piles','Other']:
            self.feature_types_input.insert(tk.END, item)

        self.member_input = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        for item in ['Stimson','Murray']:
            self.member_input.insert(tk.END, item)

        self.group_input = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        for item in ['Mount Sharp','Bradbury']:
            self.group_input.insert(tk.END, item)

        self.formation_input = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        for item in ['Jura','Pettegrove Point','Blunts Point','Sutton Island','Karasburg','Hartmann\'s Valley',
                     'Pahrump Hills','Kimberly','Yellowknife Bay']:
            self.formation_input.insert(tk.END, item)

        self.tina_unit_types_input = tk.Listbox(self, selectmode='SINGLE', exportselection=False)
        for item in ['Murray','Stimson']:
            self.tina_unit_types_input.insert(tk.END, item)

        self.sol_range_input = tk.Entry(self)
        self.ls_input = tk.Entry(self)
        self.LTST_input = tk.Entry(self)
        self.elevation_input = tk.Entry(self)
        self.tau_input = tk.Entry(self)
        self.mcam_input = tk.Entry(self)
        self.latitude_input = tk.Entry(self)
        self.longitude_input = tk.Entry(self)
        self.traverse_input = tk.Entry(self)

        self.float_input = tk.Listbox(self, selectmode = 'SINGLE', exportselection = False)
        for item in ['Float', 'No Float']: self.float_input.insert(tk.END, item)

        self.condition_list = tk.Listbox(self, selectmode = 'Single', exportselection = True, width = 50)
        self.condition_list.insert(tk.END, 'LTST:10:30:00_13:30:00') # Adds LTST limitation by default
        self.condition_list.grid(row=5, column = 4, rowspan=3,columnspan=2)

        self.condition_input_mode = "Feature Types"

        self.apply = tk.Button(self, text='Apply ->', command= lambda: self.add_condition())
        self.highlight = tk.Button(self, text = 'Highlight ->', command= lambda: self.add_highlight())
        self.remove = tk.Button(self, text='<- Remove', command=lambda: self.remove_condition())
        self.apply.grid(row = 5, column=3)
        self.highlight.grid(row=6, column=3)
        self.remove.grid(row=7,column=3)

    def remove_condition(self):
        '''Removes a condition from the condition list'''
        self.condition_list.delete(self.condition_list.curselection())

    def add_condition(self):
        '''Controls the add condition button'''
        output = self.condition_list
        mode = self.condition_input_mode
        if mode == "Feature Types":
            output.insert(tk.END, 'Feature:' + self.feature_types_input.get(self.feature_types_input.curselection()))
        if mode == 'Sol Range':
            if self.sol_range_input.get() != '':
                output.insert(tk.END, 'Sols:' + self.sol_range_input.get())
        if mode == 'LTST Range':
            if self.LTST_input.get() != '':
                output.insert(tk.END, 'LTST:' + self.LTST_input.get())
        if mode == 'Elevation Range':
            if self.elevation_input.get() != '':
                output.insert(tk.END, 'Elevation:' + self.elevation_input.get())
        if mode == 'Tau Range':
            if self.tau_input.get() != '':
                output.insert(tk.END, 'Tau:' + self.tau_input.get())
        if mode == 'Mcam IDs':
            if self.mcam_input.get() != '':
                output.insert(tk.END, 'Mcam:' + self.mcam_input.get())
        if mode == 'Member':
            output.insert(tk.END, 'Member:' + self.member_input.get(self.member_input.curselection()))
        if mode == 'Group':
            output.insert(tk.END, 'Group:' + self.group_input.get(self.group_input.curselection()))
        if mode == 'Formation':
            output.insert(tk.END, 'Formation:' + self.formation_input.get(self.formation_input.curselection()))
        if mode == 'Float':
            output.insert(tk.END,'Float:'+self.float_input.get(self.float_input.curselection()))
        if mode == 'Ls':
            if self.ls_input.get() != '':
                output.insert(tk.END, 'Ls:' + self.ls_input.get())
        if mode == 'Latitude':
            if self.latitude_input.get() != '':
                output.insert(tk.END, 'Latitude:' + self.latitude_input.get())
        if mode == 'Longitude':
            if self.longitude_input.get() != '':
                output.insert(tk.END, 'Longitude:' + self.longitude_input.get())
        if mode == 'Traverse Distance':
            if self.traverse_input.get() != '':
                output.insert(tk.END, 'Traverse:' + self.traverse_input.get())

    def add_highlight(self):
        '''Controls the add highlight button'''
        output = self.condition_list
        mode = self.condition_input_mode
        if mode == "Feature Types":
            output.insert(tk.END, 'HighlightedFeature:' + self.feature_types_input.get(self.feature_types_input.curselection()))
        if mode == 'Sol Range':
            if self.sol_range_input.get() != '':
                output.insert(tk.END, 'HighlightedSols:' + self.sol_range_input.get())
        if mode == 'LTST Range':
            if self.LTST_input.get() != '':
                output.insert(tk.END, 'HighlightedLTST:' + self.LTST_input.get())
        if mode == 'Elevation Range':
            if self.elevation_input.get() != '':
                output.insert(tk.END, 'HighlightedElevation:' + self.elevation_input.get())
        if mode == 'Tau Range':
            if self.tau_input.get() != '':
                output.insert(tk.END, 'HighlightedTau:' + self.tau_input.get())
        if mode == 'Mcam IDs':
            if self.mcam_input.get() != '':
                output.insert(tk.END, 'HighlightedMcam:' + self.mcam_input.get())
        if mode == 'Float':
            output.insert(tk.END, 'HighlightedFloat:' + self.float_input.get(self.float_input.curselection()))
        if mode == 'Tina':
            output.insert(tk.END,'HighlightedTina:' + self.tina_unit_types_input.get(self.tina_unit_types_input.curselection()).lower())
        if mode == 'Ls':
            if self.ls_input.get() != '':
                output.insert(tk.END, 'HighlightedLs:' + self.ls_input.get())
        if mode == 'Latitude':
            if self.latitude_input.get() != '':
                output.insert(tk.END, 'HighlightedLatitude:' + self.latitude_input.get())
        if mode == 'Longitude':
            if self.longitude_input.get() != '':
                output.insert(tk.END, 'HighlightedLongitude:' + self.longitude_input.get())
        if mode == 'Traverse Distance':
            if self.traverse_input.get() != '':
                output.insert(tk.END, 'HighlightedTraverse:' + self.traverse_input.get())
        if mode == 'Member':
            output.insert(tk.END, 'HighlightedMember:' + self.member_input.get(self.member_input.curselection()))
        if mode == 'Group':
            output.insert(tk.END, 'HighlightedGroup:' + self.group_input.get(self.group_input.curselection()))
        if mode == 'Formation':
            output.insert(tk.END, 'HighlightedFormation:' + self.formation_input.get(self.formation_input.curselection()))

    def condition_options(self, evt):
        '''When a constraint is clicked, this function updates the widget to the right to let the user input info
        on the constraint they are inteested in'''
        w = evt.widget
        selection = w.curselection()[0]
        self.feature_types_input.grid_forget() # Clear entries before packing new widget
        self.sol_range_input.grid_forget()
        self.LTST_input.grid_forget()
        self.elevation_input.grid_forget()
        self.tau_input.grid_forget()
        self.mcam_input.grid_forget()
        self.member_input.grid_forget()
        self.formation_input.grid_forget()
        self.group_input.grid_forget()
        self.float_input.grid_forget()
        self.tina_unit_types_input.grid_forget()
        self.ls_input.grid_forget()
        self.latitude_input.grid_forget()
        self.longitude_input.grid_forget()
        self.traverse_input.grid_forget()
        if selection == 0: #Feature Types
            self.feature_types_input.grid(row=5, column=2, rowspan=3)
            self.condition_input_mode = "Feature Types"
        if selection == 1: # Sol Range
            self.sol_range_input.grid(row = 5, column = 2, rowspan=3)
            self.condition_input_mode = 'Sol Range'
        if selection == 2: # LTST Range
            self.LTST_input.grid(row = 5, column = 2, rowspan=3)
            self.condition_input_mode = 'LTST Range'
        if selection == 3: # Tau Tange
            self.tau_input.grid(row = 5, column = 2, rowspan=3)
            self.condition_input_mode = 'Tau Range'
        if selection == 4: # Elevation Range
            self.elevation_input.grid(row = 5, column = 2, rowspan=3)
            self.condition_input_mode = 'Elevation Range'
        if selection == 5: #Mcam IDs
            self.mcam_input.grid(row = 5, column = 2, rowspan=3)
            self.condition_input_mode = 'Mcam IDs'
        if selection == 6: # Group
            self.group_input.grid(row=5, column=2, rowspan=3)
            self.condition_input_mode = 'Group'
        if selection == 7: #Formation
            self.formation_input.grid(row=5, column=2, rowspan=3)
            self.condition_input_mode = 'Formation'
        if selection == 8: #Member
            self.member_input.grid(row=5, column=2, rowspan=3)
            self.condition_input_mode = 'Member'
        if selection == 9: #Float
            self.float_input.grid(row=5, column=2, rowspan=3)
            self.condition_input_mode = 'Float'
        if selection == 10: # Ls Range
            self.ls_input.grid(row = 5, column = 2, rowspan=3)
            self.condition_input_mode = 'Ls'
        if selection == 11: # Latitude Range
            self.latitude_input.grid(row = 5, column = 2, rowspan=3)
            self.condition_input_mode = 'Latitude'
        if selection == 12: # Longitude Range
            self.longitude_input.grid(row = 5, column = 2, rowspan=3)
            self.condition_input_mode = 'Longitude'
        if selection == 13: # Traverse Range
            self.traverse_input.grid(row = 5, column = 2, rowspan=3)
            self.condition_input_mode = 'Traverse Distance'

    def save_template(self, extract):
        '''Used to save templates (obviously). Created an Extraction object and saves it as a pickle file.'''
        answer = filedialog.asksaveasfilename(parent=self,
                                              initialdir=os.getcwd()+'/Templates',
                                              title="Please select a file name for saving:",
                                              filetypes = [('Pickle File','.p')],
                                              defaultextension='.p')
        pickle.dump(extract, open(answer, 'wb')) # saves an extract to be read in later

    def load_template(self):
        '''Used to load in templates, which are Extraction objects.'''
        #Need to load in R* and all that
        chosen = filedialog.askopenfilename(parent=self,
                                            initialdir=os.getcwd()+'/Templates',
                                            title="Please select a file:",
                                            filetypes=[('Pickle File', '.p')])
        def select(widget, item): # method to select a given item in a list
            contents = widget.get(0, tk.END)
            index = contents.index(item)
            widget.select_clear(0,tk.END)
            widget.activate(index)
            widget.selection_set(index)
        template = pickle.load(open(chosen,'rb')) # loads in an extract, copies over the fields
        select(self.x_options, template.x_choice)
        select(self.x_expand1, template.x_1)
        select(self.x_expand2, template.x_2)
        select(self.x_expand3, template.x_3)
        self.x_options.event_generate('<<ListboxSelect>>') # simulates a click on listbox to update fields
        select(self.y_options, template.y_choice)
        select(self.y_expand1, template.y_1)
        select(self.y_expand2, template.y_2)
        select(self.y_expand3, template.y_3)
        self.y_options.event_generate('<<ListboxSelect>>')
        self.condition_list.select_clear(0,tk.END)
        for item in template.conditions_list:
            if not item in self.condition_list.get(0,tk.END):
                self.condition_list.insert(tk.END, item)

    def yview(self, *args):  # Function that lets scrollbar scroll all boxes
        '''Used to scroll through widget options'''
        self.x_expand1.yview(*args)
        self.x_expand2.yview(*args)
        self.x_expand3.yview(*args)

    def yview2(self, *args):  # Function that lets scrollbar scroll all boxes
        '''Used to scroll through widget options'''
        self.y_expand1.yview(*args)
        self.y_expand2.yview(*args)
        self.y_expand3.yview(*args)

    def x_show_options(self, evt):
        '''When a new X axis opion is clicked, this updates the X selection widgets accordingly.'''
        w = evt.widget
        selection = w.curselection()[0]
        self.x_pca.grid_remove()
        if selection == 0:  # Reflection
            self.x_expand1.grid(row=1, column=2)
            self.x_label1.grid(row = 0, column = 2)
            self.x_expand2.grid_remove()
            self.x_expand3.grid_remove()
            self.x_label1.config(text='Filter (nm)')
            self.x_label2.grid_remove()
            self.x_label3.grid_remove()
        elif selection == 1:  # Band depth
            self.x_expand1.grid(row=1, column=2)
            self.x_label1.grid(row=0, column=2)
            self.x_expand2.grid(row=1, column=3)
            self.x_expand3.grid(row=1, column=4)
            self.x_label1.config(text='Left Shoulder (nm)')
            self.x_label2.config(text='Middle Filter (nm)')
            self.x_label2.grid(row=0, column=3)
            self.x_label3.config(text='Right Shoulder (nm)')
            self.x_label3.grid(row=0, column=4)
        elif selection == 2:  # Slope
            self.x_expand1.grid(row=1, column=2)
            self.x_expand2.grid(row=1, column=3)
            self.x_expand3.grid_remove()
            self.x_label1.config(text='First Filter (nm)')
            self.x_label1.grid(row=0, column=2)
            self.x_label2.config(text='Second Filter (nm)')
            self.x_label2.grid(row=0, column=3)
            self.x_label3.grid_remove()
        elif selection == 3:  # Ratio
            self.x_expand1.grid(row=1, column=2)
            self.x_label1.grid(row=0, column=2)
            self.x_expand2.grid(row=1, column=3)
            self.x_expand3.grid_remove()
            self.x_label1.config(text='Numerator Filter (nm)')
            self.x_label2.config(text='Denominator Filter (nm)')
            self.x_label2.grid(row=0, column=3)
        elif selection == 9: # PCA
            self.x_expand1.grid_remove()
            self.x_expand2.grid_remove()
            self.x_expand3.grid_remove()
            self.x_label1.grid_remove()
            self.x_label2.grid_remove()
            self.x_label3.grid_remove()
            self.x_pca.grid(row=1, column=2)
        else:  # Everything else
            self.x_label3.grid_remove()
            self.x_expand1.grid_remove()
            self.x_expand2.grid_remove()
            self.x_expand3.grid_remove()
            self.x_label1.grid_remove()
            self.x_label2.grid_remove()
            self.x_label3.grid_remove()

    def y_show_options(self, evt):
        '''When a new Y axis opion is clicked, this updates the Y selection widgets accordingly.'''
        try: w = evt.widget
        except: w = evt
        selection = w.curselection()[0]
        self.y_pca.grid_remove()
        if selection == 0:  # Reflection
            self.y_expand1.grid(row=3, column=2)
            self.y_expand2.grid_remove()
            self.y_expand3.grid_remove()
            self.y_label1.config(text='Filter (nm)')
            self.y_label1.grid(row=2, column=2)
            self.y_label2.grid_remove()
            self.y_label3.grid_remove()
        elif selection == 1:  # Band depth
            self.y_expand1.grid(row=3, column=2)
            self.y_expand2.grid(row=3, column=3)
            self.y_expand3.grid(row=3, column=4)
            self.y_label1.config(text='Left Shoulder (nm)')
            self.y_label1.grid(row=2, column=2)
            self.y_label2.config(text='Middle Filter (nm)')
            self.y_label2.grid(row=2, column=3)
            self.y_label3.config(text='Right Shoulder (nm)')
            self.y_label3.grid(row=2, column=4)
        elif selection == 2:  # Slope
            self.y_expand1.grid(row=3, column=2)
            self.y_expand2.grid(row=3, column=3)
            self.y_expand3.grid_remove()
            self.y_label1.config(text='First Filter (nm)')
            self.y_label1.grid(row=2, column=2)
            self.y_label2.config(text='Second Filter (nm)')
            self.y_label2.grid(row=2, column=3)
            self.y_label3.grid_remove()
        elif selection == 3:  # Ratio
            self.y_expand1.grid(row=3, column=2)
            self.y_expand2.grid(row=3, column=3)
            self.y_expand3.grid_remove()
            self.y_label1.config(text='Numerator Filter (nm)')
            self.y_label1.grid(row=2, column=2)
            self.y_label2.config(text='Denominator Filter (nm)')
            self.y_label2.grid(row=2, column=3)
            self.y_label3.grid_remove()
        elif selection == 9: # PCA
            self.y_expand1.grid_remove()
            self.y_expand2.grid_remove()
            self.y_expand3.grid_remove()
            self.y_label1.grid_remove()
            self.y_label2.grid_remove()
            self.y_label3.grid_remove()
            self.y_pca.grid(row=3, column=2)
        else:
            self.y_expand1.grid_remove()
            self.y_expand2.grid_remove()
            self.y_expand3.grid_remove()
            self.y_label1.grid_remove()
            self.y_label2.grid_remove()
            self.y_label3.grid_remove()

    def populate_list(self, listbox):
        '''This function is used to add all the wavelength selection options to the corresponding widgets'''
        filter_list = ['445 L', '446 A', '447 R', '493 R', '494 A', '495 L', '527 L', '527 A', '527 R',
                       '551 R', '552.5','554 L', '638 R', '639 A', '640 L',
                       '676 L', '751 L', '805 R', '867 L', '908 R', '937 R', '1012 L', '1012.5 A', '1013 R']
        listbox.delete(0, tk.END)
        for item in filter_list:
            listbox.insert(tk.END, item)

    def import_metadata(self):
        '''Loading is the metadata from the server. This is called when the import data function is called'''
        if 'win32' in sys.platform:
            pathstart = '//ion.physics.wwu.edu'
        else:
            pathstart = '//Volumes'
        reference = pd.read_csv(pathstart+'/ricedata/MarsGroup/Mastcam/Metadata/Metadata.csv').values
        pickle.dump(reference, open('metadata.p', 'wb'))  # This saves the data client-side

    def get_spectra_from_directory(self, path):
        '''The main function used to load in spectra from csvs. THis one is worth looking through in detail.'''
        reference = pickle.load(open(self.home_path + '/metadata.p', 'rb'))
        ref_mcam = reference[:, 1]
        ref_mcam = np.array([int(r[4:]) for r in ref_mcam])  # Cuts out the mcam
        ref_name = reference[:, 2]
        ref_rover_elevation = reference[:, 3]
        ref_mastcam_elevation = reference[:, 4]
        ref_tau = reference[:, 5]
        ref_LTST = reference[:, 6]
        ref_focal_distance = reference[:, 7]
        ref_incident = reference[:, 8]
        ref_emission = reference[:, 9]
        ref_ls = reference[:, 10]
        ref_site = reference[:, 11]
        ref_drive = reference[:, 12]
        ref_lat = reference[:, 13]
        ref_long = reference[:, 14]
        ref_traverse = reference[:, 15]

        def get_index(mcam):
            mcam = int(mcam)
            index, = np.where(ref_mcam == mcam)
            return index[0]

        def calculate_r_star(reflectance, incident):
            if reflectance == '' or reflectance == ' ': return ''
            if incident == '' or incident == ' ': return ''
            return float(reflectance) / math.cos(math.radians(float(incident)))

        def make_float(item):
            if isinstance(item, float):
                return item
            if isinstance(item, str):
                item = item.strip(' ')
                item = item.lower()
                if item.startswith('inf'):
                    return float('inf')
                if item.startswith('-inf'):
                    return -float('inf')
                try:
                    return float(item)
                except ValueError:
                    return float('nan')

        eye2letter = {
            ' LEFT': ' L',
            ' RIGHT': ' R'
        }
        color2feature = {  # Converts color to associated feature
            'red': 'Undisturbed Soil',
            'light green': 'Dust-Cleared Rock',
            'light blue': 'Dusty Rock',
            'dark blue': 'Dusty Rock',
            'teal': 'Dusty Rock',
            'light cyan': 'Nodule-rich Rock',
            'dark green': 'Broken Rock Face',
            'light purple': 'Drill Tailings',
            'pink': 'Other',
            'sienna': 'Other',
            'bright red': 'Undisturbed Soil',
            'yellow': 'Other',
            'goldenrod': 'Veins',
            'dark red': 'Disturbed Soil',
            'dark purple': 'Dump Piles'}
        feature2color = {
            'Undisturbed Soil': 'red',
            'Dust-Cleared Rock': 'light green',
            'Dusty Rock': 'light blue',
            'Nodule-rich Rock': 'light cyan',
            'Broken Rock Face': 'dark green',
            'Drill Tailings': 'light purple',
            'Other': 'yellow',
            'Veins': 'goldenrod',
            'Disturbed Soil': 'dark red',
            'Dump Piles': 'dark purple'
        }
        color2RGB = {  # Converts color to associated feature
            'red': (212 / 255, 0, 42 / 255),
            'light green': (109 / 255, 225 / 255, 0),
            'light blue': (0, 0, 225 / 255),
            'dark blue': (0, 0, 111 / 255),
            'teal': (0, 110 / 255, 109 / 255),
            'light cyan': (0, 255 / 255, 255 / 255),
            'dark green': (16 / 255, 113 / 255, 1 / 255),
            'light purple': (254 / 255, 0, 255 / 255),
            'pink': (249 / 255, 104 / 255, 92 / 255),
            'sienna': (143 / 255, 63 / 255, 30 / 255),
            'bright red': (255 / 255, 0, 0),
            'yellow': (255 / 255, 255 / 255, 0),
            'goldenrod': (170 / 255, 117 / 255, 0),
            'dark red': (109 / 255, 0, 0),
            'dark purple': (108 / 255, 0, 110 / 255)}

        # This bit actually scans and saves all the spectra on the server
        list = np.array([])
        for file in os.listdir(path):  # loops through all files in folder
            if file.endswith('.csv') and 'spectra' in file and file.startswith('sol'):  # Finds spectra csv, ignores other files
                df = pd.read_csv(path +'/'+ file)
                a = df.values
                ROIs = df.columns.values.tolist()
                note_index = 0
                float_index = 0
                member_index = 0
                formation_index = 0
                nm = np.array(a[:, 0])
                for i, item in enumerate(nm):
                    if isinstance(item, str):
                        nm[i] = item.strip(' ').lower()
                        if nm[i] == 'notes':
                            note_index = i
                        if nm[i] == 'float':
                            float_index = i
                        if nm[i] == 'member':
                            member_index = i
                        if nm[i] == 'formation':
                            formation_index = i

                while 'notes' in nm or 'float'in nm or 'member' in nm or 'formation' in nm:
                    nm = nm[:-1]

                nm = nm.astype(np.int)
                eye = a[:, ROIs.index(' Eye')]
                mcam = int(file.split('_')[1][4:])
                if mcam in ref_mcam:  # IMPORTANT: mcam must be in in metadata file
                    for title in ROIs:
                        if 'Mean Value' in title:  # Mean value only, not error
                            col = ROIs.index(title)
                            reflectance = np.array(a[:, col])  # last value of reflectance could be notes
                            error = np.array(a[:, col + 1])
                            data = {}
                            metadata = {}
                            for num in np.arange(0, len(nm), 1):
                                data[str(nm[num]) + eye2letter[eye[num]]] = datapoint(nm[num], make_float(
                                    reflectance[num]), make_float(error[num]), eye[num])
                            metadata['Sol'] = int(file.split('_')[0][3:])
                            metadata['Mcam'] = file.split('_')[1][4:]
                            metadata['File'] = file
                            metadata['Color'] = title.split(' Mean Value')[0][1:]
                            metadata['Sequence'] = file.split('_')[0] + '_' + file.split('_')[1]
                            metadata['Feature'] = color2feature[title.split(' Mean Value')[0][1:]]
                            metadata['Name'] = ref_name[get_index(mcam)]
                            metadata['Rover Elevation'] = ref_rover_elevation[get_index(mcam)]
                            metadata['Target Elevation'] = ref_mastcam_elevation[get_index(mcam)]
                            metadata['Tau'] = make_float(ref_tau[get_index(mcam)])
                            metadata['LTST'] = ref_LTST[get_index(mcam)]
                            metadata['Incident'] = make_float(ref_incident[get_index(mcam)])
                            metadata['Emission'] = make_float(ref_emission[get_index(mcam)])
                            metadata['RGB'] = color2RGB[title.split(' Mean Value')[0][1:]]
                            metadata['Ls'] = make_float(ref_ls[get_index(mcam)])
                            metadata['Lat'] = make_float(ref_lat[get_index(mcam)])
                            metadata['Long'] = make_float(ref_long[get_index(mcam)])
                            metadata['Traverse'] = make_float(ref_traverse[get_index(mcam)])
                            metadata['Site'] = make_float(ref_site[get_index(mcam)])
                            metadata['Drive'] = make_float(ref_drive[get_index(mcam)])
                            metadata['R*'] = False
                            metadata['Normalization'] = None

                            metadata['Notes'] = ''
                            if note_index != 0:
                                if not pd.isnull(reflectance[note_index]):
                                    note = str(reflectance[note_index])
                                    if 'outcrop' in note:
                                        note = note.replace('outcrop','rock')
                                    metadata['Notes'] = note
                                    for feature, color in feature2color.items():
                                        if feature.lower() in metadata['Notes'].lower():
                                            metadata['Feature'] = feature
                                            metadata['RGB'] = color2RGB[color]

                            metadata['Float'] = 'No Float'
                            if float_index != 0:
                                if not pd.isnull(reflectance[float_index]):
                                    metadata['Float'] = 'Float'

                            metadata['Tina Member'] = ''
                            metadata['Tina Formation'] = ''
                            if member_index != 0:
                                if not pd.isnull(reflectance[member_index]):
                                    metadata['Tina Member'] = reflectance[member_index].lower()
                            if formation_index != 0:
                                if not pd.isnull(reflectance[formation_index]):
                                    metadata['Tina Formation'] = reflectance[formation_index].lower()

                            spectrum = Spectra(data, metadata)
                            list = np.append(list, spectrum)
                print(file, ' processed')
        return list.tolist()

    def import_spectra(self):
        '''Called from the import spectra menu button. Calls get_spectra_from_directory for each
        directory on the server'''
        # This is the method called when the user wants to update their data from what is currently on the server
        self.import_metadata()#This chunk of code imports the metadata table, which will be added to the spectra later

        #Must be connected to Ricedata server
        if 'win32' in sys.platform:
            path = '//ion.physics.wwu.edu/ricedata/MarsGroup/Mastcam'
        else:
            path = '//Volumes/ricedata/MarsGroup/Mastcam'
        self.spectra_list = []
        for folder in os.listdir(path):  # examines every observations
            if folder.startswith('sol'):
                if 'sol2463' in folder:
                    continue  # Skips this sol, can delete this line later
                address = path +'/'+ folder + '/Working'  # Establishes current path
                self.spectra_list += self.get_spectra_from_directory(address)
        # veins = get_spectra_from_directory(path + '/Veins')
        # veins = [x for x in veins if x.metadata['Feature'] == 'Veins']
        # self.spectra_list.append(veins)
        pickle.dump(self.spectra_list, open('spectra_list.p', 'wb'))  # This saves the data client-side
        print('Done')

    def import_user_spectra(self):
        '''This is the method called when the user wants to update their data from their own spectra. Spectra should
        be stored in the UserSpectra folder'''
        user_spec = self.get_spectra_from_directory('UserSpectra')
        self.user_spectra_list = user_spec
        pickle.dump(user_spec, open('user_spectra.p', 'wb'))  # This saves the user spectra client-side

    def export_csv(self, extract):
        '''Called by the export csv menu command. Behaves very similarly to the plotting function.'''
        def r_star(spectrum):
            new_spec = copy.deepcopy(spectrum)
            if extract.r_star == 0: return new_spec
            if not pd.isnull(new_spec.metadata['Incident']):
                for label, dp in new_spec.data.items():
                    if dp.has_data:
                        dp.reflectance = dp.reflectance / np.cos(np.radians(new_spec.metadata['Incident']))
                        dp.error = dp.error / np.cos(np.radians(new_spec.metadata['Incident']))
                return new_spec
            else:
                return 0

        def normalize(spec):
            point = extract.norm_point
            new_spec = copy.deepcopy(spec)
            if spec == 0: return 0
            left = ''
            right = ''
            if point == 0: # None
                return new_spec
            if point == 1: # 446
                left = '445 L'
                right = '447 R'
            if point == 2: # 494
                left = '495 L'
                right = '493 R'
            if point == 3: # 527
                left = '527 L'
                right = '527 R'
            if point == 4: # 552.5
                left = '554 L'
                right = '551 R'
            if point == 5: # 639
                left = '640 L'
                right = '638 R'
            if point == 6: # 1012.5
                left = '1012 L'
                right = '1013 R'
            if new_spec.has_data(left) and new_spec.has_data(right):
                middle_reflectance = (new_spec.get_reflectance(left) + new_spec.get_reflectance(right)) / 2
                left_factor = middle_reflectance / new_spec.get_reflectance(left)
                right_factor = middle_reflectance / new_spec.get_reflectance(right)
                for filter, dp in new_spec.data.items():
                    if dp.has_data():
                        if dp.eye == ' RIGHT':
                            dp.reflectance *= right_factor
                            dp.error *= right_factor
                        if dp.eye == ' LEFT':
                            dp.reflectance *= left_factor
                            dp.error *= left_factor
                return new_spec
            return 0

        def average_overlapping_bands(spec):
            new_spec = copy.deepcopy(spec)
            if spec == 0: return 0
            if new_spec.has_data('445 L') and new_spec.has_data('447 R'):
                refl = (new_spec.get_reflectance('445 L') + new_spec.get_reflectance('447 R')) / 2
                error = ((new_spec.data['445 L'].error/2)**2 + (new_spec.data['447 R'].error/2)**2)**(.5)
                new_spec.data['446 A'] = datapoint(float(446), refl, error, ' AVERAGE')
            if new_spec.has_data('495 L') and new_spec.has_data('493 R'):
                refl = (new_spec.get_reflectance('495 L') + new_spec.get_reflectance('493 R')) / 2
                error = ((new_spec.data['495 L'].error / 2) ** 2 + (new_spec.data['493 R'].error / 2) ** 2) ** (.5)
                new_spec.data['494 A'] = datapoint(float(494), refl, error, ' AVERAGE')
            if new_spec.has_data('527 L') and new_spec.has_data('527 R'):
                refl = (new_spec.get_reflectance('527 L') + new_spec.get_reflectance('527 R')) / 2
                error = ((new_spec.data['527 L'].error / 2) ** 2 + (new_spec.data['527 R'].error / 2) ** 2) ** (.5)
                new_spec.data['527 A'] = datapoint(float(527), refl, error, ' AVERAGE')
            if new_spec.has_data('554 L') and new_spec.has_data('551 R'):
                refl = (new_spec.get_reflectance('554 L') + new_spec.get_reflectance('551 R')) / 2
                error = ((new_spec.data['554 L'].error / 2) ** 2 + (new_spec.data['551 R'].error / 2) ** 2) ** (.5)
                new_spec.data['552.5 A'] = datapoint(float(552.5), refl, error, ' AVERAGE')
            if new_spec.has_data('640 L') and new_spec.has_data('638 R'):
                refl = (new_spec.get_reflectance('640 L') + new_spec.get_reflectance('638 R')) / 2
                error = ((new_spec.data['640 L'].error / 2) ** 2 + (new_spec.data['638 R'].error / 2) ** 2) ** (.5)
                new_spec.data['639 A'] = datapoint(float(639), refl, error, ' AVERAGE')
            if new_spec.has_data('1012 L') and new_spec.has_data('1013 R'):
                refl = (new_spec.get_reflectance('1012 L') + new_spec.get_reflectance('1013 R')) / 2
                error = ((new_spec.data['1012 L'].error / 2) ** 2 + (new_spec.data['1013 R'].error / 2) ** 2) ** (.5)
                new_spec.data['1012.5 A'] = datapoint(float(1012.5), refl, error, ' AVERAGE')
            return new_spec

        def get_x(spectrum):
            if extract.x_choice == 'Reflectance':
                return spectrum.get_reflectance(extract.x_1)
            if extract.x_choice == 'Slope':
                return spectrum.get_slope(extract.x_1, extract.x_2)
            if extract.x_choice == 'Ratio':
                return spectrum.get_ratio(extract.x_1, extract.x_2)
            if extract.x_choice == 'Band Depth':
                return spectrum.get_band_depth(extract.x_1, extract.x_2, extract.x_3)
            if extract.x_choice == 'Sol':
                return spectrum.metadata['Sol']
            if extract.x_choice == 'LTST':
                LTST = spectrum.metadata['LTST'].split(':')
                return (3600 * float(LTST[0])) + (60 * float(LTST[1])) + float(LTST[2])
            if extract.x_choice == 'Elevation':
                return spectrum.metadata['Elevation']
            if extract.x_choice == 'Tau':
                return spectrum.metadata['Tau']
            if extract.x_choice == 'PCA':
                index = int(extract.x_pca[3]) - 1
                return spectrum.pca_data[index]

        def get_y(spectrum):
            if extract.y_choice == 'Reflectance':
                return spectrum.get_reflectance(extract.y_1)
            if extract.y_choice == 'Slope':
                return spectrum.get_slope(extract.y_1, extract.y_2)
            if extract.y_choice == 'Ratio':
                return spectrum.get_ratio(extract.y_1, extract.y_2)
            if extract.y_choice == 'Band Depth':
                return spectrum.get_band_depth(extract.y_1, extract.y_2, extract.y_3)
            if extract.y_choice == 'Sol':
                return spectrum.metadata['Sol']
            if extract.y_choice == 'LTST':
                LTST = spectrum.metadata['LTST'].split(':')
                return (3600 * float(LTST[0])) + (60 * float(LTST[1])) + float(LTST[2])
            if extract.y_choice == 'Elevation':
                return spectrum.metadata['Elevation']
            if extract.y_choice == 'Tau':
                return spectrum.metadata['Tau']
            if extract.y_choice == 'PCA':
                index = int(extract.y_pca[3]) - 1
                return spectrum.pca_data[index]

        def spectrum_safe(spectrum):
            for wavelength in extract.get_relevant_wavelengths(): #Makes sure that spectrum has all relevant wavelengths
                if spectrum.has_data(wavelength) == False:
                    return False
            if not extract.in_sol_range(spectrum.metadata['Sol']): return False # Checks for sol restriction
            if not extract.in_features(spectrum.metadata['Feature']): return False # Checks for correct features
            if not extract.in_LTST(spectrum.metadata['LTST']): return False # Checks for proper LTST limitation
            if not extract.in_elevation(spectrum.metadata['Elevation']): return False
            if not extract.in_tau(spectrum.metadata['Tau']): return False
            if not extract.in_mcam(spectrum.metadata['Mcam']): return False
            if not extract.in_unit(spectrum.metadata['Member']) and \
                not extract.in_unit(spectrum.metadata['Formation']) and \
                not extract.in_unit(spectrum.metadata['Group']): return False
            if not extract.in_float(spectrum.metadata['Float']): return False
            if extract.x_choice == 'LTST' or extract.y_choice == 'LTST':
                if not spectrum.has_LTST(): return False
            if extract.x_choice == 'Elevation' or extract.y_choice == 'Elevation' or extract.color_type == 2:
                if not spectrum.has_elevation(): return False
            if extract.x_choice == 'Tau' or extract.y_choice == 'Tau':
                if not spectrum.has_tau(): return False
            if extract.x_choice == 'Slope':
                if spectrum.get_slope(extract.x_1, extract.x_2) == 'invalid': return False
            return True

        #The part that applies
        new_list = self.spectra_list[:]
        r_star_list = [r_star(spec) for spec in new_list]
        norm_list = [normalize(spec) for spec in r_star_list]
        avg_list = [average_overlapping_bands(spec) for spec in norm_list]
        newest_list = [value for value in avg_list if value != 0]

        saveas = filedialog.asksaveasfilename(parent=self,
                                              initialdir=os.getcwd(),
                                              title="Please select a file name for saving:",
                                              filetypes=[('CSV', '.csv')],
                                              defaultextension='.csv')

        with open(saveas,'w') as csv:
            csv.truncate()
            csv.write('Sol,Sequence,Target,' + extract.x_label + ',' + extract.y_label + ',Feature Type,Elevation,Closest Tau,LTST,Incidence Angle,Emission Angle\n')

            for spectrum in newest_list:
                if spectrum_safe(spectrum):
                    sequence = str(spectrum.metadata['Mcam'])
                    while len(sequence) < 5:
                        sequence = '0' + sequence

                    row = str(spectrum.metadata['Sol']) + ','
                    row += 'mcam' + sequence + ','
                    row += str(spectrum.metadata['Name'].replace(',','')) + ','
                    row += str(get_x(spectrum)) + ','
                    row += str(get_y(spectrum)) + ','
                    row += str(spectrum.metadata['Feature']) + ','
                    row += str(spectrum.metadata['Elevation']) + ','
                    row += str(spectrum.metadata['Tau']) + ','
                    row += str(spectrum.metadata['LTST']) + ','
                    row += str(spectrum.metadata['Incident']) + ','
                    row += str(spectrum.metadata['Emission']) + '\n'
                    csv.write(row)

    def plot(self, extract):
        '''Lord forgive me for this function. This function is called when the user presses the plot menu option or
        presses the space bar. It processes the data and pulls up a plot. This function also contains all the functions
        that the user can call on that plot, which is why this function is so big.'''

        #Relevant dictionaries used for plotting
        geologic_members={#The Top bound of each member is given, Jura is added seperatley for reasons
            'P. Pnt.':-4170,
            'B. Pnt.':-4210,
            'S. Isl.':-4280,
            'Kara.':-4370,
            'H.V.':-4420,
            'P.H.':-4435,
            'B.G.':-4460
        }
        color2feature = {  # Converts color to associated feature
            'red': 'Undisturbed Soil',
            'light green': 'Dust-Cleared Rock',
            'light blue': 'Dusty Rock',
            'dark blue': 'Dusty Rock',
            'teal': 'Dusty Rock',
            'light cyan': 'Nodule-rich Rock',
            'dark green': 'Broken Rock Face',
            'light purple': 'Drill Tailings',
            'pink': 'Other',
            'sienna': 'Other',
            'bright red': 'Undisturbed Soil',
            'yellow': 'Other',
            'goldenrod': 'Veins',
            'dark red': 'Disturbed Soil',
            'dark purple': 'Dump Piles'}
        color2RGB = {  # Converts color to associated feature
            'red': (212 / 255, 0, 42 / 255),
            'light green': (109 / 255, 225 / 255, 0),
            'light blue': (0, 0, 225 / 255),
            'dark blue': (0, 0, 111 / 255),
            'teal': (0, 110 / 255, 109 / 255),
            'light cyan': (0, 255 / 255, 255 / 255),
            'dark green': (16 / 255, 113 / 255, 1 / 255),
            'light purple': (254 / 255, 0, 255 / 255),
            'pink': (249 / 255, 104 / 255, 92 / 255),
            'sienna': (143 / 255, 63 / 255, 30 / 255),
            'bright red': (255 / 255, 0, 0),
            'yellow': (255 / 255, 255 / 255, 0),
            'goldenrod': (170 / 255, 117 / 255, 0),
            'dark red': (109 / 255, 0, 0),
            'dark purple': (108 / 255, 0, 110 / 255)}
        RGB2color = dict([[v,k] for k,v in color2RGB.items()])

        #Define functions relevant for plotting
        def r_star(spectrum):
            new_spec = copy.deepcopy(spectrum)
            if extract.r_star == 0: return new_spec
            if not pd.isnull(new_spec.metadata['Incident']):
                for label, dp in new_spec.data.items():
                    if dp.has_data:
                        dp.reflectance = dp.reflectance / np.cos(np.radians(new_spec.metadata['Incident']))
                        dp.error = dp.error / np.cos(np.radians(new_spec.metadata['Incident']))
                new_spec.metadata['R*'] = True
                return new_spec
            else:
                return 0

        def normalize(spec):
            point = extract.norm_point
            new_spec = copy.deepcopy(spec)
            if spec == 0: return 0
            left = ''
            right = ''
            if point == 0: # None
                return new_spec
            if point == 1: # 446
                left = '445 L'
                right = '447 R'
                norm_point = '446 nm'
            if point == 2: # 494
                left = '495 L'
                right = '493 R'
                norm_point = '494 nm'
            if point == 3: # 527
                left = '527 L'
                right = '527 R'
                norm_point = '527 nm'
            if point == 4: # 552.5
                left = '554 L'
                right = '551 R'
                norm_point = '552.5 nm'
            if point == 5: # 639
                left = '640 L'
                right = '638 R'
                norm_point = '649 nm'
            if point == 6: # 1012.5
                left = '1012 L'
                right = '1013 R'
                norm_point = '1012.5 nm'
            if new_spec.has_data(left) and new_spec.has_data(right):
                middle_reflectance = (new_spec.get_reflectance(left) + new_spec.get_reflectance(right)) / 2
                left_factor = middle_reflectance / new_spec.get_reflectance(left)
                right_factor = middle_reflectance / new_spec.get_reflectance(right)
                for filter, dp in new_spec.data.items():
                    if dp.has_data():
                        if dp.eye == ' RIGHT':
                            dp.reflectance *= right_factor
                            dp.error *= right_factor
                        if dp.eye == ' LEFT':
                            dp.reflectance *= left_factor
                            dp.error *= left_factor
                new_spec.metadata['Normalization'] = norm_point
                return new_spec
            return 0

        def average_overlapping_bands(spec):
            new_spec = copy.deepcopy(spec)
            if spec == 0: return 0
            if new_spec.has_data('445 L') and new_spec.has_data('447 R'):
                refl = (new_spec.get_reflectance('445 L') + new_spec.get_reflectance('447 R')) / 2
                error = ((new_spec.data['445 L'].error/2)**2 + (new_spec.data['447 R'].error/2)**2)**(.5)
                new_spec.data['446 A'] = datapoint(float(446), refl, error, ' AVERAGE')
            if new_spec.has_data('495 L') and new_spec.has_data('493 R'):
                refl = (new_spec.get_reflectance('495 L') + new_spec.get_reflectance('493 R')) / 2
                error = ((new_spec.data['495 L'].error / 2) ** 2 + (new_spec.data['493 R'].error / 2) ** 2) ** (.5)
                new_spec.data['494 A'] = datapoint(float(494), refl, error, ' AVERAGE')
            if new_spec.has_data('527 L') and new_spec.has_data('527 R'):
                refl = (new_spec.get_reflectance('527 L') + new_spec.get_reflectance('527 R')) / 2
                error = ((new_spec.data['527 L'].error / 2) ** 2 + (new_spec.data['527 R'].error / 2) ** 2) ** (.5)
                new_spec.data['527 A'] = datapoint(float(527), refl, error, ' AVERAGE')
            if new_spec.has_data('554 L') and new_spec.has_data('551 R'):
                refl = (new_spec.get_reflectance('554 L') + new_spec.get_reflectance('551 R')) / 2
                error = ((new_spec.data['554 L'].error / 2) ** 2 + (new_spec.data['551 R'].error / 2) ** 2) ** (.5)
                new_spec.data['552.5 A'] = datapoint(float(552.5), refl, error, ' AVERAGE')
            if new_spec.has_data('640 L') and new_spec.has_data('638 R'):
                refl = (new_spec.get_reflectance('640 L') + new_spec.get_reflectance('638 R')) / 2
                error = ((new_spec.data['640 L'].error / 2) ** 2 + (new_spec.data['638 R'].error / 2) ** 2) ** (.5)
                new_spec.data['639 A'] = datapoint(float(639), refl, error, ' AVERAGE')
            if new_spec.has_data('1012 L') and new_spec.has_data('1013 R'):
                refl = (new_spec.get_reflectance('1012 L') + new_spec.get_reflectance('1013 R')) / 2
                error = ((new_spec.data['1012 L'].error / 2) ** 2 + (new_spec.data['1013 R'].error / 2) ** 2) ** (.5)
                new_spec.data['1012.5 A'] = datapoint(float(1012.5), refl, error, ' AVERAGE')
            return new_spec

        def get_x(spectrum):
            if extract.x_choice == 'Reflectance':
                return spectrum.get_reflectance(extract.x_1)
            if extract.x_choice == 'Slope':
                return spectrum.get_slope(extract.x_1, extract.x_2)
            if extract.x_choice == 'Ratio':
                return spectrum.get_ratio(extract.x_1, extract.x_2)
            if extract.x_choice == 'Band Depth':
                return spectrum.get_band_depth(extract.x_1, extract.x_2, extract.x_3)
            if extract.x_choice == 'Sol':
                return spectrum.metadata['Sol']
            if extract.x_choice == 'LTST':
                LTST = spectrum.metadata['LTST'].split(':')
                return (3600 * float(LTST[0])) + (60 * float(LTST[1])) + float(LTST[2])
            if extract.x_choice == 'Elevation':
                return spectrum.metadata['Elevation']
            if extract.x_choice == 'Tau':
                return spectrum.metadata['Tau']
            if extract.x_choice == 'PCA':
                index = int(extract.x_pca[3]) - 1
                return spectrum.pca_data[index]
            if extract.x_choice == 'Ls':
                return spectrum.metadata['Ls']
            if extract.x_choice == 'Latitude':
                return spectrum.metadata['Lat']
            if extract.x_choice == 'Longitude':
                return spectrum.metadata['Long']
            if extract.x_choice == 'Traverse Distance':
                return spectrum.metadata['Traverse']

        def get_y(spectrum):
            if extract.y_choice == 'Reflectance':
                return spectrum.get_reflectance(extract.y_1)
            if extract.y_choice == 'Slope':
                return spectrum.get_slope(extract.y_1, extract.y_2)
            if extract.y_choice == 'Ratio':
                return spectrum.get_ratio(extract.y_1, extract.y_2)
            if extract.y_choice == 'Band Depth':
                return spectrum.get_band_depth(extract.y_1, extract.y_2, extract.y_3)
            if extract.y_choice == 'Sol':
                return spectrum.metadata['Sol']
            if extract.y_choice == 'LTST':
                LTST = spectrum.metadata['LTST'].split(':')
                return (3600 * float(LTST[0])) + (60 * float(LTST[1])) + float(LTST[2])
            if extract.y_choice == 'Elevation':
                return spectrum.metadata['Elevation']
            if extract.y_choice == 'Tau':
                return spectrum.metadata['Tau']
            if extract.y_choice == 'PCA':
                index = int(extract.y_pca[3]) - 1
                return spectrum.pca_data[index]
            if extract.y_choice == 'Ls':
                return spectrum.metadata['Ls']
            if extract.y_choice == 'Latitude':
                return spectrum.metadata['Lat']
            if extract.y_choice == 'Longitude':
                return spectrum.metadata['Long']
            if extract.y_choice == 'Traverse Distance':
                return spectrum.metadata['Traverse']

        def get_x_error(spectrum):
            if extract.x_choice == 'Reflectance':
                return spectrum.get_reflectance_error(extract.x_1)
            if extract.x_choice == 'Slope':
                return spectrum.get_slope_error(extract.x_1, extract.x_2)
            if extract.x_choice == 'Ratio':
                return spectrum.get_ratio_error(extract.x_1, extract.x_2)
            if extract.x_choice == 'Band Depth':
                return spectrum.get_band_depth_error(extract.x_1, extract.x_2, extract.x_3)
            return 0

        def get_y_error(spectrum):
            if extract.y_choice == 'Reflectance':
                return spectrum.get_reflectance_error(extract.y_1)
            if extract.y_choice == 'Slope':
                return spectrum.get_slope_error(extract.y_1, extract.y_2)
            if extract.y_choice == 'Ratio':
                return spectrum.get_ratio_error(extract.y_1, extract.y_2)
            if extract.y_choice == 'Band Depth':
                return spectrum.get_band_depth_error(extract.y_1, extract.y_2, extract.y_3)
            return 0

        def get_color(spectrum):
            if extract.color_type == 0:
                if spectrum_highlighted(spectrum): return 'r'
                return 'k'
            if extract.color_type == 1:
                return spectrum.metadata['RGB']
            if extract.color_type == 2:
                return cmap(norm(spectrum.metadata['Elevation']))
            if extract.color_type == 3:
                return spectrum.metadata['SimpleRGB']
            if extract.color_type == 4:
                return cmap_sol(norm_sol(spectrum.metadata['Sol']))
            if extract.color_type == 5:
                return cmap_tau(norm_tau(spectrum.metadata['Tau']))
            if extract.color_type == 6:
                return cmap_ls(norm_ls(spectrum.metadata['Ls']))

        def get_label(spectrum):
            if extract.color_type == 3:
                return spectrum.metadata['SimpleFeature']
            return spectrum.metadata['Feature']

        def spectrum_safe(spectrum):
            '''Used to determine is a given spectrum should be plotted or not. Uses the in ____ function in the extract
            object to see if a given spectrum is safe.'''
            for wavelength in extract.get_relevant_wavelengths(): #Makes sure that spectrum has all relevant wavelengths
                if spectrum.has_data(wavelength) == False:
                    return False
            if not extract.in_sol_range(spectrum.metadata['Sol']): return False # Checks for sol restriction
            if not extract.in_features(spectrum.metadata['Feature']): return False # Checks for correct features
            if not extract.in_LTST(spectrum.metadata['LTST']): return False # Checks for proper LTST limitation
            # if not extract.in_elevation(spectrum.metadata['Elevation']): return False
            if not extract.in_tau(spectrum.metadata['Tau']): return False
            if not extract.in_ls(spectrum.metadata['Ls']): return False
            if not extract.in_mcam(spectrum.metadata['Mcam']): return False
            if not extract.in_unit(spectrum.metadata['Member']) and \
                not extract.in_unit(spectrum.metadata['Formation']) and \
                not extract.in_unit(spectrum.metadata['Group']): return False
            if not extract.in_float(spectrum.metadata['Float']): return False
            if not extract.in_longitude(spectrum.metadata['Long']): return False
            if not extract.in_latitude(spectrum.metadata['Lat']): return False
            if not extract.in_traverse(spectrum.metadata['Traverse']): return False
            if extract.x_choice == 'LTST' or extract.y_choice == 'LTST':
                if not spectrum.has_LTST(): return False
            if extract.x_choice == 'Elevation' or extract.y_choice == 'Elevation' or extract.color_type == 2:
                if not spectrum.has_elevation(): return False
            if extract.x_choice == 'Tau' or extract.y_choice == 'Tau':
                if not spectrum.has_tau(): return False
            if extract.x_choice == 'Slope' or extract.y_choice == 'Slope':
                if spectrum.get_slope(extract.x_1, extract.x_2) == 'invalid': return False
            if extract.x_choice == 'Ls' or extract.y_choice == 'Ls':
                if not spectrum.has_ls(): return False
            if extract.x_choice == 'Latitude' or extract.y_choice == 'Latitude':
                if not spectrum.has_lat(): return False
            if extract.x_choice == 'Longitude' or extract.y_choice == 'Longitude':
                if not spectrum.has_long(): return False
            if extract.x_choice == 'Traverse Distance' or extract.y_choice == 'Traverse Distance':
                if not spectrum.has_traverse(): return False
            return True

        def spectrum_highlighted(spectrum):
            '''Used to determine if a plotted spectrum should be highlighted or not. Uses the extraction object
            in_highlighted_______ to determine is a spectrum should be highlighted.'''
            if extract.in_highlighted_sol_range(spectrum.metadata['Sol']): return True # Checks for sol restriction
            if extract.in_highlighted_features(spectrum.metadata['Feature']): return True # Checks for correct features
            if extract.in_highlighted_LTST(spectrum.metadata['LTST']): return True # Checks for proper LTST limitation
            # if extract.in_highlighted_elevation(spectrum.metadata['Elevation']): return True
            if extract.in_highlighted_tau(spectrum.metadata['Tau']): return True
            if extract.in_highlighted_mcam(spectrum.metadata['Mcam']): return True
            if extract.in_highlighted_unit(spectrum.metadata['Member']): return True
            if extract.in_highlighted_unit(spectrum.metadata['Formation']): return True
            if extract.in_highlighted_unit(spectrum.metadata['Group']): return True
            if extract.in_highlighted_float(spectrum.metadata['Float']): return True
            if extract.in_highlighted_ls(spectrum.metadata['Ls']): return True
            if extract.in_highlighted_longitude(spectrum.metadata['Lat']): return True
            if extract.in_highlighted_latitude(spectrum.metadata['Long']): return True
            if extract.in_highlighted_traverse(spectrum.metadata['Traverse']): return True
            return False

        def save_plot(figure):
            '''Used to save a plot as a png'''
            saveas = filedialog.asksaveasfilename(parent=self,
                                                  initialdir=os.getcwd(),
                                                  title="Please select a file name for saving:",
                                                  filetypes=[('PNG', '.png'),('JPEG','.jpg')],
                                                  defaultextension='.png')
            figure.savefig(saveas)

        def identify():
            '''Used to identify a point on the plot. Identifies the point, then pulls up another window with all the
            data and metadata.'''
            #Calculates the clicked point
            user_input = fig.ginput()
            new_x_values = x_values - user_input[0][0]
            new_y_values = y_values - user_input[0][1]
            quad = np.sqrt(np.add(np.power(new_x_values, 2), np.power(new_y_values, 2)))
            clicked_spec = plotted_spec[np.argmin(quad)]  # The clicked spectrum

            #Creates window and populates with metadata
            i = tk.Toplevel()
            i.wm_title('Point Identification')
            text = tk.Text(i, height = 10, padx = 10)
            text.insert(tk.INSERT, 'X: ' + str(get_x(clicked_spec)) + '\n')
            text.insert(tk.INSERT, 'Y: ' + str(get_y(clicked_spec)) + '\n')
            for meta, value in clicked_spec.metadata.items():
                text.insert(tk.INSERT, meta + ': '+str(value)+'\n')
            text.grid(row=1, column=1)

            images = []
            for image in os.listdir(self.home_path + '/ROI_Images'):
                if image.startswith(clicked_spec.metadata['File'].split('_spectra.csv')[0]):
                    images.append(self.home_path + '/ROI_Images/' + image)

            if len(images) >= 1:#only try to display images if they exist
                def switch_pic():
                    self.panel.grid_forget()
                    size = 128, 128
                    img = Image.open(images[self.num])
                    img = img.resize((600,600))
                    img = ImageTk.PhotoImage(img)
                    self.panel = tk.Label(i, image = img)
                    self.panel.image = img
                    self.panel.bind("<Button-1>", lambda self: switch_pic())
                    self.panel.grid(row=1, column=2, rowspan=2)
                    self.num += 1
                    if self.num == len(images): self.num = 0

                self.num = 0
                self.panel = tk.Label(i, text='Test')
                self.panel.grid(row=1, column=2, rowspan=2)
                switch_pic()

            # Plot spectrum
            def get_data_points(spec):
                if spec == 0: return 0
                overlapping_filters = [['445 L', '447 R', '446 A'],
                                       ['495 L', '493 R', '494 A'],
                                       ['527 L', '527 R', '527 A'],
                                       ['554 L', '551 R', '552.5 A'],
                                       ['640 L', '638 R', '639 A'],
                                       ['1012 L', '1013 R', '1012.5 A']]
                single_filters = ['676 L', '751 L', '805 R', '867 L', '908 R', '937 R']
                datapoints = []
                for trio in overlapping_filters:
                    left = trio[0]
                    right = trio[1]
                    average = trio[2]
                    if spec.has_data(average):
                        datapoints.append(spec.data[average])
                    else:
                        if spec.has_data(left):
                            datapoints.append(spec.data[left])
                        if spec.has_data(right):
                            datapoints.append(spec.data[right])
                for filter in single_filters:
                    if spec.has_data(filter):
                        datapoints.append(spec.data[filter])
                return sorted(datapoints, key=lambda dp: dp.filter)

            #Creates spectra plot
            fig2 = mpl.figure.Figure()
            plt2 = fig2.add_subplot(1,1,1)
            canvas2 = tkagg.FigureCanvasTkAgg(fig2, master=i)
            canvas2.get_tk_widget().grid(row = 2, column = 1)

            dp_list = get_data_points(clicked_spec)
            wavelengths = [dp.filter for dp in dp_list]
            reflectances = [dp.reflectance for dp in dp_list]
            errors = [dp.error for dp in dp_list]

            plt2.errorbar(wavelengths, reflectances, yerr=errors, color=clicked_spec.metadata['RGB'], fmt='-o')
            plt2.set_title('Spectrum', fontsize = 15)
            plt2.set_xlabel('Wavelength (nm)', fontsize = 15)
            if extract.r_star == 1:
                plt2.set_ylabel('R* Reflectance', fontsize = 15)
            else: plt2.set_ylabel('IOF', fontzise = 15)
            canvas2.draw()

        def label_point():
            '''Labels a click point. Finds out which point was clicked, then updates the plot with its label.'''
            user_input = fig.ginput()
            new_x_values = x_values - user_input[0][0]
            new_y_values = y_values - user_input[0][1]
            quad = np.sqrt(np.add(np.power(new_x_values, 2), np.power(new_y_values, 2)))
            clicked_spec = plotted_spec[np.argmin(quad)]

            plt.annotate('   ' + str(clicked_spec.metadata['Name']), xycoords='data',
                         xy=(get_x(clicked_spec), get_y(clicked_spec)), fontsize=7)
            canvas.draw()

        def resize():
            '''Used to resize the plot. User clicks two corner plots, and then the plot is updated'''
            user_input = fig.ginput(2)
            x1 = user_input[0][0]
            x2 = user_input[1][0]
            y1 = user_input[0][1]
            y2 = user_input[1][1]
            if x1 > x2:
                x1, x2 = x2, x1
            if y1 > y2:
                y1, y2, = y2, y1
            plt.set_xlim([x1, x2])
            plt.set_ylim([y1, y2])
            canvas.draw()

        def reset():
            ''' Used to reset plot to original size'''
            plt.set_xlim(x_lim)
            plt.set_ylim(y_lim)
            canvas.draw()

        def bigger_font():
            '''Increases font size of the plot.'''
            for item in ([plt.title, plt.xaxis.label, plt.yaxis.label] +
                             plt.get_xticklabels() + plt.get_yticklabels()):
                item.set_fontsize(20)
            canvas.draw()

        def show_PCs():
            '''Calculates and shows all the principle components of current data distribution.'''
            # Initialize window
            i = tk.Toplevel()
            i.wm_title('Principle Components')

            # Run PCA. FIrst the data is sorted in this loop
            data = []
            for spec in safe_list:
                wavelengths = []
                for key, dp in spec.data.items():
                    if key in pca_filters:
                        wavelengths.append([dp.filter, dp.reflectance])
                wavelengths = sorted(wavelengths, key=lambda pair: pair[0])
                wavelengths = [c[1] for c in wavelengths if not pd.isnull(c[1])]
                if len(wavelengths) == 12:
                    data.append(wavelengths)
            pca = PCA()
            pca.fit(data)

            #Make plots
            fig3, axes3 = pyplot.subplots(4, 3, figsize=(12,8), sharex=True, sharey=True)
            canvas3 = tkagg.FigureCanvasTkAgg(fig3, master=i)
            canvas3.get_tk_widget().grid(row=0, column=0)
            wavelengths = [446, 494, 527, 552, 639, 676, 751, 805, 867, 908, 937, 1012]
            variance = pca.explained_variance_ratio_
            for i, PC in enumerate(pca.components_):  # For each PC
                ax = axes3[int(i / 3), i % 3]
                ax.plot(wavelengths,PC,'-o')
                ax.set_ylim([-.65,.65])
                ax.set_title('PC ' + str(i+1) + ', Varience: ' + str(round(variance[i],8)))
                if int(i / 3) == 3:
                    ax.set_xlabel('Wavelength')
                if i % 3 == 0:
                    ax.set_ylabel('Reflectance')
            pyplot.tight_layout()
            canvas3.draw()

        def draw_line_prompt():
            '''Pulls up a dialog box to figure out what lines need to be drawn. Makes a dictionary containing
            line data, and pass it to the draw_line_plot function'''
            # Used for fit lines, diagonals, etc.
            j = tk.Toplevel()  # New Window
            j.wm_title('Line Options')

            # Horizontal line
            horizontal_bool = tk.BooleanVar(value=False)
            horizontal_button = tk.Checkbutton(j, variable=horizontal_bool,text='Horizontal Line')
            horizontal_button.grid(row=1, column=1, sticky='W')
            horizontal_input = tk.Entry(j)
            horizontal_input.insert(0, 'Y position')
            horizontal_input.grid(row=1, column=2, sticky='W')

            # Vertical line
            vertical_bool = tk.BooleanVar(value=False)
            vertical_button = tk.Checkbutton(j, variable=vertical_bool, text='Vertical Line')
            vertical_button.grid(row=2, column=1, sticky='W')
            vertical_input = tk.Entry(j)
            vertical_input.insert(0, 'X position')
            vertical_input.grid(row=2, column=2, sticky='W')

            # Diagonal line
            diagonal_bool = tk.BooleanVar(value=False)
            diagonal_button = tk.Checkbutton(j, variable=diagonal_bool, text='Diagonal Line')
            diagonal_button.grid(row=3, column=1, sticky='W')

            # Polynomial fit
            polynomial_bool = tk.BooleanVar(value=False)
            polynomial_button = tk.Checkbutton(j, variable=polynomial_bool, text='Polynomial Fit')
            polynomial_button.grid(row=4, column=1, sticky='W')
            polynomial_input = tk.Entry(j)
            polynomial_input.insert(0, 'Polynomial Degree')
            polynomial_input.grid(row=4, column=2, sticky='W')

            # Go button
            # Button build info dictionary and sends it to draw_line_plot
            def button_line_execute():
                line_dict = {}
                line_dict['PlotHor'] = horizontal_bool.get()
                line_dict['PlotVert'] = vertical_bool.get()
                line_dict['PlotDiag'] = diagonal_bool.get()
                line_dict['PlotPoly'] = polynomial_bool.get()
                line_dict['X'] = vertical_input.get()
                line_dict['Y'] = horizontal_input.get()
                line_dict['N'] = polynomial_input.get()
                draw_line_plot(line_dict)
            plot_lines = tk.Button(j, text='Plot Lines', command= lambda: button_line_execute())
            plot_lines.grid(row=5,column = 1, columnspan=2)

        def draw_line_plot(line_dict):
            '''Takes a line dictionary and plots the relevant lines'''
            if line_dict['PlotHor']:
                plt.axhline(y=float(line_dict['Y']), color='k', ls='--', zorder=0)
            if line_dict['PlotVert']:
                plt.axvline(x=float(line_dict['X']), color='k', ls='--', zorder=0)
            if line_dict['PlotDiag']:
                bounds_x = plt.get_xlim()
                plt.plot(bounds_x, bounds_x, color='k', ls='--', zorder=0)
            if line_dict['PlotPoly']:
                p = np.polyfit(x_values, y_values, int(line_dict['N'])) # Get an N degree polynomial that fits data
                polynomial = np.poly1d(p)
                bounds_x = plt.get_xlim()
                x = np.linspace(bounds_x[0], bounds_x[1], 500)
                plt.plot(x, polynomial(x), color='k', zorder=0)

                # Calculating R^2
                mean = np.average(y_values)
                mean_total = 0
                fit_total = 0
                for i, x in enumerate(x_values):
                    mean_difference = y_values[i] - mean
                    mean_total += mean_difference**2
                    fit_difference = y_values[i] - polynomial(x)
                    fit_total += fit_difference**2
                r_2 = (mean_total - fit_total) / mean_total
                weights = ', '.join(np.around(p, 4).astype(str)) # Weights of the polynomial
                label = 'R^2 Value: ' + str(r_2) + '\n' + 'Polynomial weights (high to low): ' + weights
                plt.annotate(label, xycoords='axes fraction', xy=(0.01,.95)) # Display R^2
            canvas.draw()

        spec_collect = False
        spec_collect_list = []
        def plot_multiple_spectra():
            spec_collect = not spec_collect
            while spec_collect:
                user_input = fig.ginput()
                new_x_values = x_values - user_input[0][0]
                new_y_values = y_values - user_input[0][1]
                quad = np.sqrt(np.add(np.power(new_x_values, 2), np.power(new_y_values, 2)))
                clicked_spec = plotted_spec[np.argmin(quad)]  # The clicked spectrum
                spec_collect_list.append(clicked_spec)
            if not spec_collect:
                print(spec_collect_list)


        # This is all the stuff that gets the plotting figure set up
        t = tk.Toplevel()
        t.wm_title('Plot')
        x = 11
        y = 8
        fig = mpl.figure.Figure(figsize=(x, y))
        plt = fig.add_axes([.2 * y / x, .1, .8 * y / x, .8])
        canvas = tkagg.FigureCanvasTkAgg(fig, master=t)
        canvas.get_tk_widget().pack()

        menubar = tk.Menu(t)
        menubar.add_command(label='Save Plot', command=lambda: save_plot(fig))
        menubar.add_command(label='Identify Point', command=lambda: identify())
        menubar.add_command(label='Resize Plot', command=lambda:resize())
        menubar.add_command(label='Reset Plot', command=lambda:reset())
        menubar.add_command(label='Increase Font Size', command=lambda: bigger_font())
        menubar.add_command(label='Label Point', command=lambda: label_point())
        menubar.add_command(label='Principle Components', command=lambda: show_PCs())
        menubar.add_command(label='Draw Lines', command=lambda: draw_line_prompt())
        t.config(menu=menubar)

        # Bind commands to keys. Neccesary since the menubar above doesn't show up on mac for some reason
        t.bind('i', lambda x: identify())
        t.bind('s', lambda x: save_plot(fig))
        t.bind('l', lambda x: label_point())
        t.bind('z', lambda x: resize())
        t.bind('b', lambda x: bigger_font())
        t.bind('r', lambda x: reset())
        t.bind('p', lambda x: show_PCs())
        t.bind('f', lambda x: draw_line_prompt())
        t.bind('m', lambda x: plot_multiple_spectra())

        #This section determines which spectra to consider
        self.new_list = []
        if extract.user_spectra == 1:
            self.new_list+=self.get_spectra_from_directory('UserSpectra')
        if extract.our_spectra == 1:
            self.new_list+=list(self.spectra_list)

        # This section applies corrections to spectra
        r_star_list = [r_star(spec) for spec in self.new_list]  # R* correct
        norm_list = [normalize(spec) for spec in r_star_list]  # Normalize
        avg_list = [average_overlapping_bands(spec) for spec in norm_list]  # Average overlapping
        newest_list = [value for value in avg_list if value != 0]  # Gets rid of dead ones
        safe_list = [spec for spec in newest_list if spectrum_safe(spec)]  # Ensures that spectra match user conditions

        #Defines Elevation color map
        # elevation = [spec.metadata['Elevation'] for spec in newest_list if spec.has_elevation()]
        # norm = mpl.colors.Normalize(vmin=min(elevation), vmax=max(elevation))
        # cmap = mpl.cm.nipy_spectral

        # Defines Sol color map
        sol_list = [spec.metadata['Sol'] for spec in newest_list]
        norm_sol = mpl.colors.Normalize(vmin=min(sol_list), vmax=max(sol_list))
        cmap_sol = mpl.cm.nipy_spectral

        # Defines Ls color map
        ls_list = [spec.metadata['Ls'] for spec in newest_list]
        norm_ls = mpl.colors.Normalize(vmin=min(ls_list), vmax=max(ls_list))
        cmap_ls = mpl.cm.nipy_spectral

        # Defines Tau color map
        tau_list = [spec.metadata['Tau'] for spec in newest_list if spec.has_tau()]
        norm_tau = mpl.colors.Normalize(vmin=min(tau_list), vmax=max(tau_list))
        cmap_tau = mpl.cm.nipy_spectral

        # Sets X and Y plot labels
        plt.set_xlabel(extract.x_label, fontsize = 15)
        plt.set_ylabel(extract.y_label, fontsize = 15)
        plt.tick_params(which='both', labelsize=12)

        # Assigns plotting shape, size, and priority to each spectrum
        for spectrum in safe_list:
            if spectrum_highlighted(spectrum): #If item is highlighted
                spectrum.symbol = 'D'
                spectrum.priority = 10
            else:
                if extract.userhighlight == True:
                    spectrum.symbol = '.'
                    spectrum.priority = 1
                else:
                    spectrum.symbol = 'o'
                    spectrum.priority = 1

        # This part runs PCA, if needed
        pca_filters = ['446 A', '494 A', '527 A', '552.5 A', '639 A', '676 L',
                       '751 L', '805 R', '867 L', '908 R', '937 R', '1012.5 A']
        if extract.x_choice == 'PCA' or extract.y_choice == 'PCA':
            # This organizes the data in a way the PCA code can read
            data = []
            for spec in safe_list:
                wavelengths = []
                for key, dp in spec.data.items():
                    if key in pca_filters:
                        wavelengths.append([dp.filter, dp.reflectance])
                wavelengths = sorted(wavelengths, key=lambda pair: pair[0])
                wavelengths = [c[1] for c in wavelengths if  not pd.isnull(c[1])]
                if len(wavelengths) == 12:
                    data.append(wavelengths)
            # This runs the PCA on the data
            pca = PCA()
            pca.fit(data)
            plot = pca.transform(data)
            for i, spec in enumerate(safe_list):
                spec.pca_data = plot[i]  # PC weights are assigned to each spectra object

        # The part that actually plots!
        used_colors = []
        x_values = []
        y_values = []
        plotted_spec = []
        for spectrum in safe_list:
            if not get_color(spectrum) in used_colors and extract.color_type in [1,3]:
                plt.plot(get_x(spectrum), get_y(spectrum), spectrum.symbol, color = get_color(spectrum),
                         label=get_label(spectrum), markeredgecolor = 'black',
                         markeredgewidth=.5, zorder = spectrum.priority)
                used_colors.append(get_color(spectrum))
            else: plt.plot(get_x(spectrum), get_y(spectrum), spectrum.symbol, color = get_color(spectrum),
                           markeredgecolor='black', markeredgewidth=.5, zorder=spectrum.priority)
            x_values.append(get_x(spectrum))
            y_values.append(get_y(spectrum))
            plotted_spec.append(spectrum)

        #Error bar plotting
        if extract.propogate_error == 1:
            for spectrum in newest_list:
                if spectrum_safe(spectrum):
                    plt.errorbar(get_x(spectrum), get_y(spectrum), xerr = get_x_error(spectrum),
                                 yerr=get_y_error(spectrum),fmt='none',color=get_color(spectrum))

        if extract.propogate_error == 2:
            x_avgs = []
            y_avgs = []
            for spectrum in newest_list:
                if spectrum_safe(spectrum):
                    x_avgs.append(get_x_error(spectrum))
                    y_avgs.append(get_y_error(spectrum))
            x_avg = np.median(x_avgs)
            y_avg = np.median(y_avgs)
            plt2 = fig.add_axes([y/x,.6,(x-(1.1*y))/x,.3])
            x_bounds = np.array(plt.get_xlim()) * (x-(1.1*y))/.8/y
            x_bounds -= (x_bounds[1] + x_bounds[0])/2
            plt2.set_xlim(x_bounds.tolist())
            y_bounds = np.array(plt.get_ylim())*.3/.8
            y_bounds -= (y_bounds[1] + y_bounds[0]) / 2
            plt2.set_ylim(y_bounds.tolist())
            plt2.errorbar(0, 0, xerr=x_avg, yerr=y_avg, fmt='none', color='k')
            plt2.set_title('Representative Error Bar')
            plt2.tick_params(
                axis='both',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,          # ticks along the top edge are off
                left=False,
                labelbottom=False,
                labelleft = False)

        #Deals with legends
        if extract.color_type in [1, 3]:
            plt.legend()
            plt.legend(bbox_to_anchor=(1.05, .5), loc=2, borderaxespad=0.)

        if extract.color_type == 2: # Run if color based off of elevation
            fake = fig.add_axes([1.5, 1.5, 0.075, 0.8]) # A dummy subplot needed for color bar
            im = fake.imshow([[0,0],[0,0]], cmap='nipy_spectral')
            cax = fig.add_axes([0.8, 0.1, 0.075, 0.4])
            colorbar = fig.colorbar(im,cax=cax,cmap=cmap,norm=norm,ticks=[-.1, 0, .1])
            colorbar.ax.set_yticklabels([str(int(min(elevation)))+' m',str(int((min(elevation)+max(elevation))/2))+' m',str(int(max(elevation)))+' m'])
            cax.text(0,.5,'Elevation',
                horizontalalignment='right',
                verticalalignment='center',
                rotation='vertical',
                transform=cax.transAxes,
                fontsize = 15)

        if extract.color_type == 4: # Run if color based off of Sol
            fake = fig.add_axes([1.5, 1.5, 0.075, 0.8]) # A dummy subplot needed for color bar
            im = fake.imshow([[0,0],[0,0]], cmap='nipy_spectral')
            cax = fig.add_axes([0.8, 0.1, 0.075, 0.4])
            colorbar = fig.colorbar(im,cax=cax,cmap=cmap_sol,norm=norm_sol,ticks=[-.1, .1])
            colorbar.ax.set_yticklabels([str(int(min(sol_list))),str(int(max(sol_list)))])
            cax.text(0,.5,'Sol',
                horizontalalignment='right',
                verticalalignment='center',
                rotation='vertical',
                transform=cax.transAxes,
                fontsize = 15)

        if extract.color_type == 5: # Run if color based off of Tau
            fake = fig.add_axes([1.5, 1.5, 0.075, 0.8]) # A dummy subplot needed for color bar
            im = fake.imshow([[0,0],[0,0]], cmap='nipy_spectral')
            cax = fig.add_axes([0.8, 0.1, 0.075, 0.4])
            colorbar = fig.colorbar(im,cax=cax,cmap=cmap_tau,norm=norm_tau,ticks=[-.1, .1])
            colorbar.ax.set_yticklabels([str(min(tau_list)),str(max(tau_list))])
            cax.text(0,.5,'Tau',
                horizontalalignment='right',
                verticalalignment='center',
                rotation='vertical',
                transform=cax.transAxes,
                fontsize = 15)


        if extract.color_type == 6: # Run if color based off of Ls
            fake = fig.add_axes([1.5, 1.5, 0.075, 0.8]) # A dummy subplot needed for color bar
            im = fake.imshow([[0,0],[0,0]], cmap='nipy_spectral')
            cax = fig.add_axes([0.8, 0.1, 0.075, 0.4])
            colorbar = fig.colorbar(im,cax=cax,cmap=cmap_ls,norm=norm_ls,ticks=[-.1, .1])
            colorbar.ax.set_yticklabels([str(int(min(ls_list))),str(int(max(ls_list)))])
            cax.text(0,.5,'Ls',
                horizontalalignment='right',
                verticalalignment='center',
                rotation='vertical',
                transform=cax.transAxes,
                fontsize = 15)

        #Plots geologic member boundaries
        if extract.geologic_members == 1:
            if extract.y_choice == 'Elevation':
                for boundary, elev in geologic_members.items():
                    plt.axhline(elev, dashes = (10,10), color = 'gray', zorder = 0)
                    plt.text(plt.get_xlim()[1], elev - (plt.get_ylim()[1]-plt.get_ylim()[0])*.02, boundary, rotation = 270, fontsize = 8)
                plt.text(plt.get_xlim()[1], -4150 - (plt.get_ylim()[1] - plt.get_ylim()[0]) * .02, 'Jura',
                         rotation=270, fontsize=8)

        #Shows 10:30, noon, 13:30 lines for LTST
        if extract.y_choice == 'LTST':
            plt.axhline(37800, dashes=(10, 10), color='gray', zorder=0)
            plt.axhline(43200, dashes=(10, 10), color='gray', zorder=0)
            plt.axhline(48600, dashes=(10, 10), color='gray', zorder=0)
        if extract.x_choice == 'LTST':
            plt.axvline(37800, dashes=(10, 10), color='gray', zorder=0)
            plt.axvline(43200, dashes=(10, 10), color='gray', zorder=0)
            plt.axvline(48600, dashes=(10, 10), color='gray', zorder=0)

        #This part executes the plot
        canvas.draw()

        def on_closing(): # Quits the second window
            t.destroy()
        t.protocol("WM_DELETE_WINDOW", on_closing)

        x_values = np.array(x_values)
        y_values = np.array(y_values)
        x_lim = plt.get_xlim()
        y_lim = plt.get_ylim()



app = Application(master=root)
app.mainloop()