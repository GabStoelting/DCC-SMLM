import h5py
import pandas as pd
import numpy as np
from scipy.stats import norm, ks_2samp
from sklearn.cluster import DBSCAN, OPTICS
from scipy.spatial import KDTree
from scipy.optimize import curve_fit
from scipy.special import binom

class DCCProteinOfInterest:
    def __init__(self, poi_filename, value_column):
        self.poi = pd.read_csv(poi_filename)
        self.value_column = value_column

        self.bs_mean = 0
        self.bs_lower_ci = []
        self.bs_upper_ci = []

        self.ks = []
        self.ks_stat = []

    def COM(self, x, E):
        return [1 - ((np.sum((x - np.mean(x)) ** 2)) / (np.sum((x - e) ** 2))) for e in E]

    def com_bootstrap(self, n_samples, reference, reference_bootstrap=True, save_result=True, max_n=4, only_p=False):

        return_df = np.array([])

        for i in range(0, n_samples):

            if reference_bootstrap:
                reference.reference_bootstrap(1, save_result=True, only_p=only_p)

                sample_resample = self.poi.sample(frac=1, replace=True)  # Resample from the POI samples

            # Calculate the
            bs_sample = self.COM(sample_resample[self.value_column].values, reference.E(n=max_n))

            return_df = np.append(return_df, bs_sample)  # Add the result to the return DataFrame

        return_df = return_df.reshape(-1, len(bs_sample))

        if save_result:
            self.bs_mean = return_df.mean(axis=0)
            self.bs_lower_ci = np.percentile(return_df, 2.5, axis=0)
            self.bs_upper_ci = np.percentile(return_df, 97.5, axis=0)

        return return_df

    def KS(self, reference, save_results=True):
        reference_unique = reference.reference[reference.oligo_column].unique()
        reference_values = [reference.reference[reference.reference[reference.oligo_column] == n]
                            [reference.value_column].values for n in reference_unique]
        ks_result = [ks_2samp(self.poi[self.value_column].values, r)[1] for r in reference_values]
        ks_stats = [ks_2samp(self.poi[self.value_column].values, r)[0] for r in reference_values]
        if save_results:
            self.ks = ks_result
            self.ks_stat = ks_stats

        return ks_result

class DCCReferenceProteins:
    def __init__(self, reference_filename, value_column, oligo_column):
        self.reference = pd.read_csv(reference_filename)
        self.value_column = value_column
        self.oligo_column = oligo_column
        self.reference_mean = self.reference.groupby(oligo_column).mean()[value_column].to_list()
        self.reference_std = self.reference.groupby(oligo_column).std()[value_column].to_list()

        self.p = 0
        self.m = 0
        self.p_ci = (0,0)
        self.m_ci = (0,0)

    def _calc_p(self, n, p):
        """ This returns the colocalization probability for only p (no m!) equivalent to Eq. 1"""
        return (1 - (1-p)**n)

    def _calc_p_m(self, n, m, p):
        """

        :param n: oligomeric state
        :param m: the probability of finding an intact marker
        :param p:  the probability of detecting an indicator protein if the marker is intact
        :return:
        """
        m = 1-m # We inverted the definition of m

        # The following corresponds to Eq. 6 in our paper
        def _calc(k):
            return (1 - (1 - p) ** k) * (binom(n, k) * m ** k * (1 - m) ** (n - k)) / (1 - (1 - m) ** n)

        k_list = range(1, n + 1)
        P_list = list(map(_calc, k_list))
        P_sum = sum(P_list)

        return P_list, P_sum

    def _calc_p_m_curve_fit(self, n_list, m, p):
        """ This function is just a wrapper to use calc_p with curve_fit """
        ret_s = []
        for n in n_list:
            l, s = self._calc_p_m(int(n), m, p)
            ret_s = np.append(ret_s, s)

        return ret_s

    def reference_bootstrap(self, n_samples, save_result=True, only_p=False):
        """"
        # Run a bootstrap resampling over the reference calibration values
        #
        # :param n_samples: Number of resamples
        # :param save_result: Save the result to the corresponding variables of this class (default=True)
        # :param only_p: Only fit p but not m (default=False)
        """

        return_df = np.array([])

        for i in range(0, n_samples):
            rs = self.reference.sample(frac=1, replace=True)  # Resample from the reference distribution
            # Calculate mean values per oligomeric state
            rm = rs.groupby(self.oligo_column).mean()[self.value_column].to_list()
            if not only_p:
                popt, pcov = curve_fit(self._calc_p_m_curve_fit, [1, 2, 3, 4], rm, p0=(0.76, 0.17))  # Fit Eq. 6
            elif only_p:
                popt, pcov = curve_fit(self._calc_p, [1, 2, 3, 4], rm, p0=(0.17))  # Fit Eq. 1 (only p!)
            return_df = np.append(return_df, popt)  # Add the result to the return DataFrame

        return_df = return_df.reshape(-1, popt.shape[0])

        if save_result:
            if not only_p:
                self.m = return_df.mean(axis=0)[0]
                self.m_ci = (np.percentile(return_df, 2.5, axis=0)[0], np.percentile(return_df, 97.5, axis=0)[0])

                self.p = return_df.mean(axis=0)[1]
                self.p_ci = (np.percentile(return_df, 2.5, axis=0)[1], np.percentile(return_df, 97.5, axis=0)[1])
            elif only_p:
                self.p = return_df.mean(axis=0)
                self.p_ci = (np.percentile(return_df, 2.5, axis=0), np.percentile(return_df, 97.5, axis=0))

        return return_df

    def E(self, n=4, only_p=False):
        """ This function returns the expected colocalization values based on p (and m if only_p == True)"""
        if self.p > 0:
            if not only_p:
                return [self._calc_p_m(k, self.m, self.p)[1] for k in range(1, n+1)]
            elif only_p:
                return [self._calc_p(k, self.p) for k in range(1, n + 1)]
        else:
            raise ValueError("Please fit and save the reference proteins first")

class Channel:
    def __init__(self, raw_data, x_column="x", y_column="y", intensity_column="intensity", frame_column="frame"):
        """
        Setup all necessary variables for each channel

        :param raw_data: Contains the raw data
        :param x_column: defines the name of the column containing the x-coordinate
        :param y_column: defines the name of the column containing the x-coordinate
        :param intensity_column: defines the name of the column containing the intensity value
        :param frame_column: defines the name of the column containing the frame number
        """
        self.x_column = x_column
        self.y_column = y_column
        self.intensity_column = intensity_column
        self.frame_column = frame_column
        self.raw_data = raw_data
        self.drift = pd.DataFrame([])
        self.drift_corrected = False
        # The following are the parameters for chromatic aberration correction
        self.Kx = 1
        self.Ky = 1
        self.x0 = 0
        self.y0 = 0


    def _get_roi_subselection(self, roi):
        """
        Returns a boolean list indicating the entries in self.raw_data that are within the defined ROI

        :param roi: rectangular ROI defined by the following format: (x_0, y_0, width, height)
        :return: boolean list
        """
        roi_selection = ((self.raw_data[self.x_column]>=roi[0]) & (self.raw_data[self.x_column]<(roi[0]+roi[2])) &
                        (self.raw_data[self.y_column]>=roi[1]) & (self.raw_data[self.y_column]<(roi[1]+roi[3])))

        return roi_selection

    def _linear_func(self, x, a, b):
        """ Just a linear function"""
        return a * x + b

    def _generate_prototypical_spot(self, sigma):
        """
        This creates a template of a single localization spot for the reconstruction based on a 2D gaussian

        :param sigma: Sigma value used for the 2D gaussian
        :return: numpy array with the single spot
        """
        # creates a prototype 2d gaussian spot from the product of two 1d normal distributions
        if sigma == 0:  # if sigma is 0, the spot will be a singular 1 at least
            return 1
        else:
            # create an empty array of width and height equal to 5 * sigma
            prototypicalSpot = np.zeros((sigma*5, sigma*5))

            xDist = norm(loc = sigma*2.5, scale=sigma)
            xRange = np.arange(0,len(prototypicalSpot), 1)

            yDist = norm(loc = sigma*2.5, scale=sigma)
            yRange = np.arange(0,len(prototypicalSpot), 1)

            x2d = np.repeat(xDist.pdf(xRange), len(prototypicalSpot))
            x2d = x2d.reshape((len(prototypicalSpot),len(prototypicalSpot)))

            y2d = np.repeat(yDist.pdf(yRange), len(prototypicalSpot))
            y2d = y2d.reshape((len(prototypicalSpot),len(prototypicalSpot)))

            return (np.rot90(y2d)*x2d)

    def _generate_spot(self, image, x, y, prototypicalSpot):
        """
        Create a localization spot at a specific coordinate for a reconstruction

        :param image: array for the reconstruction
        :param x: x coordinate
        :param y: y coordinate
        :param prototypicalSpot: array containing a prototypical, template spot
        :return: the full image
        """

        # Pastes a copy of a prototypical spot to the x,y coordinates in the image array
        if(isinstance(prototypicalSpot, np.ndarray)): # if the prototypical spot is an array, put that at x,y
            spotSize=len(prototypicalSpot)
            arrayShift = spotSize % 2 # one needs to add 1 to the array subselection for uneven, 0 for even numbers
            image[(y-int(spotSize/2)):(y+int(spotSize/2))+arrayShift,
            (x-int(spotSize/2)):(x+int(spotSize/2))+arrayShift]=image[(y-int(spotSize/2)):(y+int(spotSize/2))+arrayShift,
            (x-int(spotSize/2)):(x+int(spotSize/2))+arrayShift]+prototypicalSpot
        else: # if it is just a single value set x,y directly to that value
            image[x,y]=image[x,y]+prototypicalSpot
        return image

    def generate_reconstruction(self, raw_data=None, reconstructionType="binary",
                               imWidth=40960, imHeight=40960, spotSigma=16, nmPerPx=16, threshold=30, bits=8,
                               cluster_column=False, cluster_threshold=0):
        """
        Creates a reconstruction of the localizations

        :param raw_data: Raw table with the localizations
        :param reconstructionType: Currently, "binary" (each spot has a maximum o 255) or "normalized" (each spot
        is normalized to a maximum of 1 and overlapping localizations result in a summation) are possible
        (default: binary)
        :param imWidth: Width of the full image in the same unit as the raw data (typically nm)
        :param imHeight: Height of the full image in the same unit as the raw data (typically nm)
        :param spotSigma: Sigma used for the 2D gaussian of each spot (default: 16)
        :param nmPerPx: Definition of the nm per pixel of the final image (default: 16)
        :param threshold: Only used when reconstructionType is "normalized" (default: 30)
        :param bits: Bit depth of the final image (default: 8)
        :param cluster_column: Define the name of the cluster ID column in raw_data
        :param cluster_threshold: Threshold used for the cluster ID (DBSCAN, for example, assigns an ID of -1 for
        unclustered localizations)
        :return: Array with the completed reconstruction
        """

        # If raw_data is not defined, take the full raw_data of the class
        if isinstance(raw_data, pd.DataFrame):
            Data = raw_data[(raw_data[self.x_column] > 0) & (raw_data[self.y_column] > 0)].copy()
        else:
            Data = self.raw_data[(self.raw_data[self.x_column] > 0) & (self.raw_data[self.y_column] > 0)].copy()

        # If cluster_column is defined, restrict the reconstruction to those localizations
        # that have a cluster ID number above the cluster_threshold (>=0 means it's clustered)
        if cluster_column:
            Data = Data[Data[cluster_column] >= cluster_threshold]

        sigma = int(spotSigma / nmPerPx)  # round spotSigma (nm) to a number of pixels
        padsize = 2 * 5 * sigma  # the size of the padding should be a maximum of a full spot width in pixels
        imWidthPx = int(imWidth / nmPerPx)
        imHeightPx = int(imHeight / nmPerPx)
        imarray = np.zeros((imWidthPx, imHeightPx))  # create an empty array of defined width and height'
        # add padding so that spots may lay directly on the border of the image
        imarray = np.pad(imarray, padsize, "edge")
        maxValue = 2 ** bits - 1  # the maximum value for the the defined number of bits per pixel

        spot = self._generate_prototypical_spot(sigma)  # create a prototypical spot

        if (reconstructionType == "binary"):  # check if a binary reconstruction is chosen
            # normalize the prototypical spot to 1 and multiply by 255 for a binary reconstruction
            spot = (spot / np.max(spot)) * maxValue
        elif (reconstructionType == "normalized"):  # else check if a normalized reconstruction is chosen
            spot = (spot / np.max(spot))  # normalize the prototypical spot to 1
        else:  # if something else has been chosen as reconstruction type
            return False  # stop everything

        for row in Data[[self.x_column, self.y_column]].values:  # iterate through all rows in the dataset
            imarray = self._generate_spot(imarray, int(row[0] / nmPerPx) + padsize, int(row[1] / nmPerPx) + padsize,
                                   spot)  # put spots at all locations
        if (reconstructionType == "normalized"):  # check if a normalized reconstruction is chosen
            # divide everything by the threshold and multiply by the maxValue
            imarray = (imarray / threshold) * maxValue

        if (padsize >= 1):  # padding needs to be removed only if any exists
            # return the array while removeing the padding and normalizing everything to maxValue
            return np.clip(imarray[padsize:-padsize, padsize:-padsize], None, maxValue)
        else:  # otherwise just return the clipped array
            return np.clip(imarray, None, maxValue)

    def find_clusters(self, intensity_threshold=0, min_samples=100, eps=0, save_column=False):
        """
        Find clusters using DBSCAN or OPTICS

        :param intensity_threshold: Define a minimum intensity for inclusion
        :param min_samples: Define the minimum number of localizations per cluster
        :param eps: Epsilon parameter used by DBSCAN. If defined, DBSCAN will be used, if set to 0,
        OPTICS will be used instead
        :param save_column: Defines the column name which will afterwards contain the cluster ID numbers
        :return:
        """
        # ToDo: Implement search based on coordinates/ROIs
        if intensity_threshold > 0:
            marker_coordinates = self.raw_data[self.raw_data[self.intensity_column]>intensity_threshold][[
                                                                                            self.frame_column,
                                                                                             self.x_column,
                                                                                             self.y_column]].copy()
            # If no epsilon parameter is given, use OPTICS
            if eps == 0:
                # Run OPTICS on all x/y localizations with an intensity above the threshold
                bead_clusters = OPTICS(min_samples=min_samples).fit(marker_coordinates[[self.x_column,
                                                                                        self.y_column]].values)
            else:
                # Run DBSCAN on all x/y localizations with an intensity above the threshold
                bead_clusters = DBSCAN(min_samples=min_samples, eps=eps).fit(marker_coordinates[[self.x_column,
                                                                                        self.y_column]].values)

            if save_column != False:
                self.raw_data[save_column] = -1
                self.raw_data.loc[self.raw_data[self.intensity_column]>intensity_threshold,
                                  save_column] = bead_clusters.labels_

                return
            marker_coordinates["marker_id"] = bead_clusters.labels_

            return marker_coordinates

    def get_mean_clusters(self, selection_column, threshold=0, rois=None):
        """
        Get the mean values for each cluster

        :param selection_column: Defines the cluster ID column of raw_data
        :param threshold: Defines the cluster ID threshold (e.g. larger than -1 if only clustered localizations should
        be considered)
        :param rois: list of ROIs (each being a set of (x_0, y_0, width, height))
        :return: DataFrame with the mean values for each cluster
        """
        threshold_selection = self.raw_data[selection_column] >= threshold
        roi_selection = None

        if rois != None:
            for i, roi in enumerate(rois):
                # ToDo: Check that the roi coordinates are valid
                # Currently rois must be in the following, square, form: (x, y, x_width, y_width)
                if i == 0:
                    roi_selection = self._get_roi_subselection(roi)
                else:
                    roi_selection = roi_selection ^ self._get_roi_subselection(roi)

            return self.raw_data[(roi_selection) & (threshold_selection)].groupby(selection_column).mean()
        else:
            return self.raw_data[threshold_selection].groupby(selection_column).mean()


    def extract_drift(self, df, id_column, id_list=0):
        """
        Extract drift from identified clusters

        :param df: DataFrame containing the localizations to use for drift information (e.g. fiducial markers)
        :param id_column: Name of the column containing the cluster IDs for each marker
        :param id_list: List of cluster IDs to use for drift extraction
        :return: Returns a DataFrame containing the drift per frame
        """
        # if no id_list is given, create a list based on all available IDs in the id_column
        if id_list == 0:
            id_list = df[id_column].unique()

        # Create a pandas DataFrame for return
        return_df = pd.DataFrame([])

        # Iterate over all IDs
        for id in id_list:
            cluster = df[df[id_column] == id]

            # This list contains all frame numbers present in df
            frames_present = cluster.sort_values(self.frame_column).frame
            # This is a list of all frames from the first to the last frame of df
            first_frame = int(cluster.frame.min())
            last_frame = int(cluster.frame.max())
            all_frames = range(first_frame, last_frame)

            # Append empty x and y coordinates for the missing frames
            for missing_frame in set(all_frames).difference(frames_present):
                cluster = cluster.append({self.frame_column: missing_frame,
                                          self.x_column: np.nan,
                                          self.y_column: np.nan,
                                          "marker_id": id},
                                         ignore_index=True)

            # Interpolate missing coordinates using linear interpolation
            cluster = cluster.sort_values(self.frame_column).interpolate()

            # Get first x and y coordinate
            first_x = cluster[cluster.frame == first_frame][self.x_column].values[0]
            first_y = cluster[cluster.frame == first_frame][self.y_column].values[0]

            # Calculate the movement in x relative to the first_frame
            cluster["delta_x"] = cluster[self.x_column] - first_x
            # Calculate the movement in y relative to the first_frame
            cluster["delta_y"] = cluster[self.y_column] - first_y

            # Append interpolated values to the return DataFrame
            return_df = pd.concat([return_df, cluster])

        self.drift = return_df.groupby(self.frame_column).mean().copy()
        return return_df

    def correct_drift(self):
        """
        Correct drift for all localizations
        :return:
        """
        # ToDo: Check if drift is defined


        def subtractDrift(row):
            frame = int(row.frame)
            row["raw_x"] = row[self.x_column]
            row["raw_y"] = row[self.y_column]
            row[self.x_column] = row[self.x_column] - self.drift.loc[frame].delta_x
            row[self.y_column] = row[self.y_column] - self.drift.loc[frame].delta_y
            return row
        if not self.drift_corrected:
            self.raw_data = self.raw_data.apply(subtractDrift, axis=1)
            self.drift_corrected = True
        else:
            raise ValueError("This channel has already been drift subtracted. You need to load the channel \
            again to correct drift.")

    def set_chromatic_aberration_params(self, Kx, Ky, x0, y0):
        """ Define CA parameters for this channel """
        self.Kx=Kx
        self.Ky=Ky
        self.x0=x0
        self.y0=y0

    def correct_chromatic_aberration(self):
        """
        Correct chromatic aberration based on stored values

        :return:
        """
        # ToDo: Check if CA has been defined
        # ToDo: Check if CA has already been subtracted and refrain from correcting again
        def correctCA(row):
            row[self.x_column] = row[self.x_column]-self._linear_func(row[self.x_column], self.Kx, self.x0)
            row[self.y_column] = row[self.y_column]-self._linear_func(row[self.y_column], self.Ky, self.y0)
            return row

        self.raw_data = self.raw_data.apply(correctCA, axis=1)

class DCCSMLM():
    """
    This is the main class for the use of DCCSMLM. It contains all the channels and allows to perform
    all necessary operations.
    """
    def __init__(self):
        """This just initiates all necessary variables"""
        self.channel = {}

    def load_channel(self, filename, **kwargs):
        """
        Loads a channel from a given file.

        :param filename: String with the filename to load
        :type filename: str
        :param kwargs: Contains further parameters needed for specific file formats
        :return: May return a list of channels if a SciH5 file is loaded without specifying any channel to load
        """

        # If the filename ends with ".scih5", this is a SciH5 file
        # (SNSMIL + Post-Processing through Johnny Hendriks software), open
        # file using h5py and extract the channel using load_scih5()
        if filename.split(".")[-1].lower() == "scih5":
            # To load a SciH5 file, a channel within the file must be defined
            raw_file = h5py.File(filename, "r")  # Load SciH5 (HDF5) file

            if "channel" in kwargs:
                self.load_scih5(raw_file, **kwargs)
            else:
                print("Please define a channel to load. Available channels are:")
                return list(raw_file["SciDaP"]["_Microsocpy_"]["_SML_"]["_Stg1_Extract_"])

        # If this is a Matlab data file (.mat), load via the load_mat function
        # Currently, this is only meant to load results from SMAP and, as .mat files are just HDF5 files,
        # will just return the HDF file.
        elif filename.split(".")[-1].lower() == "mat":
            raw_file = h5py.File(filename, "r") # Load SMAP (.mat) file
            self.load_smap(raw_file, **kwargs)

    def _check_channels(self, channel):
        """ Return True if a channel with the name exists"""
        if channel in self.channel.keys():
            return True
        else:
            return False

    def _load_ca_file(self, ca_file, channel_id):
        """
        Loads parameters for the correction of chromatic abberation from a CSV file and applies these on the
        localizations contained in the channel

        :param ca_file: Filename of the CSV file containing information about the linear CA
        :type ca_file: str
        :param channel_id: ID of the channel to correct
        :type channel_id: str
        :return:
        """
        if not self._check_channels(channel_id):
            raise KeyError(f"The channel {channel_id} doesn't exist.")


        ca_params = pd.read_csv(ca_file)  # Load CSV file
        # Fix CA correction parameters for the channel
        self.channel[f"{channel_id}"].set_chromatic_aberration_params(
            ca_params.Kx, ca_params.Ky, ca_params.x0, ca_params.y0)

        # Apply CA on the channel
        self.channel[f"{channel_id}"].correct_chromatic_aberration()

    def load_scih5(self, raw_file, overwrite=False, **kwargs):
        """
        Loads a channel from a SciH5 file.

        :param raw_file: h5py File object
        :param overwrite: Set this to True to overwrite an existing channel
        :param kwargs: Further parameters. Currently, "channel_id" (ID for the channel to be loaded) and
        "ca_file" (CSV file containing paramers for chromatic aberration correction) may be set.

        :return:
        """

        channel_id = 0

        # If the channel ID is defined, use that. If not, it'll just be "0"
        if "channel_id" in kwargs:
            channel_id = kwargs["channel_id"]
            if not overwrite:
                if self._check_channels(channel_id):
                   raise KeyError(f"The channel {channel_id} already exists. Set overwrite=True if necessary")

        # Extract the data for the chosen channel and convert to Pandas DataFrame
        self.channel[f"{channel_id}"] = Channel(pd.DataFrame(raw_file["SciDaP"]["_Microsocpy_"]["_SML_"]
                                                ["_Stg1_Extract_"][kwargs["channel"]][:]))

        # Check if a chromatic aberration file is specified
        if "ca_file" in kwargs:
            self._load_ca_file(kwargs["ca_file"], channel_id)  # Load CA params and correct localizations accordingly

    def load_smap(self, raw_file, overwrite=False, **kwargs):
        """
        Loads a channel from a SMAP file (using the .mat Matlab file format containing a HDF5 structure)

        :param raw_file: h5py File object
        :param kwargs: Further parameters. Currently, "channel_id" (ID for the channel to be loaded) and
        "ca_file" (CSV file containing paramers for chromatic aberration correction) may be set.
        :return:
        """
        channel_id = 0

        # If the channel ID is defined, use that. If not, it'll just be "0"
        if "channel_id" in kwargs:
            channel_id = kwargs["channel_id"]
            if not overwrite:
                if self._check_channels(channel_id):
                   raise KeyError(f"The channel {channel_id} already exists. Set overwrite=True if necessary")

        # Get the columns for the data from the saveloc/loc storage
        smap_columns = list(raw_file["saveloc/loc"].keys())

        # Create a dictionary with the columns and data
        smap_data = {column: raw_file["saveloc/loc"][column][0][:] for column in smap_columns}

        # Create a new Channel
        self.channel[f"{channel_id}"] = Channel(pd.DataFrame(smap_data),
                                                x_column="xnm", y_column="ynm", intensity_column="phot")

        # Check if a chromatic aberration file is specified
        if "ca_file" in kwargs:
            self._load_ca_file(kwargs["ca_file"], channel_id)  # Load CA params and correct localizations accordingly

    def _linear_func(self, x, a, b):
        """ Just a linear function """
        return a * x + b

    def get_colocalization_ratio(self, channel_one, channel_two, cluster_col, distance_cutoff=250, rois=None):
        """
        Calculates the colocalization ratio.

        :param channel_one: ID of the first channel (marker)
        :param channel_two: ID of the second channel (indicator)
        :param cluster_col: Name of the DataFrame column that contains the cluster IDs
        :param distance_cutoff: Maximum distance to search for colocalization (default: 250)
        :param rois: List of rectangular ROIs. Each ROI is a set of the following format:
         (x_0, y_0, width, height)
        :return: (float) colocalization ratio
        """

        if not self._check_channels(channel_one):
            raise KeyError(f"Channel {channel_one} doesn't exist.")
        if not self._check_channels(channel_two):
            raise KeyError(f"Channel {channel_two} doesn't exist.")

        # The number of colocalized clusters corresponds to the number of instances with a closest cluster within
        # the distance_cutoff
        n_colocalized_clusters = self.get_closest_cluster(channel_one, channel_two, cluster_col,
                                                          distance_cutoff=distance_cutoff, rois=rois).shape[0]

        # Subtract one because all the ungrouped localizations will become group "-1" which is counted but
        # actually isn't a real cluster
        n_clusters_channel_one = self.channel[channel_one].get_mean_clusters(cluster_col, rois=rois).shape[0]-1

        # Return the ratio
        return n_colocalized_clusters/n_clusters_channel_one

    def get_closest_cluster(self, channel_one, channel_two, cluster_col, distance_cutoff=250, rois=None):
        """
        Returns a DataFrame with the closest clusters for each cluster in channel one

        :param channel_one: ID of the first channel (marker)
        :param channel_two: ID of the second channel (indicator)
        :param cluster_col: Name of the DataFrame column that contains the cluster IDs
        :param distance_cutoff: Maximum distance to search for colocalization (default: 250)
        :param rois: List of rectangular ROIs. Each ROI is a set of the following format:
         (x_0, y_0, width, height)
        :return: (Pandas DataFrame) with the closest clusters.
        """

        if not self._check_channels(channel_one):
            raise KeyError(f"Channel {channel_one} doesn't exist.")
        if not self._check_channels(channel_two):
            raise KeyError(f"Channel {channel_two} doesn't exist.")


        # Load the name definitions of x and y coordinates from the first channel
        x_col_one = self.channel[channel_one].x_column
        y_col_one = self.channel[channel_one].y_column

        # Load the name definitions of x and y coordinates from the second channel
        x_col_two = self.channel[channel_one].x_column
        y_col_two = self.channel[channel_one].y_column

        # Obtain the center (mean) coordinates for each cluster in both channels
        c1_data = self.channel[channel_one].get_mean_clusters(cluster_col,
                                                              threshold=0, rois=rois)[[x_col_one, y_col_one]].values
        c2_data = self.channel[channel_two].get_mean_clusters(cluster_col,
                                                              threshold=0, rois=rois)[[x_col_two, y_col_two]].values
        # Build KD-Tree for faster search
        c2_KDTree = KDTree(c2_data)

        return_list = []
        return_columns = ["x_distance", x_col_one, "y_distance", y_col_one]

        # Iterate over all clusters in channel one and find the closest in channel two
        for i, cluster_mean in enumerate(c1_data):
            distance, index = c2_KDTree.query(cluster_mean, distance_upper_bound=distance_cutoff, k=1)

            if distance < distance_cutoff:
                cluster = [c2_data[index][0] - cluster_mean[0],  # x_distance to nearest
                           cluster_mean[0],  # x of nearest
                           c2_data[index][1] - cluster_mean[1],  # y_distance to nearest
                           cluster_mean[1]] # y of nearest

                return_list.append(cluster)
        return pd.DataFrame(return_list, columns=return_columns)

    def get_nearby_cluster_distances(self, channel_one, channel_two, cluster_col, distance_cutoff=250, k=2, rois=None):
        """
        Returns a DataFrame with the distances to the k nearest clusters.

        :param channel_one: ID of the first channel (marker)
        :param channel_two: ID of the second channel (indicator)
        :param cluster_col: Name of the DataFrame column that contains the cluster IDs
        :param distance_cutoff: Maximum distance to search for colocalization (default: 250)
        :param k: The maximum number of closest clusters to return
        :param rois:List of rectangular ROIs. Each ROI is a set of the following format:
         (x_0, y_0, width, height)
        :return: DataFrame with the distances to the k nearest clusters
        """
        if not self._check_channels(channel_one):
            raise KeyError(f"Channel {channel_one} doesn't exist.")
        if not self._check_channels(channel_two):
            raise KeyError(f"Channel {channel_two} doesn't exist.")

        # Load the name definitions of x and y coordinates from the first channel
        x_col_one = self.channel[channel_one].x_column
        y_col_one = self.channel[channel_one].y_column

        # Load the name definitions of x and y coordinates from the second channel
        x_col_two = self.channel[channel_one].x_column
        y_col_two = self.channel[channel_one].y_column

        # Get mean positions of clusters. Set threshold to 0 to avoid including non-clustered locs (-1)
        c1_data = self.channel[channel_one].get_mean_clusters(cluster_col, threshold=0,
                                                              rois=rois)[[x_col_one, y_col_one]].values
        c2_data = self.channel[channel_two].get_mean_clusters(cluster_col, threshold=0,
                                                              rois=rois)[[x_col_two, y_col_two]].values
        # Build KD-Tree for faster search
        c2_KDTree = KDTree(c2_data)

        return_list = []
        # Create a list of k columns that will be the header for the returned DataFrame
        return_columns = [f"distance_{dist_col}" for dist_col in range(0,k)]

        for i, cluster_mean in enumerate(c1_data):
            # Search for the k nearest clusters within the distance_cutoff
            distance, index = c2_KDTree.query(cluster_mean, distance_upper_bound=distance_cutoff, k=k)
            cluster = [dist for dist in distance]
            return_list.append(cluster)

        return pd.DataFrame(return_list, columns=return_columns)

    def determine_chromatic_aberration(self, channel_one, channel_two, cluster_col, distance_cutoff=250,
                                       save_ca_params=True):
        """
        Calculate the chromatic aberration based on identified clusters in two channels.

        :param channel_one: ID of channel one
        :param channel_two: ID of channel two
        :param distance_cutoff: Maximum distance to consider for colocalization
        :param save_ca_params: Save the CA parameters to the 2nd channel
        :param return_df: Give a DataFrame to store the individual coordinates in, if desired
        :return: Returns the parameters from the fit in x and y direction (K_x, K_y, x_0, y_0)
        """

        if not self._check_channels(channel_one):
            raise KeyError(f"Channel {channel_one} doesn't exist.")
        if not self._check_channels(channel_two):
            raise KeyError(f"Channel {channel_two} doesn't exist.")

        # Load the name definitions of x and y coordinates from the first channel
        x_col_one = self.channel[channel_one].x_column
        y_col_one = self.channel[channel_one].y_column

        # Determine closest clusters (as surrogate for chromatic aberration (CA)) of clusters in channel_two
        # relative to channel_one

        ca = self.get_closest_cluster(channel_one, channel_two, cluster_col, distance_cutoff=distance_cutoff)

        # Fit chromatic aberration parameters along the x-axis
        popt_x, pcov_x = curve_fit(self._linear_func, ca[x_col_one].values, ca["x_distance"].values)
        # Fit chromatic aberration parameters along the y-axis
        popt_y, pcov_y = curve_fit(self._linear_func, ca[y_col_one].values, ca["y_distance"].values)

        if save_ca_params:
            self.channel[channel_two].set_chromatic_aberration_params(popt_x[0], popt_y[0], popt_x[1], popt_y[1])

        return pd.DataFrame(popt_x[0], popt_y[0], popt_x[1], popt_y[1]), ca

    def pair_correlation(self, channel_one, channel_two, distance_cutoff=400, k=50, delta_r=20):
        """
        Calculates the pair correlation function between two channels. The implementation is based on
        http://www.physics.emory.edu/faculty/weeks//idl/gofr2.html

        :param channel_one: ID of channel one
        :param channel_two: ID of channel two
        :param distance_cutoff: Maximum distance to consider
        :param k: Maximum number of closest clusters to consider
        :param delta_r: Size of the
        :return:
        """
        # Based on http://www.physics.emory.edu/faculty/weeks//idl/gofr2.html

        # ToDo: Make (square) ROI selection possible
        # ToDo: Check if channels exist!
        # ToDo: Deal with edge cases
        # determine cluster distances. ToDo: Check if clusters have been detected before!
        cluster_distances = self.get_nearby_cluster_distances(channel_one,
                                                                  channel_two,
                                                                  cluster_col="clusters",
                                                                  distance_cutoff=distance_cutoff,
                                                                  k=k,
                                                              rois=rois)
        # Check if all values for the k-th distances are np.inf
        # ToDo: If not, raise error saying there are too many clusters for the given k and distance_cutoff
        if (cluster_distances[f"distance_{k-1}"].unique() == np.inf).all():
            # The number of entries in cluster_distances corresponds to the number of marker proteins (M)
            N_M = cluster_distances.shape[0]
            # The number of clusters in the second channel corresponds to the number of indicator proteins (F)
            N_F = self.channel[channel_two].raw_data["clusters"].unique().max()

            cluster_distances = cluster_distances.values.reshape(-1) # Reshape into a 1D numpy array

            # Make a histogram of the (non inf) distances
            # H[0] contains the counts, H[1] the bin borders
            H = np.histogram(cluster_distances[cluster_distances!=np.inf], bins=range(0, distance_cutoff, delta_r))
            norm_H = H[0] / N_M  # Normalize the counts by the number of M
            A = 2 * np.pi * H[1][1:] * delta_r # Calculate the areas of the rings (of delta_r width) around M
            # ToDo: Make the area of the field of view variable
            density = N_F / (40960 * 40960) # Calculate the density of F across the whole field of view
            # return the bins and the pair correlation
            return H[1][1:], ((norm_H / A) / density)




