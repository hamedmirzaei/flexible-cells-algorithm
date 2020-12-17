# Flexible Cells Algorithm
Consider a set of geo-tagged events, such as the location of COVID-19 cases or car accidents. One important thing about these events is that they tend to cluster into "centers"; by center, we mean locations with high densities of events. For example, in the case of COVID-19, in which country, city, or even smaller areas does a problem seem to be occurring?
In this project we try to find centers of a set of geo-tagged events (e.g. COVID-19 related Tweets) as well as their spread and dispersion.

For all the code in this repository, parameters can be changed in the code itself.

## FCA
This folder contains our algorithm's implementation on finding multiple centers. It is much like a breadth-first search algorithm plus early termination of bigger cells.
It has three main phases
1) Searching Phase: finding candidate cells using BFS and early termination
2) Checking Phase: removing false positives from results of searching phase using brute-force checking of all possible cases
3) Parameter Tuning Phase: using Backstrom's probabilistic framework

## SDG
This folder contains our SyntheticDataGenerator algorithm which generates datasets. with a specific number of events using uniform distribution plus highly concentrated cells to serve as "centers".

## BSC
This folder contains implementation of Single-Center case of Backstrom's algorithm, plus a SyntheticDataGenerator which generates a dataset to show the weaknesses of Backstrom's work.

## data
Here, we store all the data used for this project. COVID_coords_all.csv contains the latitude and longitude of COVID-related Tweets from the continental USA. The rest are synthetically generated. Files that include "_b" at the end are generated in a way that shows the weakness of our baseline, Backstrom et al.'s algorithm, while those without this extension are generated with the default procedure. "Synthetic events" contain the entire output of generating the centers, while "synthetic coords" files only contain the latitudes and longitudes.

## notebooks
In this folder, there are four Jupyter Notebooks. To use these (which is not necessary to run our procedure), one should install Anaconda on their computer. "restrict_to_USA.ipynb" takes in a data set with coordinates and creates a feature "in_US", which is a boolean based on whether the coordinates from each row originate from the continental United States. "dataset_to_coordinates.ipynb" prepares the output of the previous notebook for using with our method. "plot_synthetic_map.ipynb" and "plot_true_map.ipynb" are used to visualize the predicted centers on a map with the true events.

## visualizations
This folder contains two sub-folders: COVID_data_16_centers and synthetic_data. In COVID_data_16_centers, there are html files that can be viewed in a browser that contain the results of a run of our method on the COVID dataset with 16 centers. synthetic_data contains the visualizations of centers from both our method and Backstrom's. These files are interactive.

### Contributors
Hamed Mirzaei<br/>
Sarah Davis
