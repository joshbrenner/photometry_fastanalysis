# photometry_fastanalysis

This MATLAB code is used to preprocess photometry data from mice for further analysis. The code loads the photometry data from specific mice and preprocesses it by downsampling, extracting specific signals, and performing linear fitting to subtract hemodynamic signals. The code then finds trial times based on a threshold and calculates the length of each trial. The code also creates an array to store the data for each mouse.

The code starts by specifying the path to the folder containing the photometry data and the names of the mice whose data will be processed. The onset and offset times of each experiment are also specified in the onset_cell and offset_cell variables. These variables are used to extract the data for each experiment.

The code defines a function named analyzethosemice that takes in the experiment number and the onset and offset times as arguments. The function performs the following steps:

    1. Loads the photometry data from the specified file.
    2. Extracts specific signals, such as the calcium signal, the hemodynamic signal, the signal from the linear encoder, and the signal from the photodiode.
    3. Downsamples the extracted signals to reduce noise.
    4. Performs linear fitting to subtract the hemodynamic signal from the calcium signal.
    5. Calculates the trial times based on a threshold value and finds the length of each trial.

The code then calls the analyzethosemice function for each mouse and stores the data for each mouse in the mice array.

Overall, this code is useful for preprocessing photometry data from mice for experiments with animal generated motion described in Figure 6 of Brenner et al. 2023 and extracting relevant signals for further analysis. 

This code is provided by Joshua Brenner and Sebastian Zahler.
