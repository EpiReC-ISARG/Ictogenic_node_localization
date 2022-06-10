# Ictogenic nodes localization

### Objective
Epilepsy surgery fails in >30% of patients with focal cortical dysplasia (FCD). The seizure persistence after surgery can be attributed to the inability 
to precisely localize the tissue with an endogenous potential to generate seizures. In this study, we aimed to identify the critical components of the 
epileptic network that were actively involved in seizure genesis.

### Methods
The directed transfer function was applied to intracranial EEG recordings and the effective connectivity was determined with a high temporal and frequency 
resolution. Pre-ictal network properties were compared with ictal epochs to identify regions actively generating ictal activity and discriminate them from 
the areas of propagation.

### Results
Analysis of 276 seizures from 30 patients revealed the existence of a seizure-related network reconfiguration in the gamma-band (25-170 Hz; p<0.005) 
– ictogenic nodes. Unlike seizure onset zone, resecting the majority of ictogenic nodes correlated with favorable outcomes (p<0.012).

### Conclusion
The prerequisite to successful epilepsy surgery is the accurate identification of brain areas from which seizures arise. We show that in FCD-related 
epilepsy, gamma-band network markers can reliably identify and distinguish ictogenic areas in macroelectrode recordings, improve intracranial EEG 
interpretation and better delineate the epileptogenic zone.

### Significance
Ictogenic nodes localize the critical parts of the epileptogenic tissue and increase the diagnostic yield of intracranial evaluation.

## Publication
Janca R, Jahodova A, Hlinka J, Jezdik P, Svobodova L, Kudr M, Kalina A, Marusic P, Krsek P, Jiruska P (2021) Ictal gamma-band interactions localize 
ictogenic nodes of the epileptic network in focal cortical dysplasia. Clin. Neurophysiol. doi: [10.1016/j.clinph.2021.04.016](https://doi.org/10.1016/j.clinph.2021.04.016)

## Test data
Files are available on [Google Drive](https://drive.google.com/drive/folders/12WbOmEdkiEHouKzqIJ-JSLEcC886w9C5?usp=sharing).

## Tips for success

The main m-file is `avrDTF_ICN_SOZ_v01.m`, and a data dictionary can be found there. Running this file leads to calling the eight other files in the same directory.

SOZ detection (ictogenic nodes), IED detection, and sub-regional masking each only require iEEG (data matrix, time vector and sampling frequency, label – names of channels) as a minimum. Specifically:

- Consider starting with 4 minutes of pre-ictal data plus up to 1 minute of data following seizure onset.
- For IED detection, a minimum of 3 hours awake plus 3 hours asleep is recommended. Even better is a continuous segment of 24- to 48 hours.

PET postprocessing requires good quality isometric T1 (<1 mm) and FDG-PET (with CT/MR).

Guide to acronyms seen in figures:

- `sot` (optional): The seizure onset time, in datenum format. It's used only for plotting of dotted red line to figure(1) and (2) for comparison with detection (dotted cyan line).
- `soz` (optional): A vector of channels marked by a clinician, also only for plotting and comparison with algorithm's result.
- `ICNs`: The result. This represents the SOZ, i.e. the channels actively generating seizure between `sot` and +6.5s. 

## Contact
Questions and error reports sent to:

**Radek Janča**: jancarad [at] fel.cvut.cz
