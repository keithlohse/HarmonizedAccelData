# R Code to Compute Sensor-Derived Upper Limb Variables from Bilateral Wrist-Worn Accelerometers

This repository includes the R code scripts used to process all upper limb accelerometry data for the project, Harmonized Upper and Lower Limb Accelerometry Data. This is a harmonized dataset of accelerometry data collected by Washington University School of Medicine investigators and collaborators. The harmonized dataset includes data from eight studies, 790 individuals with a variety of health conditions, and 2,885 recording days of accelerometry data. The dataset has been submitted to the National Institute of Child Health and Human Development (NICHD) Data and Specimen Hub (DASH) repository.  This README document will be updated with a link to the dataset, once it is available.

The purpose of sharing this code is to provide transparency with how the upper limb accelerometry data were processed and allow others to replicate these processing procedures, if they choose. In all eight studies, participants wore bilateral wrist accelerometers. During the harmonization process, all upper limb accelerometry data were processed using custom written R code scrips. All code scripts calculated the upper limb sensor variables the same way; however, variations in the number of recording days and in-lab vs. out-of-lab time for each study protocol resulted in some differences in processing procedures across studies. All upper limb sensor variables were calculated from recording times in daily life (i.e., outside the clinic/lab), with the exception of the Activity Level and Intensity of Older Adults in Skilled Nursing Rehabilitation Measured via Actigraphy in which participants were residing in a skilled nursing facility. The table below displays the code script used to process data from each study as well as details about the processing procedures.

|Study	| R Script Name |	Processing Procedures|
|-------|---------------|------------------------|
|Mapping Upper Limb Performance | MappingCohort_PMC844293.R|	One day of accelerometry data processed, per study protocol.|
|Evaluating Capacity versus Performance During Outpatient Neurorehabilitation |	CapacityVsPerformance_PMC9750113.R	| Three days of accelerometry data processed, per study protocol.|
|Does Task-Specific Training Improve Upper Limb Performance in Daily Life Poststroke?|	Dose_PMC5016233.R|	One day of accelerometry data processed, per study protocol. First 1.5 hours of recording time were removed due to in-lab time. Files were trimmed if > 24 hours after removing first 1.5 hours.|
|Activity Level and Intensity of Older Adults in Skilled Nursing Rehabilitation Measured via Actigraphy|	EMR_PMID32004240.R|	One day of accelerometry data processed, per study protocol. File trimmed if > 24 hour recording period.|
|Upper Limb Activity in Adults: Referent Values Using Accelerometry	|AdultReferent_PMC3905245.R|	One day of accelerometry data processed, per study protocol. First 2 hours of recording time were removed due to in-lab time.|
|Detection of Pediatric Upper Extremity Motor Activity and Deficits with Accelerometry|	Pediatric_PMC6487720.R	|Four days of accelerometry data were processed, per study protocol. The 30 minutes at the beginning and end of the file for each day were trimmed.|
|A Feasibility Study of Bilateral Wrist Sensors for Measuring Motor Traits in Children with Autism|	AutismPilot_PMC9974780.R	|Two, 12-hour recording days per participant were processed, per the study protocol.|
|Associations Between Coordination and Wearable Sensor Variables	|Encoding_PMID38189355.R|	Two days of accelerometry data were processed. In-lab time was removed.|

For simplicity, this repository also includes a generic code script, UpperLimbAccelerometry.R that provides an overview of each of the upper limb sensor variable calculations for a 24-hour recording period.
Executing these code scripts requires reading in four csv files: two down-sampled 1 Hz data files extracted from ActiLife 6 software (in ActiGraph activity counts, ActiGraph Inc, Pensacola, Florida) and two raw 30 Hz acceleration data files (in gravitational units), also extracted from ActiLife 6 software. Additional details regarding ActiGraph’s down-sampling procedures and quantification of activity counts can be found in the publication:
Neishabouri A, Nguyen J, Samuelsson J, et al. Quantification of acceleration as activity counts in ActiGraph wearable. Sci Rep. 2022;12(1):11958. 
The code uses these data files to compute 26 variables that reflect movement of the upper limbs in daily life. A description of each variable and the sampling frequency from which the variable was computed is displayed in the table below.
|Variable Name |	Variable Description	| Sampling Frequency|
|--------------|----------------------------|-------------------|
|recording_time |	The number of minutes of sensor recording time each day |	1 Hz|
|total_mvt_time |	Time (in hours) that either the left, right, or both upper limbs are moving	|1 Hz|
|l_time |	Time (in hours) that the left upper limb is moving |	1 Hz|
|l_only_time |	Time (in hours) that the left upper limb is moving, while the right upper limb is still|	1 Hz|
|r_time	| Time (in hours) that the right upper limb is moving	| 1 Hz|
|r_only_time|	Time (in hours) that the right upper limb is moving, while the left upper limb is still |1 Hz|
|simultaneous_time|	Time (in hours) that both upper limbs are moving	| 1 Hz|
|l_magnitude|	Median of the accelerations of the left upper limb, in activity counts|	1 Hz|
|r_magnitude|	Median of the accelerations of the right upper limb, in activity counts|	1 Hz|
|bilateral_magnitude|	Magnitude of accelerations of movement summed across both upper limbs, in activity counts|	1 Hz|
|l_magnitude_sd	| Standard deviation of the magnitude of accelerations across the left upper limb, in activity counts	| 1 Hz|
|r_magnitude_sd|	Standard deviation of the magnitude of accelerations across the right upper limb, in activity counts|	1 Hz|
|l_peak_magnitude|	The highest magnitude of accelerations of the left upper limb, in activity counts|	1 Hz|
|r_peak_magnitude|	The highest magnitude of accelerations of the right upper limb, in activity counts|	1 Hz|
|use_ratio|	Ratio of hours of non-dominant/affected upper limb movement, relative to hours od dominant/unaffected upper limb movement |	1 Hz|
|variation_ratio|	Ratio of the standard deviation of the non-dominant/affected upper limb relative to the dominant/unaffected upper limb |	1 Hz|
|simple_magnitude_ratio|	Ratio of the magnitude of accelerations of the non-dominant/affected upper limb relative to the dominant/unaffected upper limb |	1 Hz|
|l_entropy|	A measure of the time series variability from the accelerations of the left upper limb during the hour of maximum activity. Higher values indicate a more random signal. |	1 Hz|
|r_entropy|	A measure of the time series variability from the accelerations of the right upper limb during the hour of maximum activity. Higher values indicate a more random signal. |	1 Hz|
|l_jerk_ave|	The average jerk of the left upper limb. Higher values indicate less smooth movement.|	30 Hz|
|r_jerk_ave|	The average jerk of the right upper limb. Higher values indicate less smooth movement.|	30 Hz|
|jerk_aym|	The ratio of the average jerk magnitude between the non-dominant/affected upper limb and the dominant/unaffected upper limb.|	30 Hz|
|l_mean_freq|	The weighted mean of the component frequencies from the acceleration time series from the left upper limb.|	30 Hz|
|r_mean_freq|	The weighted mean of the component frequencies from the acceleration time series from the right upper limb.|	30 Hz|
|l_sd_freq|	The weighted standard deviation of the component frequencies from the acceleration time series from the left upper limb.|	30 Hz|
|r_sd_freq|	The weighted standard deviation of the component frequencies from the acceleration time series from the right upper limb.|	30 Hz|

A publication about the harmonized dataset is in preparation. The doi associated with this publication will be added to this README file, once available. We also encourage individuals who are interested in learning more about this dataset and the upper limb accelerometry variables to explore the R Shiny Object associated with the data at: https://langlab.shinyapps.io/harmonized_data/. At present, the R Shiny Object only includes data from pediatric participants. However, it is being updated to include data from adult participants.

Questions or feedback about the data or code can be directed to:
*Allison Miller, PT, DPT, PhD, NCS, Postdoctoral Research Associate, Washington University School of Medicine in St. Louis, miller.allison@wustl.edu
*Catherine Lang, PT, PhD, FASNR, FAPTA, Professor of Physical Therapy, Neurology, and Occupational Therapy, Washington University School of Medicine in St. Louis, langc@wustl.edu
*Keith Lohse, PhD, PStat, Associate Professor of Physical Therapy and Neurology, Affiliate of the Institute for Informatics, Data Science, and Biostatistics, Washington University School of Medicine in St. Louis, lohse@wustl.edu

