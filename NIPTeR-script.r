install.packages("NIPTeR")
library("NIPTeR")

##load already saved control sample and perform bin gc-correction..


NIPT_87_control_group <- readRDS( file = "/PATH/TO/FOLDER/NIPTeR_cleaned_87_control_group.rds" ) 
bingc_control_group <- gc_correct(nipt_object = NIPT_87_control_group, method = "bin")


## load raw sample and perform bin gc-correction..

NIPT_raw_sample_bin <- bin_bam_sample( bam_filepath = "/PATH/TO/FOLDER/Trisomy21.bam", separate_strands = T )
bingc_sample <- gc_correct(nipt_object = NIPT_raw_sample_bin, method = "bin")

##Now the χ2VR can be performed. Each sample needs to be corrected based on a control group. 
In our case, this is the bin GC corrected final control group. After correction both sample and control group are extracted.

NIPT_bin_chi_corrected_data <- chi_correct( nipt_sample = bingc_sample, nipt_control_group = bingc_control_group )
NIPT_bin_chi_corrected_sample <- NIPT_bin_chi_corrected_data$sample
NIPT_bin_chi_corrected_controls <- NIPT_bin_chi_corrected_data$control_group

#####Sample Match QC
"Now that the sample has been corrected in the same way as the control samples, we can calculate if the match score falls within the range of the control samples. "

"In this script control group has 87 samples. if you are using your own control group you can modify the parameter in the comand below. "n_of_samples=87" "


mean_match_sample <- mean(as.numeric(match_control_group(nipt_sample=NIPT_bin_chi_corrected_sample, nipt_control_group=NIPT_bin_chi_corrected_controls, n_of_samples=87, mode = "report")))


##Trisomy prediction
The most robust trisomy prediction is done if the same result is calculated using different algorithms (standard Z-score, Normalized Chromosome Value (NCV) and RBZ). If the different prediction methods are in agreement, the risk of a false positive result is lower. However, the sensitivity of the three  methods can differ and result in a positive result for one or two methods and a negative result for the other(s). The NCV and RBZ predictions use a subset of chromosomes to predict the number of reads expected to map on the chromosome of interest in case no trisomy is present. In general, these methods have a higher sensitivity than the standard Z-score prediction. However, they are more sensitive to bias. For instance, a maternal copy number variation in one of the chromosomes used as predictor can result in a false positive result. Therefore, to make this bias visible, for both the NCV and RBZ, we recommend performing more than one prediction using different predictor chromosomes to create the models. In case of conflicting results between models the cause of a positive result may not be a trisomy of the chromosome of interest, but rather an issue with one of the predictor chromosomes. 

#Standard Z-score
"The standard Z-score uses the fraction of reads mapped to the chromosome of interest compared with reads mapped to all other autosomes.
To calculate these Z-scores, the same function has to be called three times, each focusing on a different chromosome."

z_score_13 <- calculate_z_score(nipt_sample = NIPT_bin_chi_corrected_sample, nipt_control_group = NIPT_bin_chi_corrected_controls, chromo_focus = 13)

##to view Zscore. A Z-score higher than 3 indicates a trisomy

Z_score_13$sample_Zscore

###An important quality metric is the CV of the control groups, which can be calculated by dividing the standard deviation (sd) of reads mapped onto the chromosome of interest on the control samples by the mean of those values. A lower CV indicates a higher sensitivity.

as.numeric(as.vector(z_score_13$control_group_statistics[2]/z_score_13$control_group_statistics[1]))*100

###Normalized Chromosome Value (NCV)
"The second type of trisomy prediction included in NIPTeR is the NCV. In this method, the set of chromosomes that best predict the number of reads expected to map on the chromosome of interest are selected based on the control group. Using a higher number of denominators (maximum number of chromosomes used to predict the expected number of reads on the chromosome of interest) will mean a longer processing time. For each chromosome of interest, a different template is created. On default the function makes use of a train and a test set. This can lead to a slightly different set of predictors, coefficient of variation (CV) and sensitivity between calculations.
In our example we will use a maximum number of denominators of nine. We will not use a train and test set for NCV template creation here, so that results can be replicated. However, in practice, we recommend using a train and test set when sufficient control samples are available."

 ncv_template_13 <- prepare_ncv(nipt_control_group = NIPT_bin_chi_corrected_controls, chr_focus = 13, max_elements = 9, use_test_train_set = F)

##Now that the templates have been created, the NCVs can be calculated.

 ncv_score_13 <- calculate_ncv_score(nipt_sample = NIPT_bin_chi_corrected_sample, ncv_template = ncv_template_13)

#The results of the NCV predictions can be seen using the sample_score field. An NCV higher than 3 indicates a trisomy.
ncv_score_13$sample_score

##An important quality metric is the CV of the control groups, which can be calculated by dividing the standard deviation (sd) of reads mapped on the chromosome of interest on the control samples by the mean of those values (here shown for chromosome 13).

as.numeric(as.vector(ncv_template_13$control_group_statistics[2]/ncv_template_13$control_group_statistics[1]))*100

##Regression-based Z-score (RBZ)
"The third method for trisomy prediction that NIPTeR offers is the RBZ. This method uses a linear regression model using the most informative combinations of chromosomes to predict the expected number of reads on the chromosome of interest. On default the function makes use of a train and a test set, similar to what is done in calculating the NCV score. However, if the number of control samples is small the whole control group can be used to create the model. The user should realize that setting the use_test_train_set to ‘False’ creates a risk of overfitting the model to the control group. In this example, the complete control group is used to create the prediction models, so that results can be replicated. However, in practice we recommend always using a train and test set if the control group size allows it."

RBZ_13 <- perform_regression( nipt_sample = NIPT_bin_chi_corrected_sample, nipt_control_group = NIPT_bin_chi_corrected_controls, use_test_train_set = F, chromo_focus = 13 )

The output of the perform regression method can seem complex at first glance, because all information is presented in a single table. 

RBZ_13$prediction_statistics

"The prediction statistics show four different prediction models, each showing the RBZ and the control group statistics., -1.646 and -1.307; for chromosome 18 they are: -0.944, -0.238, -1.419 and -1.369; and for chromosome 21 they are: 11.188, 8.678, 9.698 and 8.927. For each of the chromosomes, the four models are in agreement and either all fall within the -3 to +3 range or all give a Z-score above this range. 
The second row shows the CV as compared to 1. To get the %CV, this number can be multiplied by 100 as we have done for the standard Z-score and NCV predictions (here shown for chromosome 13 prediction set 1)."

> as.numeric(as.vector(RBZ_13$prediction_statistics[2,1]))*100

"The third row in the prediction statistics states if the CV is the real practical CV as calculated using the control group, or a theoretical CV that we use to prevent overfitting of the model. Here, the practical CV is used in all cases, except for the chromosome 13 prediction set 1. This theoretical CV is based on the expected variation if no bias was present, multiplied by the overdispersion rate, which we have set to 15% by default. If the practical CV is lower than the theoretical (minimal possible) CV, the theoretical CV will be used to calculate the Z-score. This will lower the sensitivity of the prediction model slightly, but prevent false positive results. Since the predictors with the strongest correlation are selected for the first prediction model, this model has the highest risk for overfitting. This is also reflected in the increasing CVs going from prediction set 1 to 4. However, when a train and test set is used, the risk of overfitting the models is much smaller.
The fourth row shows the P value of the Shapiro Wilk test performed on the control group. A value below 0.05 indicates that the control group is not normal distributed and that the calculated Z-score should not be used. 
The fifth row shows the chromosomes used to create the prediction models. This information can be used to assess discordant results, as discussed in the NCV section above.
The sixth row shows the mean of the test set (or complete set when use_test_train_set = F is used), which should always be close to 1.
The seventh row shows the CV of the train set. This value is identical to the practical CV."

#################
This secript made by Abdullah Al-Nawfal. 
It is NOT tested yet due to sample limitation.
You do NOT need to validate the result (Like WES Benchmark), intstead of benchmark the result use the Control samples (30 min, 100 optimal)
Make sure to use Control samples form tested pouplation.
For any inquiries: a.alnawfal.1992@gmail.com


