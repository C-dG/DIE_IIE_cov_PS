README file for datasets for de Groot et al. 2025 

Date:	01.04.2025

============================================================
# General information

Publication title:  Cross-year repeatability and among-trait covariance for direct and indirect individual effects in producer-scrounger behaviour in wild house sparrows

Corresponding author information:
Name: 		Corné de Groot, 
Institution:	LMU München/Department Biologie II, Behavioural ecology group, Großhaderner Str. 2, 82152 Planegg-Martinsried
Email adress: 	c.degroot@bio.lmu.de
ORCID:		https://orcid.org/0000-0002-0355-7408

Data Collection range: 		2022-2023
Location of data collection: 	Lauvøya, Trondelag, Norway

For methods of data collection and analyses, please see details in the paper

Software used for analyses: 
R	V4.3.1 		(https://www.r-project.org/)
stan	V 2.32.2 	(https://mc-stan.org/rstan/)

All other packages and their versions are saved in the lock file

============================================================
# Description of variables in the dataset used in the analyses of the paper

Variable		Levels	Type		Description
RingNR			245	Random		The unique ringnumber combination of the focal individual
Opponent1_ID		245	Random		The unique ringnumber combination of opponent/social partner 1
Opponent2_ID		245	Random		The unique ringnumber combination of opponent/social partner 2
TrialID			1898	Random		Identifier of a given trial where individuals play the producer-scrounger game
Group_ID		53	Random		Identifier of a given group (6 individuals caught at the same farm)
Producing_events	Count	Response 	The number of producing events within a trial, where a focal individual has claimed a resource patch (well) by searching and was the first to find food in the well
Scrounging_events	Count	Response 	The number of scrounging events within a trial, where a focal individual has claimed a resource patch (well) by joining another individual directly at the well
TrialNR_Indiv		10	Fixed		The i-th trial of the focal individual within a given trial day
Trial_day_bin		2	Fixed		Binary factor describing whether the triadic trial was on the first or second trial day
Grp_tr_b4		2	Fixed		Binary factor describing whether the focal individual participated in a larger group trial before the triadic trials
bmr_YN			2	Fixed		Binary factor describing whether individuals were included in basal metabolic rate (BMR) measurments
Year			2	Fixed		Year where the trial took place
Series_ID_year		307	Random		Series combined of the focal individuals ringnumber (ID) and the year of that observation
Series_opp1_year	307	Random		Series combined of the ringnumber of opponent 1 and the year of that observation 
Series_opp2_year	307	Random		Series combined of the ringnumber of opponent 2 and the year of that observation
