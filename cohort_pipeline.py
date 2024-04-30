

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.8ab3b599-4e8b-40b5-a0ba-fe58e4a2fd4c"),
    concept_set_members=Input(rid="ri.foundry.main.dataset.e670c5ad-42ca-46a2-ae55-e917e3e161b6"),
    condition_occurrence=Input(rid="ri.foundry.main.dataset.900fa2ad-87ea-4285-be30-c6b5bab60e86"),
    location=Input(rid="ri.foundry.main.dataset.efac41e8-cc64-49bf-9007-d7e22a088318"),
    manifest=Input(rid="ri.foundry.main.dataset.b1e99f7f-5dcd-4503-985a-bbb28edc8f6f"),
    measurement=Input(rid="ri.foundry.main.dataset.d6054221-ee0c-4858-97de-22292458fa19"),
    microvisit_to_macrovisit_lds=Input(rid="ri.foundry.main.dataset.5af2c604-51e0-4afa-b1ae-1e5fa2f4b905"),
    person=Input(rid="ri.foundry.main.dataset.50cae11a-4afb-457d-99d4-55b4bc2cbe66")
)
# Confirmed Covid Cohort (677f370e-7059-4c78-87fa-3f1b1eb0b265): v51
#* COHORT:
#*   inputs: [measurement, concept_set_members, person, location, manifest, condition_occurrence, microvisit_to_macrovisit_lds]
#*   ext: py 
#*   desc: "identifies cohort of patients"

#Purpose - The purpose of this pipeline is to produce a day level and a persons level fact table for the confirmed (positive COVID PCR or AG lab test or U07.1 diagnosis charted) COVID positive patients in the N3C enclave. More information can be found in the README linked here (https://unite.nih.gov/workspace/notepad/view/ri.notepad.main.notepad.eb0ce7bb-a20b-4e81-bcef-aeb919bd350a).
#Creator/Owner/contact - Andrea Zhou
#Last Update - 12/6/23
#Description - This node identifies all patients with positive results from a PCR or AG COVID-19 lab test and the date of the patients' first instance of this type of COVID-19+ test.  It also identifies all patients with a COVID-19 diagnosis charted and the date of the patients’ first instance of this type of diagnosis (when available).  The earlier of the two is considered the index date for downstream calculations.  This transform then gathers some commonly used facts about these patients from the "person" and "location" tables, as well as some facts about the patient's institution (from the "manifest" table).  Available age, race, ethnicity, and locations data (including SDOH variables for L3 only) is gathered at this node.  The patient’s number of visits before and after covid as well as the number of days in their observation period before and after covid is calculated from the “microvisits_to_macrovisits” table in this node after filtering for plausible visit dates based on a site's latest data_extraction_date.  These facts will eventually be joined with the final patient-level table in the final node.

def COHORT(measurement, concept_set_members, person, location, manifest, condition_occurrence, microvisit_to_macrovisit_lds):
        
    """
    Select proportion of enclave patients to use: A value of 1.0 indicates the pipeline will use all patients in the persons table.  
    A value less than 1.0 takes a random sample of the patients with a value of 0.001 (for example) representing a 0.1% sample of the persons table will be used.
    """
    proportion_of_patients_to_use = 1.0

    concepts_df = concept_set_members
    person_sample = (person
        .select('person_id','year_of_birth','month_of_birth','day_of_birth','ethnicity_concept_name','race_concept_name','gender_concept_name','location_id','data_partner_id')
        .withColumnRenamed('gender_concept_name', 'sex')
        .distinct()
        .sample(False, proportion_of_patients_to_use, 111))
    measurement_df = (measurement
        .select('person_id', 'measurement_date', 'measurement_concept_id', 'value_as_concept_id')
        .where(measurement.measurement_date.isNotNull())
        .join(person_sample, 'person_id', 'inner'))
    conditions_df = (condition_occurrence
        .select('person_id', 'condition_start_date', 'condition_concept_id')
        .where(condition_occurrence.condition_start_date.isNotNull())
        .join(person_sample, 'person_id','inner'))
    manifest_df = manifest.select('data_partner_id','run_date','cdm_name','cdm_version','shift_date_yn','max_num_shift_days').withColumnRenamed("run_date", "data_extraction_date")   
    visits_df = microvisit_to_macrovisit_lds.select("person_id", "data_partner_id", "macrovisit_start_date", "visit_start_date") \
        .join(manifest_df.select('data_partner_id', 'data_extraction_date'), on='data_partner_id', how='left') \
        .where(
        (F.col('visit_start_date') >= "2018-01-01") &
        (F.col('visit_start_date') < (F.col('data_extraction_date')+(365*2)))
    ).drop('data_partner_id', 'data_extraction_date')
    location_df = location.dropDuplicates(subset=['location_id']).select('location_id','city','state','zip','county').withColumnRenamed('zip','postal_code')   
     
    """
    make list of concept IDs for Covid tests and filter measurements table for only these concept IDs
    then make list of concept IDs for POSITIVE Covid tests and label covid test measurements table as 1 for pos covid tests concept IDs and 0 for not
    """
    covid_measurement_test_ids = list(concepts_df.where((concepts_df.concept_set_name=="ATLAS SARS-CoV-2 rt-PCR and AG") & (concepts_df.is_most_recent_version=='true')).select('concept_id').toPandas()['concept_id'])
    covid_positive_measurement_ids = list(concepts_df.where((concepts_df.concept_set_name=="ResultPos") & (concepts_df.is_most_recent_version=='true')).select('concept_id').toPandas()['concept_id'])
    measurements_of_interest = measurement_df.where(measurement_df.measurement_concept_id.isin(covid_measurement_test_ids)) 
    measurements_of_interest = measurements_of_interest.where(measurements_of_interest.value_as_concept_id.isin(covid_positive_measurement_ids)).withColumnRenamed("measurement_date","covid_measurement_date").dropDuplicates(subset=['person_id','covid_measurement_date']).select('person_id','covid_measurement_date')
    first_covid_pos_lab = measurements_of_interest.groupBy('person_id').agg(F.min('covid_measurement_date').alias('COVID_first_PCR_or_AG_lab_positive'))
    
    # add flag for first date of COVID-19 diagnosis code if available
    COVID_concept_ids = list(concepts_df.where((concepts_df.concept_set_name=="N3C Covid Diagnosis") & (concepts_df.is_most_recent_version=='true')).select('concept_id').toPandas()['concept_id'])
    conditions_of_interest = conditions_df.where(conditions_df.condition_concept_id.isin(COVID_concept_ids)).withColumnRenamed("condition_start_date","covid_DIAGNOSIS_date").dropDuplicates(subset=['person_id','covid_DIAGNOSIS_date']).select('person_id','covid_DIAGNOSIS_date')
    first_covid_DIAGNOSIS = conditions_of_interest.groupBy('person_id').agg(F.min('covid_DIAGNOSIS_date').alias('COVID_first_diagnosis_date'))
    
    #join lab positive with diagnosis positive to create all confirmed covid patients cohort
    df = first_covid_pos_lab.join(first_covid_DIAGNOSIS, 'person_id', 'outer')
    #add a column for the earlier of the diagnosis or the lab test dates for all confirmed covid patients
    df = df.withColumn("COVID_first_poslab_or_diagnosis_date", F.least(df.COVID_first_PCR_or_AG_lab_positive, df.COVID_first_diagnosis_date))

    #add in demographics+locations data for all confirmed covid patients
    df = df.join(person_sample, 'person_id', 'inner')
    #join in location_df data for all confirmed covid patients
    df = df.join(location_df, 'location_id','left')

    #join in manifest_df information
    df = df.join(manifest_df, 'data_partner_id','inner')
    df = df.withColumn('max_num_shift_days', F.when(F.col('max_num_shift_days')=="", F.lit('0')).otherwise(F.regexp_replace(F.lower('max_num_shift_days'), 'na', '0')))
    
    #calculate date of birth for all confirmed covid patients
    df = df.withColumn("new_year_of_birth", F.when(F.col('year_of_birth').isNull(),1)
                                                .otherwise(F.col('year_of_birth')))
    df = df.withColumn("new_month_of_birth", F.when(F.col('month_of_birth').isNull(), 7)
                                                .when(F.col('month_of_birth')==0, 7)
                                                .otherwise(F.col('month_of_birth')))
    df = df.withColumn("new_day_of_birth", F.when(F.col('day_of_birth').isNull(), 1)
                                                .when(F.col('day_of_birth')==0, 1)
                                                .otherwise(F.col('day_of_birth')))

    df = df.withColumn("date_of_birth", F.concat_ws("-", F.col("new_year_of_birth"), F.col("new_month_of_birth"), F.col("new_day_of_birth")))
    df = df.withColumn("date_of_birth", F.to_date("date_of_birth", format=None)) 

    #convert date of birth string to date and apply min and max reasonable birthdate filter parameters, inclusive
    max_shift_as_int = df.withColumn("shift_days_as_int", F.col('max_num_shift_days').cast(IntegerType())) \
        .select(F.max('shift_days_as_int')) \
        .head()[0]
    min_reasonable_dob = "1902-01-01"
    max_reasonable_dob = F.date_add(F.current_date(), max_shift_as_int)
    df = df.withColumn("date_of_birth", F.when(F.col('date_of_birth').between(min_reasonable_dob, max_reasonable_dob), F.col('date_of_birth')).otherwise(None))

    #df = df.withColumn("age", F.floor(F.months_between(max_reasonable_dob, "date_of_birth", roundOff=False)/12))
    df = df.withColumn("age_at_covid", F.floor(F.months_between("COVID_first_poslab_or_diagnosis_date", "date_of_birth", roundOff=False)/12))

    H = ['Hispanic']
    A = ['Asian', 'Asian Indian', 'Bangladeshi', 'Bhutanese', 'Burmese', 'Cambodian', 'Chinese', 'Filipino', 'Hmong', 'Indonesian', 'Japanese', 'Korean', 'Laotian', 'Malaysian', 'Maldivian', 'Nepalese', 'Okinawan', 'Pakistani', 'Singaporean', 'Sri Lankan', 'Taiwanese', 'Thai', 'Vietnamese']
    B_AA = ['African', 'African American', 'Barbadian', 'Black', 'Black or African American', 'Dominica Islander', 'Haitian', 'Jamaican', 'Madagascar', 'Trinidadian', 'West Indian']
    W = ['White']
    NH_PI = ['Melanesian', 'Micronesian', 'Native Hawaiian or Other Pacific Islander', 'Other Pacific Islander', 'Polynesian']
    AI_AN = ['American Indian or Alaska Native']
    O = ['More than one race', 'Multiple race', 'Multiple races', 'Other', 'Other Race']
    U = ['Asian or Pacific islander', 'No Information', 'No matching concept', 'Refuse to Answer', 'Unknown', 'Unknown racial group']

    df = df.withColumn("race", F.when(F.col("race_concept_name").isin(H), "Hispanic or Latino")
                        .when(F.col("race_concept_name").isin(A), "Asian")
                        .when(F.col("race_concept_name").isin(B_AA), "Black or African American")
                        .when(F.col("race_concept_name").isin(W), "White")
                        .when(F.col("race_concept_name").isin(NH_PI), "Native Hawaiian or Other Pacific Islander") 
                        .when(F.col("race_concept_name").isin(AI_AN), "American Indian or Alaska Native")
                        .when(F.col("race_concept_name").isin(O), "Other")
                        .when(F.col("race_concept_name").isin(U), "Unknown")
                        .otherwise("Unknown"))

    df = df.withColumn("race_ethnicity", F.when(F.col("ethnicity_concept_name") == 'Hispanic or Latino', "Hispanic or Latino Any Race")
                        .when(F.col("race_concept_name").isin(H), "Hispanic or Latino Any Race")
                        .when(F.col("race_concept_name").isin(A), "Asian Non-Hispanic")
                        .when(F.col("race_concept_name").isin(B_AA), "Black or African American Non-Hispanic")
                        .when(F.col("race_concept_name").isin(W), "White Non-Hispanic")
                        .when(F.col("race_concept_name").isin(NH_PI), "Native Hawaiian or Other Pacific Islander Non-Hispanic") 
                        .when(F.col("race_concept_name").isin(AI_AN), "American Indian or Alaska Native Non-Hispanic")
                        .when(F.col("race_concept_name").isin(O), "Other Non-Hispanic")
                        .when(F.col("race_concept_name").isin(U), "Unknown")
                        .otherwise("Unknown"))
                        
    #create visit counts/obs period for before and post COVID 
    hosp_visits = visits_df.where(F.col("macrovisit_start_date").isNotNull()) \
        .orderBy("visit_start_date") \
        .coalesce(1) \
        .dropDuplicates(["person_id", "macrovisit_start_date"]) #hospital
    non_hosp_visits = visits_df.where(F.col("macrovisit_start_date").isNull()) \
        .dropDuplicates(["person_id", "visit_start_date"]) #non-hospital
    visits_df = hosp_visits.union(non_hosp_visits) #join the two

    """
    join in earliest index date value and use to calculate datediff between lab and visit 
    if positive then date is before the PCR/AG+ date
    if negative then date is after the PCR/AG+ date
    """
    visits_df = visits_df \
        .join(df.select('person_id','COVID_first_poslab_or_diagnosis_date','shift_date_yn','max_num_shift_days'), 'person_id', 'inner') \
        .withColumn('earliest_index_minus_visit_start_date', F.datediff('COVID_first_poslab_or_diagnosis_date','visit_start_date'))

    #counts for visits before
    visits_before = visits_df.where(F.col('earliest_index_minus_visit_start_date') > 0) \
        .groupBy("person_id") \
        .count() \
        .select("person_id", F.col('count').alias('number_of_visits_before_covid')) 
    #obs period in days before, where earliest_index_minus_visit_start_date = 0 means the pt_max_visit_date is the index date
    observation_before = visits_df.where(F.col('earliest_index_minus_visit_start_date') >= 0) \
        .groupby('person_id').agg(
        F.max('visit_start_date').alias('pt_max_visit_date'),
        F.min('visit_start_date').alias('pt_min_visit_date')) \
        .withColumn('observation_period_before_covid', F.datediff('pt_max_visit_date', 'pt_min_visit_date')) \
        .select('person_id', 'observation_period_before_covid')

    #counts for visits after
    visits_post = visits_df.where(F.col('earliest_index_minus_visit_start_date') < 0) \
        .groupBy("person_id") \
        .count() \
        .select("person_id", F.col('count').alias('number_of_visits_post_covid'))
    #obs period in days after, where earliest_index_minus_visit_start_date = 0 means the pt_min_visit_date is the index date
    observation_post = visits_df.where(F.col('earliest_index_minus_visit_start_date') <= 0) \
        .groupby('person_id').agg(
        F.max('visit_start_date').alias('pt_max_visit_date'),
        F.min('visit_start_date').alias('pt_min_visit_date')) \
        .withColumn('observation_period_post_covid', F.datediff('pt_max_visit_date', 'pt_min_visit_date')) \
        .select('person_id', 'observation_period_post_covid')
    
    #join visit counts/obs periods dataframes with main dataframe
    df = df.join(visits_before, "person_id", "left")
    df = df.join(observation_before, "person_id", "left")
    df = df.join(visits_post, "person_id", "left")
    df = df.join(observation_post, "person_id", "left")

    #LEVEL 2 ONLY
    #df = df.withColumn('max_num_shift_days', F.concat(F.col('max_num_shift_days'), F.lit(" + 180"))).withColumn('shift_date_yn', F.lit('Y'))

    df = df.select(
        'person_id',
        'COVID_first_PCR_or_AG_lab_positive',
        'COVID_first_diagnosis_date',
        'COVID_first_poslab_or_diagnosis_date',
        'number_of_visits_before_covid',
        'observation_period_before_covid',
        'number_of_visits_post_covid',
        'observation_period_post_covid',
        'sex',
        'city',
        'state',
        'postal_code',
        'county',
        'age_at_covid',
        'race',
        'race_ethnicity',
        'data_partner_id',
        'data_extraction_date',
        'cdm_name',
        'cdm_version',
        'shift_date_yn',
        'max_num_shift_days')

    return df
    

#################################################
## Global imports and functions included below ##
#################################################
from pyspark.sql import functions as F
from pyspark.sql.types import IntegerType

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.41ea332f-c59a-4b3d-9ccf-529854028b78"),
    COHORT=Input(rid="ri.foundry.main.dataset.8ab3b599-4e8b-40b5-a0ba-fe58e4a2fd4c"),
    COVID_deaths=Input(rid="ri.foundry.main.dataset.2ef7ea5d-3bb3-44b0-be55-ceb61cce404a"),
    Sdoh_variables_all_patients=Input(rid="ri.foundry.main.dataset.0ba44b7c-2167-4845-a647-97decc41e41c"),
    cohort_all_facts_table=Input(rid="ri.foundry.main.dataset.eaf478f9-7ad2-4118-8ae5-7aac86457cf8"),
    customize_concept_sets=Input(rid="ri.foundry.main.dataset.c3113ddd-37b3-42ed-84d6-d5c17088d238"),
    visits_of_interest=Input(rid="ri.foundry.main.dataset.ba5af22d-7022-450e-9652-f196fea9aaa2")
)
# Covid-19 patient summary fact table (60eab80f-463d-43bb-a48e-481fb0a05db7): v42
#* COVID_Patient_Summary_Table_LDS:
#*   inputs: [cohort_all_facts_table, COHORT, visits_of_interest, COVID_deaths, customize_concept_sets, Sdoh_variables_all_patients]
#*   ext: py 
#*   desc: "creates summary table with one row per patient"

#Purpose - The purpose of this pipeline is to produce a day level and a persons level fact table for the confirmed (positive COVID PCR or AG lab test or U07.1 diagnosis charted) COVID positive patients in the N3C enclave.
#Creator/Owner/contact - Andrea Zhou
#Last Update - 12/6/23
#Description - The final step is to aggregate information to create a data frame that contains a single row of data for each patient in the cohort.  This node aggregates all information from the cohort_all_facts_table and summarizes each patient's facts in a single row.  The patient’s hospitalization length of stay is calculated in this node.  For patients with ED visits and/or hospitalizations concurrent with their positive COVID-19 index date, indicators are created in this node. This transformation then joins the before or day of COVID, during hospitalization, and post COVID indicator data frames on the basis of unique patients.

def COVID_Patient_Summary_Table(cohort_all_facts_table, COHORT, visits_of_interest, COVID_deaths, customize_concept_sets, Sdoh_variables_all_patients):

    SDOH_vars = Sdoh_variables_all_patients #dataset not joined in, just referenced here so that it will import with template for easy access if user wants to join to their summary tables

    visits_df = visits_of_interest
    deaths_df = COVID_deaths.select('person_id','COVID_patient_death')
    all_facts = cohort_all_facts_table
    fusion_sheet = customize_concept_sets

    pre_columns = list(
        fusion_sheet.filter(fusion_sheet.pre_during_post.contains('pre'))
        .select('indicator_prefix')
        .distinct().toPandas()['indicator_prefix'])
    pre_columns.extend(['person_id', 'BMI_rounded', 'Antibody_Pos', 'Antibody_Neg', 'had_vaccine_administered'])
    during_columns = list(
        fusion_sheet.filter(fusion_sheet.pre_during_post.contains('during'))
        .select('indicator_prefix')
        .distinct().toPandas()['indicator_prefix'])
    during_columns.extend(['person_id', 'COVID_patient_death'])
    post_columns = list(
        fusion_sheet.filter(fusion_sheet.pre_during_post.contains('post'))
        .select('indicator_prefix')
        .distinct().toPandas()['indicator_prefix'])
    post_columns.extend(['person_id', 'BMI_rounded', 'PCR_AG_Pos', 'PCR_AG_Neg', 'Antibody_Pos', 'Antibody_Neg', 'is_first_reinfection', 'had_vaccine_administered'])

    df_pre_COVID = all_facts \
        .where(all_facts.pre_COVID==1) \
        .select(list(set(pre_columns) & set(all_facts.columns)))
    df_during_strong_COVID_hospitalization = all_facts \
        .where(all_facts.during_first_strong_COVID_hospitalization==1) \
        .select(list(set(during_columns) & set(all_facts.columns)))
    df_during_weak_COVID_hospitalization = all_facts \
        .where(all_facts.during_first_weak_COVID_hospitalization==1) \
        .select(list(set(during_columns) & set(all_facts.columns)))
    df_post_COVID = all_facts \
        .where(all_facts.post_COVID==1) \
        .select(list(set(post_columns) & set(all_facts.columns)))
   
    df_pre_COVID = df_pre_COVID.groupby('person_id').agg(
        F.max('BMI_rounded').alias('BMI_max_observed_or_calculated_before_or_day_of_covid'),
        *[F.max(col).alias(col + '_before_or_day_of_covid_indicator') for col in df_pre_COVID.columns if col not in ('person_id', 'BMI_rounded', 'had_vaccine_administered')],
        F.sum('had_vaccine_administered').alias('number_of_COVID_vaccine_doses_before_or_day_of_covid'))
    
    df_during_strong_COVID_hospitalization = df_during_strong_COVID_hospitalization.groupby('person_id').agg(
        *[F.max(col).alias(col + '_during_strong_covid_hospitalization_indicator') for col in df_during_strong_COVID_hospitalization.columns if col not in ('person_id')])

    df_during_weak_COVID_hospitalization = df_during_weak_COVID_hospitalization.groupby('person_id').agg(
        *[F.max(col).alias(col + '_during_weak_covid_hospitalization_indicator') for col in df_during_weak_COVID_hospitalization.columns if col not in ('person_id')])

    df_post_COVID = df_post_COVID.groupby('person_id').agg(
        F.max('BMI_rounded').alias('BMI_max_observed_or_calculated_post_covid'),
        *[F.max(col).alias(col + '_post_covid_indicator') for col in df_post_COVID.columns if col not in ('person_id', 'BMI_rounded', 'is_first_reinfection', 'had_vaccine_administered')],
        F.sum('had_vaccine_administered').alias('number_of_COVID_vaccine_doses_post_covid'),
        F.max('is_first_reinfection').alias('had_at_least_one_reinfection_post_covid_indicator'))

    #join above four tables on patient ID 
    df = df_pre_COVID.join(df_during_strong_COVID_hospitalization, 'person_id', 'outer')
    df = df.join(df_during_weak_COVID_hospitalization, 'person_id', 'outer')
    df = df.join(df_post_COVID, 'person_id', 'outer')
    
    df = df.join(visits_df,'person_id', 'outer')

    #already dependent on decision made in visits of interest node, no changes necessary here
    df = df.withColumn('strong_COVID_hospitalization_length_of_stay', 
        F.datediff("first_strong_COVID_hospitalization_end_date", "first_strong_COVID_hospitalization_start_date"))
    
    #join back in generic death flag for any patient with or without a date
    df = df.join(deaths_df, 'person_id', 'left').withColumnRenamed('COVID_patient_death', 'COVID_patient_death_indicator')
    #join back in death within fixed window post covid for patients with a date to use in severity of index infection
    df = df.join(all_facts.select('person_id','death_within_specified_window_post_covid').where(F.col('death_within_specified_window_post_covid')==1), 'person_id', 'left')
    #join in demographics and manifest data from cohort node
    df = COHORT.join(df, 'person_id','left')
    
    df = df.na.fill(value=0, subset = [col for col in df.columns if col not in ('BMI_max_observed_or_calculated_before_or_day_of_covid','BMI_max_observed_or_calculated_post_covid', 'postal_code', 'age_at_covid')])

    df = df.withColumn("Severity_Type", 
        F.when((df.COVID_first_PCR_or_AG_lab_positive.isNull() & df.COVID_first_diagnosis_date.isNull()), "No_COVID_index")
        .when((df.death_within_specified_window_post_covid == 1), "Death_within_n_days_after_COVID_index")
        .when((df.LL_ECMO_during_strong_covid_hospitalization_indicator == 1) | (df.LL_IMV_during_strong_covid_hospitalization_indicator == 1), "Severe_ECMO_IMV_in_Hosp_around_strong_signal_COVID_index")
        .when(df.first_strong_COVID_hospitalization_start_date.isNotNull(), "Moderate_Hosp_around_strong_signal_COVID_index")
        .when(df.first_strong_COVID_ED_only_start_date.isNotNull(), "Mild_ED_around_strong_signal_COVID_index")
        .when((df.LL_ECMO_during_weak_covid_hospitalization_indicator == 1) | (df.LL_IMV_during_weak_covid_hospitalization_indicator == 1), "Severe_ECMO_IMV_in_Hosp_around_weak_signal_COVID_index")
        .when(df.first_weak_COVID_hospitalization_start_date.isNotNull(), "Moderate_Hosp_around_weak_signal_COVID_index")
        .when(df.first_weak_COVID_ED_only_start_date.isNotNull(), "Mild_ED_around_weak_signal_COVID_index")
        .otherwise("Mild_No_ED_or_Hosp_around_COVID_index"))

    
    return df
        
#################################################
## Global imports and functions included below ##
#################################################

from pyspark.sql import functions as F

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.2ef7ea5d-3bb3-44b0-be55-ceb61cce404a"),
    COHORT=Input(rid="ri.foundry.main.dataset.8ab3b599-4e8b-40b5-a0ba-fe58e4a2fd4c"),
    concept_set_members=Input(rid="ri.foundry.main.dataset.e670c5ad-42ca-46a2-ae55-e917e3e161b6"),
    death=Input(rid="ri.foundry.main.dataset.d8cc2ad4-215e-4b5d-bc80-80ffb3454875"),
    microvisit_to_macrovisit_lds=Input(rid="ri.foundry.main.dataset.5af2c604-51e0-4afa-b1ae-1e5fa2f4b905")
)
# Covid deaths (f583218a-7f0c-4c9a-b500-8330a8ca99c9): v12
#* COVID_deaths:
#*   inputs: [death, COHORT, microvisit_to_macrovisit_lds, concept_set_members]
#*   ext: py 
#*   desc: "finds death concepts"

#Purpose - The purpose of this pipeline is to produce a day level and a persons level fact table for the confirmed (positive COVID PCR or AG lab test or U07.1 diagnosis charted) COVID positive patients in the N3C enclave.
#Creator/Owner/contact - Andrea Zhou
#Last Update - 12/6/23
#Description - This node filters the visits table for rows that have a discharge_to_concept_id that corresponds with the DECEASED or HOSPICE concept sets and combines these records with the patients in the deaths table. Death dates are taken from the deaths table and from the visits table if the patient has a discharge_to_concept_id that corresponds with the DECEASED concept set. Death dates are prioritized such that we filter to only plausible death dates from the OMOP death table and take a min prior to adding in any additional plausible death dates from the OMOP visits table with concept id belonging to the DECEASED concept set and take a max for patients who do not already have a death date from the OMOP death table. No date is retained for patients who were discharged to hospice.

def COVID_deaths(death, COHORT, microvisit_to_macrovisit_lds, concept_set_members):
 
    persons = COHORT.select('person_id', 'COVID_first_poslab_or_diagnosis_date', 'data_extraction_date')
    death_df = death \
        .select('person_id', 'death_date') \
        .distinct() \
        .join(persons, on = 'person_id', how = 'inner')
    visits_df = microvisit_to_macrovisit_lds \
        .select('person_id', 'visit_start_date', 'visit_end_date', 'discharge_to_concept_id') \
        .distinct() \
        .join(persons, on = 'person_id', how = 'inner')
    concepts_df = concept_set_members \
        .select('concept_set_name', 'is_most_recent_version', 'concept_id') \
        .where(F.col('is_most_recent_version')=='true')

    #create lists of concept ids to look for in the discharge_to_concept_id column of the visits_df
    death_from_visits_ids = list(concepts_df.where(F.col('concept_set_name') == "DECEASED").select('concept_id').toPandas()['concept_id'])
    hospice_from_visits_ids = list(concepts_df.where(F.col('concept_set_name') == "HOSPICE").select('concept_id').toPandas()['concept_id'])

    #filter visits table to patient and date rows that have DECEASED that matches list of concept_ids
    death_from_visits_df = visits_df \
        .where(F.col('discharge_to_concept_id').isin(death_from_visits_ids)) \
        .drop('discharge_to_concept_id') \
        .distinct()
    #filter visits table to patient rows that have HOSPICE that matches list of concept_ids
    hospice_from_visits_df = visits_df.drop('visit_start_date', 'visit_end_date', 'data_extraction_date') \
        .where(F.col('discharge_to_concept_id').isin(hospice_from_visits_ids)) \
        .drop('discharge_to_concept_id') \
        .distinct()

    #####################################################################
    ###combine relevant visits sourced deaths with deaths table deaths###
    #####################################################################
   
    #keep rows where death_date is plausible 
    #then take earliest recorded date per person
    death_w_dates_df = death_df.where(
        (F.col('death_date') >= "2018-01-01") &
        (F.col('death_date') < (F.col('data_extraction_date')+(365*2)))
    ).groupby('person_id').agg(F.min('death_date').alias('death_date')) \
    .drop('data_extraction_date')

    #from rows of visit table that have concept_id belonging to DECEASED concept set,
    #create new column visit_death_date that is plauisble visit_end_date when available,
    #and plausible visit_start_date when visit_end_date is Null or not plausible
    #then take latest recorded date per person
    death_from_visits_w_dates_df = death_from_visits_df \
        .withColumn('visit_death_date', F.when(
            (F.col('visit_end_date') >= "2018-01-01") &
            (F.col('visit_end_date') < (F.col('data_extraction_date')+(365*2))), F.col('visit_end_date')
            ).when(
            (F.col('visit_start_date') >= "2018-01-01") &
            (F.col('visit_start_date') < (F.col('data_extraction_date')+(365*2))), F.col('visit_start_date')
            ).otherwise(None)
        ).groupby('person_id').agg(F.max('visit_death_date').alias('visit_death_date')) \
        .drop('data_extraction_date')

    #join deaths with dates from both domains
    df = death_w_dates_df.join(death_from_visits_w_dates_df, on = 'person_id', how = 'outer')
    #prioritize death_dates from the deaths table over from visits table
    df = df.withColumn('date', F.when(F.col('death_date').isNotNull(), F.col('death_date')).otherwise(F.col('visit_death_date'))) \
        .drop('death_date', 'visit_death_date')

    #join in patients, without any date, from deaths table
    #join in patients, without any date, from visits table
    #join in patients, without any date, for HOSPICE from visits table
    #inner join with cohort node patients to keep only confirmed covid patients
    df = df.join(death_df.select('person_id', 'death_date'), on='person_id', how='outer') \
        .join(death_from_visits_df.select('person_id'), on='person_id', how='outer') \
        .join(hospice_from_visits_df, on='person_id', how='outer') \
        .join(persons, on = 'person_id', how = 'inner') \
        .dropDuplicates()
    #flag all patients as having died regardless of date
    df = df.withColumn("COVID_patient_death", F.lit(1))
    df = df.withColumn("30d_mort", F.lit())
    
    return df

#################################################
## Global imports and functions included below ##
#################################################

from pyspark.sql import functions as F

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.eaf478f9-7ad2-4118-8ae5-7aac86457cf8"),
    COHORT=Input(rid="ri.foundry.main.dataset.8ab3b599-4e8b-40b5-a0ba-fe58e4a2fd4c"),
    COVID_deaths=Input(rid="ri.foundry.main.dataset.2ef7ea5d-3bb3-44b0-be55-ceb61cce404a"),
    conditions_of_interest=Input(rid="ri.foundry.main.dataset.66e11f26-eab2-4ac9-a04d-e76561655fd6"),
    devices_of_interest=Input(rid="ri.foundry.main.dataset.2dfc7935-760c-43e0-a25d-18d1a2ac9335"),
    drugs_of_interest=Input(rid="ri.foundry.main.dataset.48b61b37-b9dc-4404-b8c6-f9f853a61197"),
    measurements_of_interest=Input(rid="ri.foundry.main.dataset.6c6d90fe-a756-49ff-a53f-676c73ab9d3d"),
    microvisit_to_macrovisit_lds=Input(rid="ri.foundry.main.dataset.5af2c604-51e0-4afa-b1ae-1e5fa2f4b905"),
    observations_of_interest=Input(rid="ri.foundry.main.dataset.2ab56165-057d-4989-8ec2-77889f70e9fc"),
    procedures_of_interest=Input(rid="ri.foundry.main.dataset.ffa4d627-676e-4f75-99eb-97b829c15402"),
    vaccines_of_interest=Input(rid="ri.foundry.main.dataset.6942249e-6273-4f3c-b428-1b5b845accef"),
    visits_of_interest=Input(rid="ri.foundry.main.dataset.ba5af22d-7022-450e-9652-f196fea9aaa2")
)
# Covid-19 patients visit day facts table (d87ac29f-6486-410d-97f3-8fed33d805f9): v36
#* cohort_all_facts_table:
#*   inputs: [conditions_of_interest, measurements_of_interest, visits_of_interest, procedures_of_interest, observations_of_interest, drugs_of_interest, COVID_deaths, COHORT, devices_of_interest, microvisit_to_macrovisit_lds, vaccines_of_interest]
#*   ext: py 
#*   desc: "creates date level table with one row per patient per date"

#Purpose - The purpose of this pipeline is to produce a day level and a persons level fact table for the confirmed (positive COVID PCR or AG lab test or U07.1 diagnosis charted) COVID positive patients in the N3C enclave.
#Creator/Owner/contact - Andrea Zhou
#Last Update - 12/6/23
#Description - All facts collected in the previous steps are combined in this cohort_all_facts_table on the basis of unique days for each patient and logic is applied to see if the instance of the coded fact appeared in the EHR prior to or after the patient's first COVID-19 positive PCR or AG lab test.  Indicators are created for the presence or absence of events, medications, conditions, measurements, device exposures, observations, procedures, and outcomes, either occurring before COVID index date, during the patient’s hospitalization, or in the period after COVID index date.  It also creates an indicator for whether the date where a fact was noted occurred during any hospitalization, not just the COVID associated ones found in the visits of interest node.  A flag for the visit in which the patient is noted to have their first covid reinfection is also thrown in this node. The default time range is 60 days after index date, but the number of days can be specified by parameter input based on researcher interests. This table is useful if the analyst needs to use actual dates of events as it provides more detail than the final patient-level table.  Use the max and min functions to find the first and last occurrences of any events.

def cohort_all_facts_table(conditions_of_interest, measurements_of_interest, visits_of_interest, procedures_of_interest, observations_of_interest, drugs_of_interest, COVID_deaths, COHORT, devices_of_interest, microvisit_to_macrovisit_lds, vaccines_of_interest):

    persons_df = COHORT.select('person_id', 'COVID_first_PCR_or_AG_lab_positive', 'COVID_first_diagnosis_date', 'COVID_first_poslab_or_diagnosis_date').dropDuplicates()
    macrovisits_df = microvisit_to_macrovisit_lds
    vaccines_df = vaccines_of_interest
    procedures_df = procedures_of_interest
    devices_df = devices_of_interest
    observations_df = observations_of_interest
    conditions_df = conditions_of_interest
    drugs_df = drugs_of_interest
    measurements_df = measurements_of_interest
    visits_df = visits_of_interest
    deaths_df = COVID_deaths.where(
        (COVID_deaths.date.isNotNull()) 
        & (COVID_deaths.date >= "2018-01-01") 
        & (COVID_deaths.date < (F.col('data_extraction_date')+(365*2)))) \
        .drop('data_extraction_date')

    df = macrovisits_df.select('person_id','visit_start_date').withColumnRenamed('visit_start_date','date')
    df = df.join(vaccines_df, on=list(set(df.columns)&set(vaccines_df.columns)), how='outer')
    df = df.join(procedures_df, on=list(set(df.columns)&set(procedures_df.columns)), how='outer')
    df = df.join(devices_df, on=list(set(df.columns)&set(devices_df.columns)), how='outer')
    df = df.join(observations_df, on=list(set(df.columns)&set(observations_df.columns)), how='outer')
    df = df.join(conditions_df, on=list(set(df.columns)&set(conditions_df.columns)), how='outer')
    df = df.join(drugs_df, on=list(set(df.columns)&set(drugs_df.columns)), how='outer')
    df = df.join(measurements_df, on=list(set(df.columns)&set(measurements_df.columns)), how='outer')    
    df = df.join(deaths_df, on=list(set(df.columns)&set(deaths_df.columns)), how='outer')

    df = df.na.fill(value=0, subset = [col for col in df.columns if col not in ('BMI_rounded')])
   
    #add F.max of all indicator columns to collapse all cross-domain flags to unique person and visit rows
    df = df.groupby('person_id', 'date').agg(*[F.max(col).alias(col) for col in df.columns if col not in ('person_id','date')])
   
    #join persons
    df = persons_df.join(df, 'person_id', 'left')
    df = visits_df.join(df, 'person_id', 'outer') 

    #create reinfection indicator, minimum 90 day window from index date to subsequent positive test
    reinfection_wait_time = 90
    reinfection_df = df.withColumn('is_reinfection', 
        F.when(( (F.col('PCR_AG_Pos')==1) & (F.datediff(F.col('date'), F.col('COVID_first_poslab_or_diagnosis_date')) > reinfection_wait_time) ), 1).otherwise(0)) \
        .where(F.col('is_reinfection')==1) \
        .groupby('person_id') \
        .agg(F.min('date').alias('date'), 
            F.max('is_reinfection').alias('is_first_reinfection'))
    df = df.join(reinfection_df, on=['person_id','date'], how='left')

    #create new death within fixed window post COVID flag to be used for severity calculation, same window as looking for reinfection 
    df = df.withColumn('death_within_specified_window_post_covid', F.when(
        (F.col('COVID_patient_death') ==1) &
        (F.datediff("date", "COVID_first_poslab_or_diagnosis_date") > 0) 
        & (F.datediff("date", "COVID_first_poslab_or_diagnosis_date") < 30.1), 1).otherwise(0))
    
    #defaulted to find the lesser date value of the first lab positive result date and the first diagnosis date, could be adjusted to only "COVID_first_diagnosis_date" or only "COVID_first_PCR_or_AG_lab_positive" based on desired index event definition
    df = df.withColumn('pre_COVID', F.when(F.datediff("COVID_first_poslab_or_diagnosis_date","date")>=0, 1).otherwise(0))
    df = df.withColumn('post_COVID', F.when(F.datediff("COVID_first_poslab_or_diagnosis_date","date")<0, 1).otherwise(0))

    #dependent on the definition chosen in the visits of interest node, no changes necessary here
    df = df.withColumn('during_first_strong_COVID_hospitalization', F.when((F.datediff("first_strong_COVID_hospitalization_end_date","date")>=0) & (F.datediff("first_strong_COVID_hospitalization_start_date","date")<=0), 1).otherwise(0))
    df = df.withColumn('during_first_strong_COVID_ED_visit', F.when(F.datediff("first_strong_COVID_ED_only_start_date","date")==0, 1).otherwise(0))

    #weak ED and hospitalization from visits of interest node, no changes necessary here
    df = df.withColumn('during_first_weak_COVID_hospitalization', F.when((F.datediff("first_weak_COVID_hospitalization_end_date","date")>=0) & (F.datediff("first_weak_COVID_hospitalization_start_date","date")<=0), 1).otherwise(0))
    df = df.withColumn('during_first_weak_COVID_ED_visit', F.when(F.datediff("first_weak_COVID_ED_only_start_date","date")==0, 1).otherwise(0))

    #drop dates for all facts table once indicators are created for 'during_first_COVID_hospitalization'
    df = df.drop('first_strong_COVID_hospitalization_end_date','first_strong_COVID_ED_only_start_date', 'first_weak_COVID_hospitalization_end_date', 'first_weak_COVID_ED_only_start_date', 'macrovisit_start_date', 'macrovisit_end_date')

    #create and join in flag that indicates whether the visit was during a macrovisit (1) or not (0)
    #any conditions, observations, procedures, devices, drugs, measurements, and/or death flagged 
    #with a (1) on that particular visit date would then be considered to have happened during a macrovisit    
    macrovisits_df = macrovisits_df \
        .select('person_id', 'macrovisit_start_date', 'macrovisit_end_date') \
        .where(F.col('macrovisit_start_date').isNotNull() & F.col('macrovisit_end_date').isNotNull()) \
        .distinct()
    df_hosp = df.select('person_id', 'date').distinct() \
        .join(macrovisits_df, on='person_id', how= 'outer')
    df_hosp = df_hosp.withColumn('during_macrovisit_hospitalization', F.when(F.col('date').between(F.col('macrovisit_start_date'), F.col('macrovisit_end_date')), 1).otherwise(0)) \
        .where(F.col('during_macrovisit_hospitalization') == 1) \
        .dropDuplicates(['person_id', 'date'])
    df = df.join(df_hosp, on=['person_id','date'], how="left")
    
    #final fill of null non-continuous variables with 0
    df = df.na.fill(value=0, subset = [col for col in df.columns if col not in ('BMI_rounded')])

    return df
    
#################################################
## Global imports and functions included below ##
#################################################

from pyspark.sql import functions as F

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.66e11f26-eab2-4ac9-a04d-e76561655fd6"),
    COHORT=Input(rid="ri.foundry.main.dataset.8ab3b599-4e8b-40b5-a0ba-fe58e4a2fd4c"),
    concept_set_members=Input(rid="ri.foundry.main.dataset.e670c5ad-42ca-46a2-ae55-e917e3e161b6"),
    condition_occurrence=Input(rid="ri.foundry.main.dataset.900fa2ad-87ea-4285-be30-c6b5bab60e86"),
    customize_concept_sets=Input(rid="ri.foundry.main.dataset.c3113ddd-37b3-42ed-84d6-d5c17088d238")
)
# Conditions of Interest (a48a1f27-e400-411d-bee8-b4d3f0e2faa2): v25
#* conditions_of_interest:
#*   inputs: [COHORT, concept_set_members, condition_occurrence, customize_concept_sets]
#*   ext: py 
#*   desc: "finds condition concepts"

#Purpose - The purpose of this pipeline is to produce a day level and a persons level fact table for the confirmed (positive COVID PCR or AG lab test or U07.1 diagnosis charted) COVID positive patients in the N3C enclave.
#Creator/Owner/contact - Andrea Zhou
#Last Update - 12/7/22
#Description - This node filters the condition_occurences table for rows that have a condition_concept_id associated with one of the concept sets described in the data dictionary in the README through the use of a fusion sheet.  Indicator names for these conditions are assigned, and the indicators are collapsed to unique instances on the basis of patient and date.

def conditions_of_interest(COHORT, concept_set_members, condition_occurrence, customize_concept_sets):

    #bring in only cohort patient ids
    persons = COHORT.select('person_id')
    #filter observations table to only cohort patients    
    conditions_df = condition_occurrence \
        .select('person_id', 'condition_start_date', 'condition_concept_id') \
        .where(F.col('condition_start_date').isNotNull()) \
        .withColumnRenamed('condition_start_date','date') \
        .withColumnRenamed('condition_concept_id','concept_id') \
        .join(persons,'person_id','inner')

    #filter fusion sheet for concept sets and their future variable names that have concepts in the conditions domain
    fusion_df = customize_concept_sets \
        .filter(customize_concept_sets.domain.contains('condition')) \
        .select('concept_set_name','indicator_prefix')
    #filter concept set members table to only concept ids for the conditions of interest
    concepts_df = concept_set_members \
        .select('concept_set_name', 'is_most_recent_version', 'concept_id') \
        .where(F.col('is_most_recent_version')=='true') \
        .join(fusion_df, 'concept_set_name', 'inner') \
        .select('concept_id','indicator_prefix')

    #find conditions information based on matching concept ids for conditions of interest
    df = conditions_df.join(concepts_df, 'concept_id', 'inner')
    #collapse to unique person and visit date and pivot on future variable name to create flag for rows associated with the concept sets for conditions of interest    
    df = df.groupby('person_id','date').pivot('indicator_prefix').agg(F.lit(1)).na.fill(0)
   
    return df
    

#################################################
## Global imports and functions included below ##
#################################################

from pyspark.sql import functions as F

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.c3113ddd-37b3-42ed-84d6-d5c17088d238"),
    LL_DO_NOT_DELETE_REQUIRED_concept_sets_confirmed=Input(rid="ri.foundry.main.dataset.794f85ce-746c-4eb7-a12d-2c2718b5dc4a"),
    molnu_concept_sets_fusion=Input(rid="ri.foundry.main.dataset.4241f1d5-bed2-4b3e-a94b-be71d36244a4")
)
# customize_concept_sets (9849ee10-fa74-44a5-b20d-c8c6e39d7314): v13
#* customize_concept_sets:
#*   inputs: [LL_concept_sets_fusion, LL_DO_NOT_DELETE_REQUIRED_concept_sets_confirmed]
#*   ext: py 
#*   desc: "combines the two fusion inputs"

#Description - The purpose of this node is to optimize the user's experience connecting a customized concept set "fusion sheet" input data frame to replace molnu_concept_sets_fusion.

def customize_concept_sets( LL_DO_NOT_DELETE_REQUIRED_concept_sets_confirmed, molnu_concept_sets_fusion):

    required = LL_DO_NOT_DELETE_REQUIRED_concept_sets_confirmed
    customizable = molnu_concept_sets_fusion
    
    df = required.join(customizable, on = required.columns, how = 'outer')
    
    return df
    

#################################################
## Global imports and functions included below ##
#################################################

from pyspark.sql import functions as F

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.2dfc7935-760c-43e0-a25d-18d1a2ac9335"),
    COHORT=Input(rid="ri.foundry.main.dataset.8ab3b599-4e8b-40b5-a0ba-fe58e4a2fd4c"),
    concept_set_members=Input(rid="ri.foundry.main.dataset.e670c5ad-42ca-46a2-ae55-e917e3e161b6"),
    customize_concept_sets=Input(rid="ri.foundry.main.dataset.c3113ddd-37b3-42ed-84d6-d5c17088d238"),
    device_exposure=Input(rid="ri.foundry.main.dataset.d685db48-6583-43d6-8dc5-a9ebae1a827a")
)
# Devices of interest (3e6d32f1-7a08-4d8d-9d63-aad0e46287d4): v13
#* devices_of_interest:
#*   inputs: [device_exposure, COHORT, concept_set_members, customize_concept_sets]
#*   ext: py 
#*   desc: "finds device concepts"

#Purpose - The purpose of this pipeline is to produce a day level and a persons level fact table for the confirmed (positive COVID PCR or AG lab test or U07.1 diagnosis charted) COVID positive patients in the N3C enclave.
#Creator/Owner/contact - Andrea Zhou
#Last Update - 12/7/22
#Description - This nodes filter the source OMOP tables for rows that have a standard concept id associated with one of the concept sets described in the data dictionary in the README through the use of a fusion sheet.  Indicator names for these variables are assigned, and the indicators are collapsed to unique instances on the basis of patient and date.

def devices_of_interest(device_exposure, COHORT, concept_set_members, customize_concept_sets):

    #bring in only cohort patient ids
    persons = COHORT.select('person_id')
    #filter device exposure table to only cohort patients
    devices_df = device_exposure \
        .select('person_id','device_exposure_start_date','device_concept_id') \
        .where(F.col('device_exposure_start_date').isNotNull()) \
        .withColumnRenamed('device_exposure_start_date','date') \
        .withColumnRenamed('device_concept_id','concept_id') \
        .join(persons,'person_id','inner')

    #filter fusion sheet for concept sets and their future variable names that have concepts in the devices domain
    fusion_df = customize_concept_sets \
        .filter(customize_concept_sets.domain.contains('device')) \
        .select('concept_set_name','indicator_prefix')
    #filter concept set members table to only concept ids for the devices of interest
    concepts_df = concept_set_members \
        .select('concept_set_name', 'is_most_recent_version', 'concept_id') \
        .where(F.col('is_most_recent_version')=='true') \
        .join(fusion_df, 'concept_set_name', 'inner') \
        .select('concept_id','indicator_prefix')
        
    #find device exposure information based on matching concept ids for devices of interest
    df = devices_df.join(concepts_df, 'concept_id', 'inner')
    #collapse to unique person and visit date and pivot on future variable name to create flag for rows associated with the concept sets for devices of interest
    df = df.groupby('person_id','date').pivot('indicator_prefix').agg(F.lit(1)).na.fill(0)

    return df
    

#################################################
## Global imports and functions included below ##
#################################################

from pyspark.sql import functions as F
    

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.48b61b37-b9dc-4404-b8c6-f9f853a61197"),
    COHORT=Input(rid="ri.foundry.main.dataset.8ab3b599-4e8b-40b5-a0ba-fe58e4a2fd4c"),
    concept_set_members=Input(rid="ri.foundry.main.dataset.e670c5ad-42ca-46a2-ae55-e917e3e161b6"),
    customize_concept_sets=Input(rid="ri.foundry.main.dataset.c3113ddd-37b3-42ed-84d6-d5c17088d238"),
    drug_exposure=Input(rid="ri.foundry.main.dataset.ec252b05-8f82-4f7f-a227-b3bb9bc578ef")
)
# Drugs of interest (2f62711b-0797-4d28-933b-48a5c7d67add): v14
#* drugs_of_interest:
#*   inputs: [concept_set_members, drug_exposure, COHORT, customize_concept_sets]
#*   ext: py 
#*   desc: "finds drug concepts"

#Purpose - The purpose of this pipeline is to produce a day level and a persons level fact table for the confirmed (positive COVID PCR or AG lab test or U07.1 diagnosis charted) COVID positive patients in the N3C enclave.
#Creator/Owner/contact - Andrea Zhou
#Last Update - 12/7/22
#Description - This nodes filter the source OMOP tables for rows that have a standard concept id associated with one of the concept sets described in the data dictionary in the README through the use of a fusion sheet.  Indicator names for these variables are assigned, and the indicators are collapsed to unique instances on the basis of patient and date.

def drugs_of_interest(concept_set_members, drug_exposure, COHORT, customize_concept_sets):
  
    #bring in only cohort patient ids
    persons = COHORT.select('person_id')
    #filter drug exposure table to only cohort patients    
    drug_df = drug_exposure \
        .select('person_id','drug_exposure_start_date','drug_concept_id') \
        .where(F.col('drug_exposure_start_date').isNotNull()) \
        .withColumnRenamed('drug_exposure_start_date','date') \
        .withColumnRenamed('drug_concept_id','concept_id') \
        .join(persons,'person_id','inner')

    #filter fusion sheet for concept sets and their future variable names that have concepts in the drug domain
    fusion_df = customize_concept_sets \
        .filter(customize_concept_sets.domain.contains('drug')) \
        .select('concept_set_name','indicator_prefix')
    #filter concept set members table to only concept ids for the drugs of interest
    concepts_df = concept_set_members \
        .select('concept_set_name', 'is_most_recent_version', 'concept_id') \
        .where(F.col('is_most_recent_version')=='true') \
        .join(fusion_df, 'concept_set_name', 'inner') \
        .select('concept_id','indicator_prefix')
        
    #find drug exposure information based on matching concept ids for drugs of interest
    df = drug_df.join(concepts_df, 'concept_id', 'inner')
    #collapse to unique person and visit date and pivot on future variable name to create flag for rows associated with the concept sets for drugs of interest
    df = df.groupby('person_id','date').pivot('indicator_prefix').agg(F.lit(1)).na.fill(0)

    return df
    

#################################################
## Global imports and functions included below ##
#################################################

from pyspark.sql import functions as F

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.6c6d90fe-a756-49ff-a53f-676c73ab9d3d"),
    COHORT=Input(rid="ri.foundry.main.dataset.8ab3b599-4e8b-40b5-a0ba-fe58e4a2fd4c"),
    concept_set_members=Input(rid="ri.foundry.main.dataset.e670c5ad-42ca-46a2-ae55-e917e3e161b6"),
    measurement=Input(rid="ri.foundry.main.dataset.d6054221-ee0c-4858-97de-22292458fa19")
)
# Measurements of Interest (b37105f2-072b-4808-a11a-ecacc9fb6c8e): v13
#* measurements_of_interest:
#*   inputs: [measurement, concept_set_members, COHORT]
#*   ext: py 
#*   desc: "finds measurement concepts and values"

#Purpose - The purpose of this pipeline is to produce a day level and a persons level fact table for the confirmed (positive COVID PCR or AG lab test or U07.1 diagnosis charted) COVID positive patients in the N3C enclave.
#Creator/Owner/contact - Andrea Zhou
#Last Update - 12/7/22
#Description - This node filters the measurements table for rows that have a measurement_concept_id associated with one of the concept sets described in the data dictionary in the README.  It finds the harmonized value as a number for the quantitative measurements and collapses these values to unique instances on the basis of patient and date.  It also finds the value as concept id for the qualitative measurements (covid labs) and collapses these to unique instances on the basis of patient and date.  Measurement BMI cutoffs included are intended for adults. Analyses focused on pediatric measurements should use different bounds for BMI measurements.

def measurements_of_interest(measurement, concept_set_members, COHORT):
    
    #bring in only cohort patient ids
    persons = COHORT.select('person_id')
    #filter procedure occurrence table to only cohort patients    
    df = measurement \
        .select('person_id','measurement_date','measurement_concept_id','harmonized_value_as_number', 'value_as_concept_id') \
        .where(F.col('measurement_date').isNotNull()) \
        .withColumnRenamed('measurement_date','date') \
        .join(persons,'person_id','inner')
        
    concepts_df = concept_set_members \
        .select('concept_set_name', 'is_most_recent_version', 'concept_id') \
        .where(F.col('is_most_recent_version')=='true')
          
    #Find BMI closest to COVID using both reported/observed BMI and calculated BMI using height and weight.  Cutoffs for reasonable height, weight, and BMI are provided and can be changed by the template user.
    lowest_acceptable_BMI = 10
    highest_acceptable_BMI = 100
    lowest_acceptable_weight = 5 #in kgs
    highest_acceptable_weight = 300 #in kgs
    lowest_acceptable_height = .6 #in meters
    highest_acceptable_height = 2.43 #in meters

    bmi_codeset_ids = list(concepts_df.where(
        (concepts_df.concept_set_name=="body mass index") 
        & (concepts_df.is_most_recent_version=='true')
        ).select('concept_id').toPandas()['concept_id'])
    weight_codeset_ids = list(concepts_df.where(
        (concepts_df.concept_set_name=="Body weight (LG34372-9 and SNOMED)") 
        & (concepts_df.is_most_recent_version=='true')
        ).select('concept_id').toPandas()['concept_id'])
    height_codeset_ids = list(concepts_df.where(
        (concepts_df.concept_set_name=="Height (LG34373-7 + SNOMED)") 
        & (concepts_df.is_most_recent_version=='true')
        ).select('concept_id').toPandas()['concept_id'])
    
    pcr_ag_test_ids = list(concepts_df.where(
        (concepts_df.concept_set_name=="ATLAS SARS-CoV-2 rt-PCR and AG") 
        & (concepts_df.is_most_recent_version=='true')
        ).select('concept_id').toPandas()['concept_id'])
    antibody_test_ids = list(concepts_df.where(
        (concepts_df.concept_set_name=="Atlas #818 [N3C] CovidAntibody retry") 
        & (concepts_df.is_most_recent_version=='true')
        ).select('concept_id').toPandas()['concept_id'])
    covid_positive_measurement_ids = list(concepts_df.where(
        (concepts_df.concept_set_name=="ResultPos") 
        & (concepts_df.is_most_recent_version=='true')
        ).select('concept_id').toPandas()['concept_id'])
    covid_negative_measurement_ids = list(concepts_df.where(
        (concepts_df.concept_set_name=="ResultNeg") 
        & (concepts_df.is_most_recent_version=='true')
        ).select('concept_id').toPandas()['concept_id'])

    #add value columns for rows associated with the above concept sets, but only include BMI or height or weight when in reasonable range
    BMI_df = df.where(F.col('harmonized_value_as_number').isNotNull()) \
        .withColumn('Recorded_BMI', F.when(df.measurement_concept_id.isin(bmi_codeset_ids) & df.harmonized_value_as_number.between(lowest_acceptable_BMI, highest_acceptable_BMI), df.harmonized_value_as_number).otherwise(0)) \
        .withColumn('height', F.when(df.measurement_concept_id.isin(height_codeset_ids) & df.harmonized_value_as_number.between(lowest_acceptable_height, highest_acceptable_height), df.harmonized_value_as_number).otherwise(0)) \
        .withColumn('weight', F.when(df.measurement_concept_id.isin(weight_codeset_ids) & df.harmonized_value_as_number.between(lowest_acceptable_weight, highest_acceptable_weight), df.harmonized_value_as_number).otherwise(0)) 
        
    labs_df = df.withColumn('PCR_AG_Pos', F.when(df.measurement_concept_id.isin(pcr_ag_test_ids) & df.value_as_concept_id.isin(covid_positive_measurement_ids), 1).otherwise(0)) \
        .withColumn('PCR_AG_Neg', F.when(df.measurement_concept_id.isin(pcr_ag_test_ids) & df.value_as_concept_id.isin(covid_negative_measurement_ids), 1).otherwise(0)) \
        .withColumn('Antibody_Pos', F.when(df.measurement_concept_id.isin(antibody_test_ids) & df.value_as_concept_id.isin(covid_positive_measurement_ids), 1).otherwise(0)) \
        .withColumn('Antibody_Neg', F.when(df.measurement_concept_id.isin(antibody_test_ids) & df.value_as_concept_id.isin(covid_negative_measurement_ids), 1).otherwise(0))
     
    #collapse all reasonable values to unique person and visit rows
    BMI_df = BMI_df.groupby('person_id', 'date').agg(
    F.max('Recorded_BMI').alias('Recorded_BMI'),
    F.max('height').alias('height'),
    F.max('weight').alias('weight'))
    labs_df = labs_df.groupby('person_id', 'date').agg(
    F.max('PCR_AG_Pos').alias('PCR_AG_Pos'),
    F.max('PCR_AG_Neg').alias('PCR_AG_Neg'),
    F.max('Antibody_Pos').alias('Antibody_Pos'),
    F.max('Antibody_Neg').alias('Antibody_Neg'))

    #add a calculated BMI for each visit date when height and weight available.  Note that if only one is available, it will result in zero
    #subsequent filter out rows that would have resulted from unreasonable calculated_BMI being used as best_BMI for the visit 
    BMI_df = BMI_df.withColumn('calculated_BMI', (BMI_df.weight/(BMI_df.height*BMI_df.height)))
    BMI_df = BMI_df.withColumn('BMI', F.when(BMI_df.Recorded_BMI>0, BMI_df.Recorded_BMI).otherwise(BMI_df.calculated_BMI)) \
        .select('person_id','date','BMI')
    BMI_df = BMI_df.filter((BMI_df.BMI<=highest_acceptable_BMI) & (BMI_df.BMI>=lowest_acceptable_BMI)) \
        .withColumn('BMI_rounded', F.round(BMI_df.BMI)) \
        .drop('BMI')
    BMI_df = BMI_df.withColumn('OBESITY', F.when(BMI_df.BMI_rounded>=30, 1).otherwise(0))

    #join BMI_df with labs_df to retain all lab results with only reasonable BMI_rounded and OBESITY flags
    df = labs_df.join(BMI_df, on=['person_id', 'date'], how='left')

    return df

#################################################
## Global imports and functions included below ##
#################################################

from pyspark.sql import functions as F

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.2ab56165-057d-4989-8ec2-77889f70e9fc"),
    COHORT=Input(rid="ri.foundry.main.dataset.8ab3b599-4e8b-40b5-a0ba-fe58e4a2fd4c"),
    concept_set_members=Input(rid="ri.foundry.main.dataset.e670c5ad-42ca-46a2-ae55-e917e3e161b6"),
    customize_concept_sets=Input(rid="ri.foundry.main.dataset.c3113ddd-37b3-42ed-84d6-d5c17088d238"),
    observation=Input(rid="ri.foundry.main.dataset.b998b475-b229-471c-800e-9421491409f3")
)
# Observations of Interest (3669cbbf-59a0-40c5-af56-e95c99181a94): v16
#* observations_of_interest:
#*   inputs: [observation, concept_set_members, COHORT, customize_concept_sets]
#*   ext: py 
#*   desc: "finds observation concepts"

#Purpose - The purpose of this pipeline is to produce a day level and a persons level fact table for the confirmed (positive COVID PCR or AG lab test or U07.1 diagnosis charted) COVID positive patients in the N3C enclave.
#Creator/Owner/contact - Andrea Zhou
#Last Update - 12/7/22
#Description - This nodes filter the source OMOP tables for rows that have a standard concept id associated with one of the concept sets described in the data dictionary in the README through the use of a fusion sheet.  Indicator names for these variables are assigned, and the indicators are collapsed to unique instances on the basis of patient and date.

def observations_of_interest(observation, concept_set_members, COHORT, customize_concept_sets):

    #bring in only cohort patient ids
    persons = COHORT.select('person_id')
    #filter observations table to only cohort patients    
    observations_df = observation \
        .select('person_id','observation_date','observation_concept_id') \
        .where(F.col('observation_date').isNotNull()) \
        .withColumnRenamed('observation_date','date') \
        .withColumnRenamed('observation_concept_id','concept_id') \
        .join(persons,'person_id','inner')

    #filter fusion sheet for concept sets and their future variable names that have concepts in the observations domain
    fusion_df = customize_concept_sets \
        .filter(customize_concept_sets.domain.contains('observation')) \
        .select('concept_set_name','indicator_prefix')
    #filter concept set members table to only concept ids for the observations of interest
    concepts_df = concept_set_members \
        .select('concept_set_name', 'is_most_recent_version', 'concept_id') \
        .where(F.col('is_most_recent_version')=='true') \
        .join(fusion_df, 'concept_set_name', 'inner') \
        .select('concept_id','indicator_prefix')

    #find observations information based on matching concept ids for observations of interest
    df = observations_df.join(concepts_df, 'concept_id', 'inner')
    #collapse to unique person and visit date and pivot on future variable name to create flag for rows associated with the concept sets for observations of interest    
    df = df.groupby('person_id','date').pivot('indicator_prefix').agg(F.lit(1)).na.fill(0)

    return df

#################################################
## Global imports and functions included below ##
#################################################

from pyspark.sql import functions as F

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.ffa4d627-676e-4f75-99eb-97b829c15402"),
    COHORT=Input(rid="ri.foundry.main.dataset.8ab3b599-4e8b-40b5-a0ba-fe58e4a2fd4c"),
    concept_set_members=Input(rid="ri.foundry.main.dataset.e670c5ad-42ca-46a2-ae55-e917e3e161b6"),
    customize_concept_sets=Input(rid="ri.foundry.main.dataset.c3113ddd-37b3-42ed-84d6-d5c17088d238"),
    procedure_occurrence=Input(rid="ri.foundry.main.dataset.f6f0b5e0-a105-403a-a98f-0ee1c78137dc")
)
# Procedures of interest (54a01e81-ff63-4271-9e87-181af5de0d73): v13
#* procedures_of_interest:
#*   inputs: [COHORT, concept_set_members, procedure_occurrence, customize_concept_sets]
#*   ext: py 
#*   desc: "finds procedure concepts"

#Purpose - The purpose of this pipeline is to produce a day level and a persons level fact table for the confirmed (positive COVID PCR or AG lab test or U07.1 diagnosis charted) COVID positive patients in the N3C enclave.
#Creator/Owner/contact - Andrea Zhou
#Last Update - 12/7/22
#Description - This nodes filter the source OMOP tables for rows that have a standard concept id associated with one of the concept sets described in the data dictionary in the README through the use of a fusion sheet.  Indicator names for these variables are assigned, and the indicators are collapsed to unique instances on the basis of patient and date.

def procedures_of_interest(COHORT, concept_set_members, procedure_occurrence, customize_concept_sets):
  
    #bring in only cohort patient ids
    persons = COHORT.select('person_id')
    #filter procedure occurrence table to only cohort patients    
    procedures_df = procedure_occurrence \
        .select('person_id','procedure_date','procedure_concept_id') \
        .where(F.col('procedure_date').isNotNull()) \
        .withColumnRenamed('procedure_date','date') \
        .withColumnRenamed('procedure_concept_id','concept_id') \
        .join(persons,'person_id','inner')

    #filter fusion sheet for concept sets and their future variable names that have concepts in the procedure domain
    fusion_df = customize_concept_sets \
        .filter(customize_concept_sets.domain.contains('procedure')) \
        .select('concept_set_name','indicator_prefix')
    #filter concept set members table to only concept ids for the procedures of interest
    concepts_df = concept_set_members \
        .select('concept_set_name', 'is_most_recent_version', 'concept_id') \
        .where(F.col('is_most_recent_version')=='true') \
        .join(fusion_df, 'concept_set_name', 'inner') \
        .select('concept_id','indicator_prefix')
 
    #find procedure occurrence information based on matching concept ids for procedures of interest
    df = procedures_df.join(concepts_df, 'concept_id', 'inner')
    #collapse to unique person and visit date and pivot on future variable name to create flag for rows associated with the concept sets for procedures of interest    
    df = df.groupby('person_id','date').pivot('indicator_prefix').agg(F.lit(1)).na.fill(0)

    return df
    

#################################################
## Global imports and functions included below ##
#################################################

from pyspark.sql import functions as F

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.6942249e-6273-4f3c-b428-1b5b845accef"),
    COHORT=Input(rid="ri.foundry.main.dataset.8ab3b599-4e8b-40b5-a0ba-fe58e4a2fd4c"),
    Vaccine_fact_lds=Input(rid="ri.foundry.main.dataset.7482e426-55a2-4a0b-9976-4cb3aa35788d")
)
# vaccines_of_interest (45f42c26-bd5b-432b-becc-031e1fb09511): v5
#* vaccines_of_interest:
#*   inputs: [COHORT, Vaccine_fact_lds]
#*   ext: py 
#*   desc: "finds vaccination information"

#Purpose - The purpose of this pipeline is to produce a day level and a persons level fact table for the confirmed (positive COVID PCR or AG lab test or U07.1 diagnosis charted) COVID positive patients in the N3C enclave.
#Creator/Owner/contact - Andrea Zhou
#Last Update - 12/7/22
#Description - This node converts the vaccine_fact dataset from columns of dates for each dose number to a simple flag for dates on which a vaccine dose was received by the patient. The indicator is collapsed to unique instances on the basis of patient and date.

def vaccines_of_interest(COHORT, Vaccine_fact_lds):
    
    persons = COHORT.select('person_id')
    vax_df = Vaccine_fact_lds.select('person_id', '1_vax_date', '2_vax_date', '3_vax_date', '4_vax_date') \
        .join(persons, 'person_id', 'inner')

    first_dose = vax_df.select('person_id', '1_vax_date') \
        .withColumnRenamed('1_vax_date', 'date') \
        .where(F.col('date').isNotNull())
    second_dose = vax_df.select('person_id', '2_vax_date') \
        .withColumnRenamed('2_vax_date', 'date') \
        .where(F.col('date').isNotNull())        
    third_dose = vax_df.select('person_id', '3_vax_date') \
        .withColumnRenamed('3_vax_date', 'date') \
        .where(F.col('date').isNotNull())
    fourth_dose = vax_df.select('person_id', '4_vax_date') \
        .withColumnRenamed('4_vax_date', 'date') \
        .where(F.col('date').isNotNull())

    df = first_dose.join(second_dose, on=['person_id', 'date'], how='outer') \
        .join(third_dose, on=['person_id', 'date'], how='outer') \
        .join(fourth_dose, on=['person_id', 'date'], how='outer') \
        .distinct()

    df = df.withColumn('had_vaccine_administered', F.lit(1))

    return df

#################################################
## Global imports and functions included below ##
#################################################

from pyspark.sql import functions as F

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.ba5af22d-7022-450e-9652-f196fea9aaa2"),
    COHORT=Input(rid="ri.foundry.main.dataset.8ab3b599-4e8b-40b5-a0ba-fe58e4a2fd4c"),
    concept_set_members=Input(rid="ri.foundry.main.dataset.e670c5ad-42ca-46a2-ae55-e917e3e161b6"),
    microvisit_to_macrovisit_lds=Input(rid="ri.foundry.main.dataset.5af2c604-51e0-4afa-b1ae-1e5fa2f4b905")
)
# Visits of interest (99c478ba-8671-4b6c-9183-45d8a1059606): v25
#* visits_of_interest:
#*   inputs: [COHORT, microvisit_to_macrovisit_lds, concept_set_members]
#*   ext: py 
#*   desc: "finds ED and hosp visit concepts within desired timeframe"

#Purpose - The purpose of this pipeline is to produce a visit day level and a persons level fact table for the confirmed (positive COVID PCR or AG lab test or U07.1 diagnosis charted) COVID positive patients in the N3C enclave.
#Creator/Owner/contact - Andrea Zhou
#Last Update - 12/6/23
#Description - This node queries the microvisits_to_macrovisits table to identify hospitalizations.  The input table can be changed to the visits table if indicated, but the code in the transform would need to be modified accordingly.  The parameter called covid_associated_hospitalization_requires_lab_AND_diagnosis is created and allows the user to easily change whether they define COVID-19 associated ED visits and hospitalizations using the CDC definition (lab positive with a COVID-19 diagnosis charted) OR using anyone who is either lab positive or has a COVID-19 diagnosis charted.  Number of days between a patient’s diagnosis date and their positive lab result is also calculated in this node. As of 12/6/23, ED visits and hospitalizations that meet the CDC definition are considered "strong_COVID_associated" while ED visits and hospitalizations that only have a positive lab or a covid diagnosis are considered "weak_COVID"associated" visits. 

def visits_of_interest(COHORT, microvisit_to_macrovisit_lds, concept_set_members):
    
    #select test/dx date columns for cohort patients and add column for date diff between positive lab test and COVID diagnosis when available
    persons = COHORT \
        .select('person_id', 'COVID_first_PCR_or_AG_lab_positive', 'COVID_first_diagnosis_date', 'COVID_first_poslab_or_diagnosis_date') \
        .withColumn('lab_minus_diagnosis_date', F.datediff('COVID_first_PCR_or_AG_lab_positive','COVID_first_diagnosis_date'))
    #filter macrovisit table to only cohort patients    
    df = microvisit_to_macrovisit_lds \
        .select('person_id','visit_start_date','visit_concept_id','macrovisit_start_date','macrovisit_end_date') \
        .join(persons,'person_id','inner')  

    concepts_df = concept_set_members \
        .select('concept_set_name', 'is_most_recent_version', 'concept_id') \
        .where(F.col('is_most_recent_version')=='true')  

    # use macrovisit table to find ED only visits (that do not lead to hospitalization)   
    ED_concept_ids = list(concepts_df.where((concepts_df.concept_set_name=="[PASC] ED Visits") & (concepts_df.is_most_recent_version=='true')).select('concept_id').toPandas()['concept_id'])
    df_ED = df.where(df.macrovisit_start_date.isNull()&(df.visit_concept_id.isin(ED_concept_ids)))
    df_ED = df_ED.withColumn('lab_minus_ED_visit_start_date', F.datediff('COVID_first_PCR_or_AG_lab_positive','visit_start_date'))
    df_ED = df_ED.withColumn('diagnosis_minus_ED_visit_start_date', F.datediff('COVID_first_diagnosis_date','visit_start_date'))

    """
    create parameter for toggling COVID-19 related ED only visit and hospital admission definitions
    when parameter =True: Per CDC definitions of a COVID-19 associated ED or hospital admission visit, ensure that a COVID-19 diagnosis and ED/hospital admission occurred in the 16 days after or 1 day prior to the PCR or AG positive test (index event).
    when parameter =False: ED or hospital admission visits flagged based on the first instance of a positive COVID-19 PCR or AG lab result OR the first instance of a charted COVID-19 diagnosis when there is no positive lab result within specified timeframe of ED/hospital admission.
    """
    covid_associated_ED_or_hosp_requires_lab_AND_diagnosis = True
    num_days_before_index = -1
    num_days_after_index = 30
    

    if covid_associated_ED_or _hosp_requires_lab_AND_diagnosis:
        df_ED = (df_ED.withColumn('covid_pcr_or_ag_associated_ED_only_visit', F.when(F.col('lab_minus_ED_visit_start_date').between(-num_days_after_index,num_days_before_index), 1).otherwise(0))
                .withColumn('COVID_lab_positive_and_diagnosed_ED_visit', F.when((F.col('covid_pcr_or_ag_associated_ED_only_visit')==1) & (F.col('lab_minus_diagnosis_date').between(-num_days_after_index,num_days_before_index)), 1).otherwise(0))
                .withColumn('strong_covid_ED_only_start_date', F.when(F.col('COVID_lab_positive_and_diagnosed_ED_visit')==1, F.col('visit_start_date')).otherwise(None))
                .dropDuplicates())
    else:
        df_ED = df_ED.withColumn("earliest_index_minus_ED_start_date", F.datediff("COVID_first_poslab_or_diagnosis_date","visit_start_date"))
        #first lab or diagnosis date based, ED only visit
        df_ED = (df_ED.withColumn("covid_lab_or_dx_associated_ED_only_visit", F.when(F.col('earliest_index_minus_ED_start_date').between(-num_days_after_index,num_days_before_index), 1).otherwise(0))
                .withColumn('strong_covid_ED_only_start_date', F.when(F.col('covid_lab_or_dx_associated_ED_only_visit')==1, F.col('visit_start_date')).otherwise(None))
                .dropDuplicates())

    # flag visits where positive covid lab or diagnosis recorded on or the day after an ED only visit
    # limiting duration of interest to single day and overnight ED stays since less than 0.014% of covid positive patients with an ED only visit have duration between 2 to 7 days
    df_ED = df_ED.withColumn('weak_covid_associated_ED_only_visit', 
            F.when((F.col('lab_minus_ED_visit_start_date').between(-1,0)) | (F.col('diagnosis_minus_ED_visit_start_date').between(-1,0)), 1).otherwise(0)) \
        .withColumn('weak_covid_ED_only_start_date', 
            F.when(F.col('weak_covid_associated_ED_only_visit')==1, F.col('visit_start_date'))
            .otherwise(None))

    # use macrovisit table to find visits associated with hospitalization
    df_hosp = df.where(df.macrovisit_start_date.isNotNull())
    df_hosp = df_hosp.withColumn("lab_minus_hosp_start_date", F.datediff("COVID_first_PCR_or_AG_lab_positive","macrovisit_start_date")) 
    
    if covid_associated_ED_or_hosp_requires_lab_AND_diagnosis:
        df_hosp = (df_hosp.withColumn("covid_pcr_or_ag_associated_hospitalization", F.when(F.col('lab_minus_hosp_start_date').between(-num_days_after_index,num_days_before_index), 1).otherwise(0))
                .withColumn("covid_lab_positive_and_diagnosed_hospitalization", F.when((F.col('covid_pcr_or_ag_associated_hospitalization')==1) & (F.col('lab_minus_diagnosis_date').between(-num_days_after_index,num_days_before_index)), 1).otherwise(0))
                .withColumn('strong_covid_hospitalization_start_date', F.when(F.col('covid_lab_positive_and_diagnosed_hospitalization')==1, F.col('macrovisit_start_date')).otherwise(None))
                .withColumn('strong_covid_hospitalization_end_date', F.when(F.col('covid_lab_positive_and_diagnosed_hospitalization')==1, F.col('macrovisit_end_date')).otherwise(None))
                .dropDuplicates())
    else:
        df_hosp = df_hosp.withColumn("earliest_index_minus_hosp_start_date", F.datediff("COVID_first_poslab_or_diagnosis_date","macrovisit_start_date")) 

        #first lab or diagnosis date based, hospitalization visit
        df_hosp = (df_hosp.withColumn("covid_lab_or_diagnosis_associated_hospitilization", F.when(F.col('earliest_index_minus_hosp_start_date').between(-num_days_after_index,num_days_before_index), 1).otherwise(0))
                .withColumn('strong_covid_hospitalization_start_date', F.when(F.col('covid_lab_or_diagnosis_associated_hospitilization')==1, F.col('macrovisit_start_date')).otherwise(None))
                .withColumn('strong_covid_hospitalization_end_date', F.when(F.col('covid_lab_or_diagnosis_associated_hospitilization')==1, F.col('macrovisit_end_date')).otherwise(None))
                .dropDuplicates())

    # flag visits where positive covid lab or diagnosis recorded between macrovisit start and end date
    df_hosp = df_hosp.withColumn('weak_covid_associated_hospitalization', 
            F.when((F.col('COVID_first_diagnosis_date').between(F.col('macrovisit_start_date'), F.col('macrovisit_end_date'))) | (F.col('COVID_first_PCR_or_AG_lab_positive').between(F.col('macrovisit_start_date'), F.col('macrovisit_end_date'))), 1).otherwise(0)) \
        .withColumn('weak_covid_hospitalization_start_date', 
            F.when(F.col('weak_covid_associated_hospitalization')==1, F.col('macrovisit_start_date'))
            .otherwise(None)) \
        .withColumn('weak_covid_hospitalization_end_date', 
            F.when(F.col('weak_covid_associated_hospitalization')==1, F.col('macrovisit_end_date'))
            .otherwise(None))
 
    #join ED and hosp dataframes
    df = df.join(df_ED,'person_id', 'outer')
    df = df.join(df_hosp,'person_id', 'outer')
    
    #collapse all values to one row per person
    df = df.groupby('person_id').agg(
    F.min('strong_covid_ED_only_start_date').alias('first_strong_COVID_ED_only_start_date'),
    F.min('strong_covid_hospitalization_start_date').alias('first_strong_COVID_hospitalization_start_date'),
    F.min('strong_covid_hospitalization_end_date').alias('first_strong_COVID_hospitalization_end_date'),
    F.min('weak_covid_ED_only_start_date').alias('first_weak_COVID_ED_only_start_date'), 
    F.min('weak_covid_hospitalization_start_date').alias('first_weak_COVID_hospitalization_start_date'),
    F.min('weak_covid_hospitalization_end_date').alias('first_weak_COVID_hospitalization_end_date'))

    return df
    
#################################################
## Global imports and functions included below ##
#################################################

from pyspark.sql import functions as F

