

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.fcfa3b2d-8509-48df-8e6d-c6ec0cdb3130"),
    molnu_std=Input(rid="ri.foundry.main.dataset.338c64a0-0274-4587-aca6-f6e45d2e196a")
)
# table1_ITM_other (7d2527b0-09d0-439c-9558-13083dae73f7): v2
Table1_molnu_v_std <- function(molnu_std) {

    library(tableone)
    library(dplyr)
    library(tibble)

    options(width=1000)

    all_variables_in_table1_rows <- c("sex", "age_at_covid", 'race_ethnicity', 'COVID_first_poslab_or_diagnosis_date', 'mildliver', 'rheum', 'chf', 'mi', 'cad', 'htn', 'ckd', 'cancer', 'metcancer', 'dmcomp', 'dmuncomp', 'obesity', 'lungdisease', 'transplant', 'hiv', 'immunocomp', 'systemicsteroids', 'biologic', 'tobacco', 'subuse', 'vaccine_number', 'number_of_visits_before_covid', 'number_of_visits_post_covid', 'observation_period_before_covid', 'observation_period_post_covid')

    categorical_variables_in_table1_rows <- c('sex','race_ethnicity', 'mildliver', 'rheum', 'chf', 'mi', 'cad', 'htn', 'ckd', 'cancer', 'metcancer', 'dmcomp', 'dmuncomp', 'obesity', 'lungdisease', 'transplant', 'hiv', 'immunocomp', 'systemicsteroids', 'biologic', 'tobacco', 'subuse')

    table_1 <- CreateTableOne(vars = all_variables_in_table1_rows, strata = "MOLNUPIRAVIR", data = molnu_std, factorVars = categorical_variables_in_table1_rows, addOverall = FALSE)

    df = as.data.frame(print(table_1, smd = TRUE, showAllLevels = FALSE, explain = FALSE, pDigits = 4, formatOptions = list(big.mark = ","))) %>% rownames_to_column("Name")

    drops <- c("test")

    df <- df[ , !(names(df) %in% drops)]

    return(df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.3e192625-95f0-4393-987a-6ad8efb54084"),
    m_ps=Input(rid="ri.foundry.main.dataset.ed7d0600-aece-4d29-8f2a-252033739f2f")
)
# table2
Table2_molnu_v_std <- function(m_ps) {

    library(tableone)
    library(dplyr)
    library(survey)

m_ps_survey <- svydesign(ids=~1, data = m_ps, weights = ~ipw)

    options(width=1000)

    all_variables_in_table1_rows <- c("sex", "age_at_covid", 'race_ethnicity', 'COVID_first_poslab_or_diagnosis_date', 'mildliver', 'rheum', 'chf', 'mi', 'cad', 'htn', 'ckd', 'cancer', 'metcancer', 'dmcomp', 'dmuncomp', 'obesity', 'lungdisease', 'transplant', 'hiv', 'immunocomp', 'systemicsteroids', 'biologic', 'tobacco', 'subuse', 'vaccine_number', 'number_of_visits_before_covid', 'number_of_visits_post_covid', 'observation_period_before_covid', 'observation_period_post_covid')

    categorical_variables_in_table1_rows <- c('sex','race_ethnicity', 'mildliver', 'rheum', 'chf', 'mi', 'cad', 'htn', 'ckd', 'cancer', 'metcancer', 'dmcomp', 'dmuncomp', 'obesity', 'lungdisease', 'transplant', 'hiv', 'immunocomp', 'systemicsteroids', 'biologic', 'tobacco', 'subuse')

    table_1 <- svyCreateTableOne(vars = all_variables_in_table1_rows, strata = "MOLNUPIRAVIR", data = m_ps_survey, factorVars = categorical_variables_in_table1_rows, addOverall = FALSE)

    df = as.data.frame(print(table_1, smd = TRUE, showAllLevels = FALSE, explain = FALSE, pDigits = 4, formatOptions = list(big.mark = ","))) %>% tibble::rownames_to_column("Name")

    drops <- c("test")

    df <- df[ , !(names(df) %in% drops)]

    return(df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.94ddd233-4d76-448b-86ed-2c57f9253f1f"),
    m_ps=Input(rid="ri.foundry.main.dataset.ed7d0600-aece-4d29-8f2a-252033739f2f")
)
balance <- function(m_ps) {
    library(ggplot2)
    library(cobalt)
    data('m_ps', package= 'cobalt')
    set.cobalt.options(binary = "std")

    fig1 <- love.plot(MOLNUPIRAVIR ~ age_at_covid + I(sex) + I(race_ethnicity) + mildliver + cad + chf+ mi + htn + cancer + metcancer + ckd + dmuncomp + dmcomp + obesity + hiv + lungdisease+ immunocomp +  tobacco + vaccine_number, weights = m_ps$ipw, var.order = 'unadjusted', abs = TRUE, thresholds = c(m = 0.1), data = m_ps)

print(fig1)

return(NULL)
}

@transform_pandas(
    Output(rid="ri.vector.main.execute.e2b406a6-4372-43f7-9764-d75b723df234"),
    m_iptw_mvp=Input(rid="ri.foundry.main.dataset.3ad1e778-3c7e-44d9-8e75-8376cde14e36")
)
balance_mp <- function(m_iptw_mvp) {

    library(ggplot2)
    library(cobalt)
    data('m_iptw_mvp', package= 'cobalt')
    set.cobalt.options(binary = "std")

    fig2 <- love.plot(MOLNUPIRAVIR ~ age_at_covid + I(sex) + I(race_ethnicity) + mildliver + cad + chf+ mi + htn + cancer + metcancer + ckd + dmuncomp + dmcomp + obesity + hiv + lungdisease+ immunocomp +  tobacco + vaccine_number, weights = m_iptw_mvp$ipw, var.order = 'unadjusted', abs = TRUE, thresholds = c(m = 0.1), data = m_iptw_mvp)

print(fig2)

return(NULL)
}

@transform_pandas(
    Output(rid="ri.vector.main.execute.bd17bed8-4657-48a2-9c69-2cb09aa98053"),
    m_iptw_d3=Input(rid="ri.vector.main.execute.28994976-ae91-4f60-b019-136e6126ac0d")
)
d3_outcomes <- function(m_iptw_d3) {
    mort30w <- glm(mort_30 ~ MOLNUPIRAVIR, family = quasibinomial (link = logit), weights = ipw, data = m_iptw_d3)
    mort30u <- glm(mort_30 ~ MOLNUPIRAVIR, family = binomial (link = logit), data = m_iptw_d3)
    comp_w <- glm(d30_composite ~ MOLNUPIRAVIR, family = quasibinomial (link = logit), weights = ipw, data = m_iptw_d3)
    comp_u <- glm(d30_composite ~ MOLNUPIRAVIR, family = binomial (link = logit), data = m_iptw_d3)

    or_mort_w <- exp(cbind(OR = coef(mort30w), confint(mort30w)))
    or_mort_u <- exp(cbind(OR = coef(mort30u), confint(mort30u)))
    or_comp_w <- exp(cbind(OR = coef(comp_w), confint(comp_w)))
    or_comp_u <- exp(cbind(OR = coef(comp_u), confint(comp_u)))

    d3_outcomes <- data.frame(or_mort_w, or_mort_u, or_comp_w, or_comp_u)
    
    return(d3_outcomes)
}

@transform_pandas(
    Output(rid="ri.vector.main.execute.28994976-ae91-4f60-b019-136e6126ac0d"),
    molnu_std=Input(rid="ri.foundry.main.dataset.338c64a0-0274-4587-aca6-f6e45d2e196a")
)
m_iptw_d3 <- function(molnu_std) {
    library(dplyr)
    library(ggplot2)
    library(cobalt)

# Generate Propensity Score Model Based on Pre-specified DAG
    model_d3 <- glm(d3_MOLNUPIRAVIR ~ age_at_covid + I(sex) + I(race_ethnicity) + mildliver + cad + chf+ mi + htn + cancer + metcancer + ckd + dmuncomp + dmcomp + obesity + hiv + lungdisease+ immunocomp +  systemicsteroids + biologic + tobacco + vaccine_number, family=binomial (link=logit), data=molnu_std)

# Generate Propensity Scores, IPTW, and Weight Stabilized IPTW for Participants in Dataframe
    molnu_std$pr_score <- predict(model_d3, type= 'response')
    molnu_std$iptw <- 1/(predict(model_d3, type = 'response'))
    molnu_std$wt_stab <- mean(molnu_std$d3_MOLNUPIRAVIR)/(predict(model_d3, type = 'response'))
    molnu_std$ipw <- ifelse(molnu_std$d3_MOLNUPIRAVIR == 1, mean(molnu_std$d3_MOLNUPIRAVIR)/molnu_std$pr_score, (1-mean(molnu_std$d3_MOLNUPIRAVIR))/(1-molnu_std$pr_score))

# Calculated 1st and 99th percentil of PS and Trim Dataset to exclude those < 1st and > 99th percentiles

p05 <- quantile(molnu_std$pr_score, 0.05)
p95 <- quantile(molnu_std$pr_score, 0.95)

molnu_std <- subset(molnu_std, molnu_std$pr_score >p05 & molnu_std$pr_score <p95)

# Examine the region of common support
labs <- paste("Actual treatment:", c("Molnupiravir", "Supportive Care"))
molnu_std %>%
  mutate(d3_MOLNUPIRAVIR = ifelse(d3_MOLNUPIRAVIR==1, labs[1], labs[2]))

p <- ggplot(data=molnu_std, aes(x= ipw, colour = d3_MOLNUPIRAVIR)) +
        geom_density() +
        facet_wrap(~d3_MOLNUPIRAVIR) +
        xlab("Wt Stabilized IPTW of Molnupiravir Reciept")

    return(molnu_std)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.3ad1e778-3c7e-44d9-8e75-8376cde14e36"),
    molnu_paxlovid=Input(rid="ri.foundry.main.dataset.9967d4fe-c8da-4e64-a7d7-503047a4623f")
)
m_iptw_mvp <- function(molnu_paxlovid) {
library(dplyr)
library(ggplot2)

# Generate Propensity Score Model Based on Pre-specified DAG
    model2 <- glm(MOLNUPIRAVIR ~ age_at_covid + I(sex) + I(race_ethnicity) + mildliver + cad + chf+ mi + htn + cancer + metcancer + ckd + dmuncomp + dmcomp + obesity + hiv + lungdisease+ immunocomp + systemicsteroids + biologic + tobacco + vaccine_number + number_of_visits_before_covid + observation_period_before_covid, family=binomial (link=logit), data=molnu_paxlovid)

# Generate Propensity Scores, IPTW, and Weight Stabilized IPTW for Participants in Dataframe
    molnu_paxlovid$pr_score <- predict(model2, type= 'response')
    molnu_paxlovid$iptw <- 1/(predict(model2, type = 'response'))
    molnu_paxlovid$wt_stab <- mean(molnu_paxlovid$MOLNUPIRAVIR)/(predict(model2, type = 'response'))
    molnu_paxlovid$ipw <- ifelse(molnu_paxlovid$MOLNUPIRAVIR == 1, mean(molnu_paxlovid$MOLNUPIRAVIR)/molnu_paxlovid$pr_score, (1-mean(molnu_paxlovid$MOLNUPIRAVIR))/(1-molnu_paxlovid$pr_score))

# Calculated 1st and 99th percentil of PS and Trim Dataset to exclude those < 1st and > 99th percentiles

p01 <- quantile(molnu_paxlovid$pr_score, 0.01)
p99 <- quantile(molnu_paxlovid$pr_score, 0.99)

molnu_paxlovid <- subset(molnu_paxlovid, molnu_paxlovid$pr_score >p01 & molnu_paxlovid$pr_score <p99)

# Examine the region of common support
labs <- paste("Actual treatment:", c("Molnupiravir", "Paxlovid"))
molnu_paxlovid %>%
  mutate(MOLNUPIRAVIR = ifelse(MOLNUPIRAVIR==1, labs[1], labs[2]))

p <- ggplot(data=molnu_paxlovid, aes(x= ipw, colour = MOLNUPIRAVIR)) +
        geom_density() +
        facet_wrap(~MOLNUPIRAVIR) +
        xlab("Wt Stabilized IPTW of Molnupiravir Reciept")
        
    
# Assess for balance

    return(molnu_paxlovid)
    
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.ed7d0600-aece-4d29-8f2a-252033739f2f"),
    molnu_std=Input(rid="ri.foundry.main.dataset.338c64a0-0274-4587-aca6-f6e45d2e196a")
)
m_ps <- function(molnu_std) {
library(dplyr)
library(ggplot2)
library(cobalt)

# Generate Propensity Score Model Based on Pre-specified DAG
    model1 <- glm(MOLNUPIRAVIR ~ age_at_covid + I(sex) + I(race_ethnicity) + mildliver + cad + chf+ mi + htn + cancer + metcancer + ckd + dmuncomp + dmcomp + obesity + hiv + lungdisease+ immunocomp + systemicsteroids + biologic + tobacco + vaccine_number + number_of_visits_before_covid + observation_period_before_covid, family=binomial (link=logit), data=molnu_std)

# Generate Propensity Scores, IPTW, and Weight Stabilized IPTW for Participants in Dataframe
    molnu_std$pr_score <- predict(model1, type= 'response')
    molnu_std$iptw <- 1/(predict(model1, type = 'response'))
    molnu_std$wt_stab <- mean(molnu_std$MOLNUPIRAVIR)/(predict(model1, type = 'response'))
    molnu_std$ipw <- ifelse(molnu_std$MOLNUPIRAVIR == 1, mean(molnu_std$MOLNUPIRAVIR)/molnu_std$pr_score, (1-mean(molnu_std$MOLNUPIRAVIR))/(1-molnu_std$pr_score))

# Calculated 1st and 99th percentil of PS and Trim Dataset to exclude those < 1st and > 99th percentiles

p01 <- quantile(molnu_std$pr_score, 0.01)
p99 <- quantile(molnu_std$pr_score, 0.99)

molnu_std <- subset(molnu_std, molnu_std$pr_score >p01 & molnu_std$pr_score <p99)

# Examine the region of common support
labs <- paste("Actual treatment:", c("Molnupiravir", "Supportive Care"))
molnu_std %>%
  mutate(MOLNUPIRAVIR = ifelse(MOLNUPIRAVIR==1, labs[1], labs[2]))

p <- ggplot(data=molnu_std, aes(x= ipw, colour = MOLNUPIRAVIR)) +
        geom_density() +
        facet_wrap(~MOLNUPIRAVIR) +
        xlab("Wt Stabilized IPTW of Molnupiravir Reciept")

    return(molnu_std)
}

@transform_pandas(
    Output(rid="ri.vector.main.execute.7c0954a4-ddc9-48a8-ab20-059f0c4ec6c5"),
    pasc=Input(rid="ri.vector.main.execute.30a2b936-6f4f-4c83-8f3d-fef917f2fde5")
)
moln_pasc <- function(pasc) {
    library(dplyr)
    library(ggplot2)
    library(cobalt)

# Generate Propensity Score Model Based on Pre-specified DAG
    m_pasc <- glm(MOLNUPIRAVIR ~ age_at_covid + I(sex) + I(race_ethnicity) + mildliver + cad + chf+ mi + htn + cancer + metcancer + ckd + dmuncomp + dmcomp + obesity + hiv + lungdisease+ immunocomp + systemicsteroids + biologic + tobacco + vaccine_number + number_of_visits_before_covid + observation_period_before_covid, family=binomial (link=logit), data=pasc)

# Generate Propensity Scores, IPTW, and Weight Stabilized IPTW for Participants in Dataframe
    pasc$pr_score <- predict(m_pasc, type= 'response')
    pasc$iptw <- 1/(predict(m_pasc, type = 'response'))
    pasc$wt_stab <- mean(pasc$MOLNUPIRAVIR)/(predict(m_pasc, type = 'response'))
    pasc$ipw <- ifelse(pasc$MOLNUPIRAVIR == 1, mean(pasc$MOLNUPIRAVIR)/pasc$pr_score, (1-mean(pasc$MOLNUPIRAVIR))/(1-pasc$pr_score))

# Calculated 1st and 99th percentil of PS and Trim Dataset to exclude those < 1st and > 99th percentiles

p01 <- quantile(pasc$pr_score, 0.01)
p99 <- quantile(pasc$pr_score, 0.99)

pasc <- subset(pasc, pr_score >p01 & pr_score <p99)

return(pasc)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.f51f2a9d-a91c-4956-9104-118b27f2a8ba"),
    molnu_paxlovid=Input(rid="ri.foundry.main.dataset.9967d4fe-c8da-4e64-a7d7-503047a4623f")
)
# table1_ITM_other (7d2527b0-09d0-439c-9558-13083dae73f7): v10
molnu_pax_Table1 <- function(molnu_paxlovid) {

    library(tableone)
    library(dplyr)
    library(tibble)

    options(width=1000)

    all_variables_in_table1_rows <- c("sex", "age_at_covid", 'race_ethnicity', 'COVID_first_poslab_or_diagnosis_date', 'mildliver', 'rheum', 'chf', 'mi', 'cad', 'htn', 'ckd', 'cancer', 'metcancer', 'dmcomp', 'dmuncomp', 'obesity', 'lungdisease', 'transplant', 'hiv', 'immunocomp', 'systemicsteroids', 'biologic', 'tobacco', 'subuse', 'vaccine_number', 'number_of_visits_before_covid', 'number_of_visits_post_covid', 'observation_period_before_covid', 'observation_period_post_covid')

    categorical_variables_in_table1_rows <- c('sex','race_ethnicity', 'mildliver', 'rheum', 'chf', 'mi', 'cad', 'htn', 'ckd', 'cancer', 'metcancer', 'dmcomp', 'dmuncomp', 'obesity', 'lungdisease', 'transplant', 'hiv', 'immunocomp', 'systemicsteroids', 'biologic', 'tobacco', 'subuse')

    table_1 <- CreateTableOne(vars = all_variables_in_table1_rows, strata = "MOLNUPIRAVIR", data = molnu_paxlovid, factorVars = categorical_variables_in_table1_rows, addOverall = FALSE)

    df = as.data.frame(print(table_1, smd = TRUE, showAllLevels = FALSE, explain = FALSE, pDigits = 4, formatOptions = list(big.mark = ","))) %>% rownames_to_column("Name")

    drops <- c("test")

    df <- df[ , !(names(df) %in% drops)]

    return(df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.e42fae0e-904d-4366-b162-3f597d9bde20"),
    m_iptw_mvp=Input(rid="ri.foundry.main.dataset.3ad1e778-3c7e-44d9-8e75-8376cde14e36")
)
molnu_pax_Table2 <- function(m_iptw_mvp) {
    # table2

    library(tableone)
    library(dplyr)
    library(survey)

m_iptw_mvp_survey <- svydesign(ids=~1, data = m_iptw_mvp, weights = ~ipw)

    options(width=1000)

    all_variables_in_table1_rows <- c("sex", "age_at_covid", 'race_ethnicity', 'COVID_first_poslab_or_diagnosis_date', 'mildliver', 'rheum', 'chf', 'mi', 'cad', 'htn', 'ckd', 'cancer', 'metcancer', 'dmcomp', 'dmuncomp', 'obesity', 'lungdisease', 'transplant', 'hiv', 'immunocomp', 'systemicsteroids', 'biologic', 'tobacco', 'subuse', 'vaccine_number', 'number_of_visits_before_covid', 'number_of_visits_post_covid', 'observation_period_before_covid', 'observation_period_post_covid')

    categorical_variables_in_table1_rows <- c('sex','race_ethnicity', 'mildliver', 'rheum', 'chf', 'mi', 'cad', 'htn', 'ckd', 'cancer', 'metcancer', 'dmcomp', 'dmuncomp', 'obesity', 'lungdisease', 'transplant', 'hiv', 'immunocomp', 'systemicsteroids', 'biologic', 'tobacco', 'subuse')

    table_1 <- svyCreateTableOne(vars = all_variables_in_table1_rows, strata = "MOLNUPIRAVIR", data = m_iptw_mvp_survey, factorVars = categorical_variables_in_table1_rows, addOverall = FALSE)

    df = as.data.frame(print(table_1, smd = TRUE, showAllLevels = FALSE, explain = FALSE, pDigits = 4, formatOptions = list(big.mark = ","))) %>% tibble::rownames_to_column("Name")

    drops <- c("test")

    df <- df[ , !(names(df) %in% drops)]

    return(df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.9967d4fe-c8da-4e64-a7d7-503047a4623f"),
    iptw=Input(rid="ri.foundry.main.dataset.c3763cfc-0f9b-449e-b619-347a75f6b04d")
)
molnu_paxlovid <- function(iptw) {
    library(dplyr)
        local_df <- iptw %>%
        filter((PAXLOVID == 1 & MOLNUPIRAVIR == 0)| (MOLNUPIRAVIR == 1 & PAXLOVID ==0)) %>%
        filter(hosp180_pre_dx == 0) %>%
    return(local_df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.338c64a0-0274-4587-aca6-f6e45d2e196a"),
    iptw=Input(rid="ri.foundry.main.dataset.c3763cfc-0f9b-449e-b619-347a75f6b04d")
)
molnu_std <- function(iptw) {
    library(dplyr)
        local_df <- iptw %>%
        filter(PAXLOVID == 0) %>%
        filter(hosp180_pre_dx == 0) %>%
    return(local_df)
}

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.7f12d4b8-d1ed-46df-b73e-05cedba8a518"),
    m_iptw_mvp=Input(rid="ri.foundry.main.dataset.3ad1e778-3c7e-44d9-8e75-8376cde14e36")
)
outcome_m_p <- function(m_iptw_mvp) {
    
    mort30w <- glm(mort_30 ~ MOLNUPIRAVIR, family = quasibinomial (link = logit), weights = ipw, data = m_iptw_mvp)
    mort30u <- glm(mort_30 ~ MOLNUPIRAVIR, family = binomial (link = logit), data = m_iptw_mvp)
    comp_w <- glm(d30_composite ~ MOLNUPIRAVIR, family = quasibinomial (link = logit), weights = ipw, data = m_iptw_mvp)
    comp_u <- glm(d30_composite ~ MOLNUPIRAVIR, family = binomial (link = logit), data = m_iptw_mvp)

    or_mort_w <- exp(cbind(OR = coef(mort30w), confint(mort30w)))
    or_mort_u <- exp(cbind(OR = coef(mort30u), confint(mort30u)))
    or_comp_w <- exp(cbind(OR = coef(comp_w), confint(comp_w)))
    or_comp_u <- exp(cbind(OR = coef(comp_u), confint(comp_u)))

    outcomes <- data.frame(or_mort_w, or_mort_u, or_comp_w, or_comp_u)
    
    return(outcomes)
}

@transform_pandas(
    Output(rid="ri.vector.main.execute.6a2e6384-8fcd-427c-a8e7-cde35a3d2b99"),
    m_ps=Input(rid="ri.foundry.main.dataset.ed7d0600-aece-4d29-8f2a-252033739f2f")
)
outcome_m_s <- function(m_ps) {
    mort30w <- glm(mort_30 ~ MOLNUPIRAVIR, family = quasibinomial (link = logit), weights = ipw, data = m_ps)
    mort30wa <- glm(mort_30 ~ MOLNUPIRAVIR + age_at_covid, family = quasibinomial (link = logit), weights = ipw, data = m_ps)
    mort30u <- glm(mort_30 ~ MOLNUPIRAVIR, family = binomial (link = logit), data = m_ps)
    comp_w <- glm(d30_composite ~ MOLNUPIRAVIR, family = quasibinomial (link = logit), weights = ipw, data = m_ps)
    comp_wa <- glm(d30_composite ~ MOLNUPIRAVIR + age_at_covid, family = quasibinomial (link = logit), weights = ipw, data = m_ps)
    comp_u <- glm(d30_composite ~ MOLNUPIRAVIR, family = binomial (link = logit), data = m_ps)

    or_mort_w <- exp(cbind(OR = coef(mort30w), confint(mort30w)))
    or_mort_wa <- exp(cbind(OR = coef(mort30wa), confint(mort30wa)))
    or_mort_u <- exp(cbind(OR = coef(mort30u), confint(mort30u)))
    or_comp_w <- exp(cbind(OR = coef(comp_w), confint(comp_w)))
    or_comp_wa <- exp(cbind(OR = coef(comp_wa), confint(comp_wa)))
    or_comp_u <- exp(cbind(OR = coef(comp_u), confint(comp_u)))

    outcomes <- data.frame(or_mort_w, or_comp_w, or_mort_u, or_comp_u)
    
    return(outcomes)
}

@transform_pandas(
    Output(rid="ri.vector.main.execute.30a2b936-6f4f-4c83-8f3d-fef917f2fde5"),
    iptw=Input(rid="ri.foundry.main.dataset.c3763cfc-0f9b-449e-b619-347a75f6b04d")
)
pasc <- function(iptw) {
    library(dplyr)
    local_df <- iptw %>%
    filter(d180 == 1) %>%
return(local_df)
}

@transform_pandas(
    Output(rid="ri.vector.main.execute.073fe9be-901a-4834-9b50-a280ed7b6a70"),
    moln_pasc=Input(rid="ri.vector.main.execute.7c0954a4-ddc9-48a8-ab20-059f0c4ec6c5")
)
pasc_outcome <- function(moln_pasc) {
    pasc_w <- glm(pasc ~ MOLNUPIRAVIR, family = quasibinomial (link = logit), weights = ipw, data = moln_pasc)
    pasc_u <- glm(pasc ~ MOLNUPIRAVIR, family = binomial (link = logit), data = moln_pasc)

    or_pasc_w <- exp(cbind(OR = coef(pasc_w), confint(pasc_w)))
    or_pasc_u <- exp(cbind(OR = coef(pasc_u), confint(pasc_u)))

    outcomes <- data.frame(or_pasc_w, or_pasc_u)
    
    return(outcomes)
}

