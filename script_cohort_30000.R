########## Objectif du script est de faire des stats sur les donnees d'une cohorte ##########
########## de 30000 personnes ayant eu le COVID, cohort issu de l'IHU Mediterranne ##########

library(data.table)
library(tidyverse)
library(DataExplorer)
library(easystats)
library(emmeans)

my_directory = setwd("C:/Users/mcout/Desktop/Statistique/stat_script_bazard/covid_19/covid_cohort_30000")

##### import de la table #####
dt_cohort = fread("Database_Cohort_30423_COVID-19_IHU_032923.txt") |> 
  _[, AGE := as.factor(AGE)]

##### relalisation du model sur la proba de mourir #####
plot_missing(dt_cohort) # beaucoup de valeur manquante pour certaines variable importante

## model simple de la proba de mourir en fonction du sexe
model_1 = glm(DEATH ~ SEX, 
               data = dt_cohort, 
               family = binomial(link = "logit"))
model_parameters(model_1) # il y a bien un effet SEXE sur la proba de mourir
model_performance(model_1) 

# on calcul la difference entre les 2 sexes
model_1_result = emmeans(model_1, specs = pairwise ~ SEX, type = "response")
model_1_result ; plot(model_1_result)
# le sexe 1 (male) a 3% de proba de mourir 
# le sexe 2 (femelle) a 1.8% de proba de mourir
# le sexe 1 (male) a 1.7 fois plus de proba de mourir que le sexe 2 (femelle)

## complexification du model pour integrer d'autre variable
model_2 = glm(DEATH ~ SEX + AGE + HCQ + AZ + IVM + OBESITY + VACCINATION + HBP + 
                 DIABETE + ASTHMA + CANCER + IMMUNODEFICIENCY + ChronicCardiacDiseases +
                 AutoImmuneDiseases + COPD, 
               data = dt_cohort, , 
               family = binomial(link = "logit"))
model_parameters(model_2) |> plot()

## on va enlever les variables clairement non significative
model_3 = glm(DEATH ~ SEX + AGE + HCQ + IVM + AZ + OBESITY + VACCINATION + 
                DIABETE + CANCER + IMMUNODEFICIENCY + COPD, 
              data = dt_cohort, , 
              family = binomial(link = "logit"))
model_parameters(model_3) ; model_parameters(model_3) |> plot()

AIC(model_1, model_2, model_3) # model 1 le meilleur sans doute possible

# on calcul la difference entre les differents level des differentes variables
model_3_sexe = emmeans(model_3, specs = pairwise ~ SEX, type = "response")
model_3_sexe ; plot(model_3_sexe) # risque * 1.9 pour les hommes 

model_3_age = emmeans(model_3, specs = pairwise ~ AGE, type = "response")
model_3_age ; plot(model_3_age) # la classe d'age 1 a 350 * moins de risque que la 4

model_3_vaccination = emmeans(model_3, specs = pairwise ~ VACCINATION, type = "response")
model_3_vaccination ; plot(model_3_vaccination) # les vaccines on 3 * moins de risque


