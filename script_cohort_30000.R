########## Objectif du script est de faire des stats sur les donnees d'une cohorte ##########
########## de 30000 personnes ayant eu le COVID, cohort issu de l'IHU Mediterranne ##########

library(data.table)
library(tidyverse)
library(DataExplorer)
library(easystats)
library(emmeans)
library(tidymodels)
library(here)

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


## realisation de la courbe ROC et autre tests pour juger de la qualite

# prediction theorique
dt_cohort = dt_cohort[!is.na(DEATH),]
pred_prob = predict(model_3, newdata = dt_cohort, type = "response")

# realisation d'une table specifique pour la courbe ROC
df_ROC = data.table(truth = dt_cohort$DEATH |> as.factor(),
                    pred_0 = 1 - pred_prob,
                    pred_1 = pred_prob,
                    prediction = as.factor(ifelse(pred_prob > 0.5, 1, 0) ) )

# application de la fonction calculant l'AUC + la courbe
roc_auc(df_ROC, truth, pred_0 )
roc_curve(df_ROC, truth, pred_0 ) |> 
  ggplot(aes(1 - specificity, sensitivity )) + 
  geom_path() + 
  geom_abline(linetype = 4) +
  annotate(geom = "label", 0.33, 0.75,
           label = paste("AUC = ",roc_auc(df_ROC, truth, pred_0 )[3] |> round(3))) +
  theme_bw()
# Le score est très bon, mais cela cache quelque chose !

# Créer l'objet PRROC (l'idee est d'avoir un score de 1)
pr_auc(df_ROC, truth, pred_0)
pr_curve(df_ROC, truth, pred_0 ) |> 
  ggplot(aes(recall, precision)) +
  annotate(geom = "label", 0.33, 0.998,
           label = paste("AUC = ",pr_auc(df_ROC, truth, pred_0)[3] |> round(3))) +
  geom_path() +
  theme_bw()
# encore une fois, c'est trop beau pour etre vrai
  
# Matrice de confusion
df_ROC |> 
  conf_mat(truth = truth, estimate = prediction) |>
  summary(event_level = "second")
# kap = 0.017 signifie que le modele fait a peine mieux que le hazard (qui est de 0)
# sens = 0.00909 signifie que ~1% des vrais positifs sont trouves (c'est catastrophique)
# spec = 1 signifie que 100% des negatifs sont trouves mais c'est parce que le modele ne trouve que ca
# ppv = 0.17 signifie que 17% des positifs du model le sont pour de vrai (c'est tres mauvais)
# npv = 0.993 signifie que 99.3% des negatifs du model le sont mais comme iel predit que des negatifs...
# bal_accuracy = 0.504, un model random aurait 0.5 donc ...
## CL : ce modele est usless car il ne peut pas predire les positifs du au fait qu'il n'y en a 
## presque pas dans les donnees de base, il vas donc falloir faire quelque chose !!!


##### Changer la valeur du seuil afin d'augmenter la detection des 1 #####
## 2 methodes : la simulation et la formule

## methode 1 : simualtion pour determiner la valeur de seuil qui maximise le Fscore
fscore = NA
for(i in seq(0, 5000, 1)){
  
  fscore[i] = df_ROC |> 
    _[, prediction_bis := factor(ifelse(pred_prob > i / 10000, 1, 0), levels = c(0,1))] |> 
    conf_mat(truth = truth, estimate = prediction_bis) |>
    summary(event_level = "second") |>
    _[, ".estimate"][13,]
  }
fscore[sapply(fscore, is.null)] = NA

df_fscore = data.table(threshold = seq_along(fscore) / 10000, f1_score = unlist(fscore)) |> 
  _[order(-f1_score),]
# d'apres ma simulation, la meilleur valeur de score etant 0.0673
# CL : il faut que je considere comme "1" des que la prediction me donne une valeur de 0.065 et plus
# on est bien loin du classique 0.5

# confirmation du Fscore au nouveau seuil
df_ROC |> 
  _[, prediction_bis := factor(ifelse(pred_prob > 0.0673, 1, 0), levels = c(0,1))] |> 
  conf_mat(truth = truth, estimate = prediction_bis) |>
  summary(event_level = "second") |>
  _[, ".estimate"][13,]


## methode 2 : formule mathematique

# creation d'une varialbe "index" base sur l'ordre decroissant des valeurs "pred_1"
df_ROC_bis = df_ROC[order(-pred_1), .(truth, pred_1)] |> 
  _[, index := 1:.N][!is.na(pred_1),]

# nb de positif en realite (110 dans mes data)
npos = sum(df_ROC_bis$truth == "1")

# Somme cumule des vrais positifs (110 a la fin)
df_ROC_bis[, ':='(
  tp_cum = cumsum(truth == "1") )] 

# Calcul de rappel et precision, puis du Fscore qui necessite c'est 2 parametres
df_ROC_bis[, ':='( rappel = tp_cum / npos,
                    precision = tp_cum / (1:.N))] |> 
  _[, fscore := 2 * (precision * rappel) / (precision + rappel)]
df_ROC_bis[order(-fscore),] #le meilleur seuil serait 0.067899 proche de 0.0673

## AU final les 2 methodes donne un resultat tres proche soucis, la performance n'est pas assez bonne car
df_ROC_bis[, prediction := as.factor(ifelse(pred_1 > 0.067899, 1, 0) ) ] |> 
  conf_mat(truth = truth, estimate = prediction) |>
  summary(event_level = "second")

# kap = 0.193 signifie que le modele fait mieux que le hazard (qui est de 0) mais bon...
# sens = 0.364 signifie que ~36% des vrais positifs sont trouves (c'est trop faible)
# spec = 0.984 signifie que 98% des negatifs sont trouves (c'est pas le plsu utile)
# ppv = 0.139  signifie que 14% des positifs du model le sont pour de vrai (c'est mauvais mais pas grave)
# npv = 0.996 signifie que 99.6% des negatifs du model le sont mais comme iel predit que des negatifs...
# bal_accuracy = 0.674, un model random aurait 0.5 donc c'est pas terrible
## CL : ce modele est faible car il ne predit pas bien les positifs qui sont les plus important a predire


##### Augmenter artificiellement le nombre de 1 #####
## on se doute bien que cela peut engendrer des problemes mais on verra 

# creation de la table avec les "1" augmente
dt_cohort_oversampled = copy(dt_cohort)[, .(DEATH, SEX, AGE, HCQ, AZ, IVM, OBESITY, VACCINATION, HBP,
                                          DIABETE, ASTHMA, CANCER, IMMUNODEFICIENCY, ChronicCardiacDiseases,
                                          AutoImmuneDiseases, COPD)] |> na.omit()

# creation d'une table contenant tout les "1" duplique "146 fois" pour match le nb de "0"
dt_positif = dt_cohort_oversampled[DEATH == 1,] |> 
  _[rep(1:.N, each = nrow(dt_cohort_oversampled) / .N),]

# fusion de cette table remplis de "1" et ma table d'origine quasi remplis de "0"
dt_cohort_oversampled = rbind(dt_cohort_oversampled, dt_positif)

# partion du data en data train et test
split = dt_cohort_oversampled |> 
  initial_split(prop = 0.8)

dt_train = training(split)
dt_test = testing(split)

# modelisation sur les donnees train
m_train = glm(DEATH ~ SEX + AGE + HCQ + IVM + AZ + OBESITY + VACCINATION + 
               DIABETE + CANCER + IMMUNODEFICIENCY + COPD, 
             data = dt_train, 
             family = binomial(link = "logit"))
model_parameters(m_train) ; model_parameters(m_train) |> plot()


## test de la qualite du model
pred_prob_train = predict(m_train, newdata = dt_train, type = "response")

# realisation d'une table specifique pour la courbe ROC
df_ROC_train = data.table(truth = dt_train$DEATH |> as.factor(),
                    pred_0 = 1 - pred_prob_train,
                    pred_1 = pred_prob_train,
                    prediction = as.factor(ifelse(pred_prob_train > 0.5, 1, 0) ) )

# application de la fonction calculant l'AUC + la courbe
roc_auc(df_ROC_train, truth, pred_0 )
roc_curve(df_ROC_train, truth, pred_0 ) |> 
  ggplot(aes(1 - specificity, sensitivity )) + 
  geom_path() + 
  geom_abline(linetype = 4) +
  annotate(geom = "label", 0.33, 0.75,
           label = paste("AUC = ",roc_auc(df_ROC_train, truth, pred_0 )[3] |> round(3))) +
  theme_bw()
# Le score est bon

# Matrice de confusion
df_ROC_train |> 
  conf_mat(truth = truth, estimate = prediction) |>
  summary(event_level = "second")
# kap = 0.705 tres bon
# f_meas = 0852 tres bon 

## rtester le model sur les donnees test

pred_prob_test = predict(m_train, newdata = dt_test, type = "response")

# realisation d'une table specifique pour la courbe ROC
df_ROC_test = data.table(truth = dt_test$DEATH |> as.factor(),
                          pred_0 = 1 - pred_prob_test,
                          pred_1 = pred_prob_test,
                          prediction = as.factor(ifelse(pred_prob_test > 0.5, 1, 0) ) )

# application de la fonction calculant l'AUC + la courbe
roc_auc(df_ROC_test, truth, pred_0 )
roc_curve(df_ROC_test, truth, pred_0 ) |> 
  ggplot(aes(1 - specificity, sensitivity )) + 
  geom_path() + 
  geom_abline(linetype = 4) +
  annotate(geom = "label", 0.33, 0.75,
           label = paste("AUC = ",roc_auc(df_ROC_test, truth, pred_0 )[3] |> round(3))) +
  theme_bw()
# Le score est bon

# Matrice de confusion
df_ROC_test |> 
  conf_mat(truth = truth, estimate = prediction) |>
  summary(event_level = "second")
# kap = 0.705 tres bon
# f_meas = 0852 tres bon





## test le model sur les donnees d'origine

pred_prob_origin = predict(m_train, newdata = dt_cohort, type = "response")

# realisation d'une table specifique pour la courbe ROC
df_ROC_origin = data.table(truth = dt_cohort$DEATH |> as.factor(),
                         pred_0 = 1 - pred_prob_origin,
                         pred_1 = pred_prob_origin,
                         prediction = as.factor(ifelse(pred_prob_origin > 0.5, 1, 0) ) )

# application de la fonction calculant l'AUC + la courbe
roc_auc(df_ROC_origin, truth, pred_0 )
roc_curve(df_ROC_origin, truth, pred_0 ) |> 
  ggplot(aes(1 - specificity, sensitivity )) + 
  geom_path() + 
  geom_abline(linetype = 4) +
  annotate(geom = "label", 0.33, 0.75,
           label = paste("AUC = ",roc_auc(df_ROC_origin, truth, pred_0 )[3] |> round(3))) +
  theme_bw()
# Le score est bon

# Matrice de confusion
df_ROC_origin |> 
  conf_mat(truth = truth, estimate = prediction) |>
  summary(event_level = "second")
## a priori c'est mieux car il ne rate que 17 mort sur les 110 mais je veux 
## qu'il n'en rate aucun, quand je fais varier le seuil, malheureusement
## il met des 1 partout pour detecter les 110 vrai mort