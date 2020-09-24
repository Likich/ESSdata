library(foreign)
library(nFactors)
library(GPArotation)
library(lavaan)
library(psych)
library(semPlot)
library(kableExtra)
library(LittleHelpers)
library(DiagrammeR)
library(car)
library(devtools)
?semPlot
full_ess <- read.spss(
  file = "ESS4e04_5.sav",
  use.value.labels = F,
  use.missings = T,
  to.data.frame = T)

#First of all, let's perfrom the EFA on the other country, it was decided to choose Ukraine

UA.values <- full_ess[full_ess$cntry=="UA", c("ppltrst","pplfair","pplhlp",
                                              "trstlgl","trstplc","trstplt","trstep","trstun","trstprt","trstprl")]
#Scree plot
UA.values <-  na.omit(UA.values)
save(UA.values, file="UA.values.R")
UA.cov <- cov(UA.values)
UA.scree <-nScree(UA.cov, model="factors")
plot(UA.scree)

#We have supposedly 2 or 3 factors

twofactor <- fa(UA.values,nfactors = 2,rotate = "varimax",fm="minres")
print(twofactor$loadings,cutoff = 0.3)

#Let's check for 3 factors

threefactor <- fa(UA.values,nfactors = 3,rotate = "varimax",fm="minres")
print(threefactor$loadings,cutoff = 0.3)

#Which model is better?
formulaEFA <- as.formula("~ppltrst  + pplfair  + pplhlp + trstlgl + trstplc + trstplt + trstep + trstun + trstprt + trstprl")
EFA2 <- factanal(formulaEFA, factors=2, data=UA.values)
EFA3 <- factanal(formulaEFA, factors=3, data=UA.values)
chisq.difference <- EFA2$STATISTIC - EFA3$STATISTIC
df.difference <- EFA2$dof - EFA3$dof
1-pchisq(chisq.difference, df.difference)

#We should go for 3 factors model

RUfull.values <- full_ess[full_ess$cntry=="RU", ] 
RU.values <- full_ess[full_ess$cntry=="RU", c("ppltrst","pplfair","pplhlp",
                                              "trstlgl","trstplc","trstplt","trstep","trstun","trstprt","trstprl")]
#CFA
#Let's confirm there are 3 factors:

# Trust in people: ppltrst  + pplfair  + pplhlp
# Trust in government: trstlgl + trstplc + trstplt + trstprt + trstprl
# Trust in international policy: trstep + trstun 

#All indicators are measured in 11 point scale, so we don't need to run the standartization procedure


cfa1 <- cfa(model = "TrustInter =~ trstep + trstun;
                    TrustPeople =~ ppltrst  + pplfair  + pplhlp;
                    TrustGov =~ trstlgl + trstplc + trstplt + trstprt + trstprl",
            data=as.matrix(RU.values))

summary(cfa1, fit.measures=TRUE)

#RMSEA is high so we need to find another parameters to include

modindices(cfa1, sort = T)


cfa2 <- cfa(model = "TrustInter =~ trstep + trstun;
                    TrustPeople =~ ppltrst  + pplfair  + pplhlp;
                    TrustGov =~ trstlgl + trstplc + trstplt + trstprt + trstprl
            trstplt ~~ trstprt",
            data=as.matrix(RU.values))

summary(cfa2, fit.measures=TRUE)

#RMSEA is high so we need to find another parameters to include

modindices(cfa2, sort = T)

cfa3 <- cfa(model = "TrustInter =~ trstep + trstun;
                    TrustPeople =~ ppltrst  + pplfair  + pplhlp;
                    TrustGov =~ trstlgl + trstplc + trstplt + trstprt + trstprl;
            trstplt ~~ trstprt;
             trstlgl ~~ trstplc",
            data=as.matrix(RU.values))

summary(cfa3, fit.measures=TRUE)

#All tests are satisfied

lavTestLRT(cfa1, cfa2, cfa3)
#Chi-squared is significantly lower with each next model

semPaths(cfa3, whatLabels="est", layout = "tree" , rotation=2, esize = 3, shapeLat = "ellipse", nCharNodes=5)

fit.index<-data.frame(fitmeasures(cfa1),
                      fitmeasures(cfa2),
                      fitmeasures(cfa3)
)
names(fit.index) <- c("Трехфакторная модель", "Трехфакторная модель с ковариацией  2 остатков", "Трехфакторная с ковариацией 4 остатков") 
fit.index<- fit.index[c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "bic"),] 
tb2 <- knitr::kable(fit.index,  digits=3)
kableExtra::kable_styling(tb2)

      #Second order factor model

SecondOrder <- sem("TrustInter =~ trstep + trstun;
                TrustPeople =~ ppltrst  + pplfair  + pplhlp;
                TrustGov =~ trstlgl + trstplc + trstplt + trstprt + trstprl;
                trstplt ~~ trstprt;
                trstlgl ~~ trstplc;
                Trust =~ TrustInter + TrustPeople + TrustGov", RU.values)
summary(SecondOrder, fit.measures = TRUE)
semPaths(SecondOrder, whatLabels="est", layout = "tree" , label.cex = 1.5, 
         edge.label.cex = 1, edge.color = "darkgreen", style = "lisrel", 
         rotation=2, esize = 3, shapeLat = "ellipse", nCharNodes=5)

#The impact of forst order factors into the second order one is significant
#Let's build a model of structured equation checking for indirect effect

RUfull.values$impfree.r <- car::Recode(RUfull.values$impfree, "1=6; 2=5; 3=4; 4=3; 5=2; 6=1; else=NA")

semindirect <- sem("TrustInter =~ trstep + trstun;
             TrustPeople =~ ppltrst  + pplfair  + pplhlp;
             TrustGov =~ trstlgl + trstplc + trstplt + trstprt + trstprl;
             trstplt ~~ trstprt;
             trstlgl ~~ trstplc;
             TrustInter ~ stfeco;
             TrustPeople ~  stfeco;
             TrustGov ~ stfeco + I*impfree.r + TrustPeople;
             impfree.r ~ ED*eduyrs;
             indirect_ED := I*ED;", RUfull.values)


semlatentindirect <- sem("TrustInter =~ trstep + trstun;
             TrustPeople =~ ppltrst  + pplfair  + pplhlp;
             TrustGov =~ trstlgl + trstplc + trstplt + trstprt + trstprl;
             trstplt ~~ trstprt;
             trstlgl ~~ trstplc;
             TrustInter ~ stfeco;
             TrustPeople ~  STFE*stfeco;
             TrustGov ~ I*impfree.r;
             impfree.r ~ ED*eduyrs;
             stfeco ~ TG*TrustGov;
             indirect_STFE := TG*STFE;
             indirect_ED := I*ED;", RUfull.values);

newCFA <- cfa(semindirect, RUfull.values)
summary(newCFA, fit.measures=TRUE)
newCFA1 <- cfa(semlatentindirect, RUfull.values)
summary(newCFA1, fit.measures=TRUE)
semPaths(newCFA, whatLabels="est", layout = "tree", label.cex = 1.5, 
         edge.label.cex = 1, edge.color = "darkgreen", style = "ram", 
         rotation=2, esize = 3, shapeLat = "ellipse", nCharNodes=5)

fit.index<-data.frame(fitmeasures(cfa1),
                      fitmeasures(cfa2),
                      fitmeasures(cfa3),
                      fitmeasures(newCFA),
                      fitmeasures (newCFA1)
)

names(fit.index) <- c("Three factor model", "Three factor model with covariation", "Three factor model with 2 covariations", "Model with structure component", "Model with indirect effect between latent variables") 
fit.index<- fit.index[c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "bic"),] 
tb2 <- knitr::kable(fit.index,  digits=3)
kableExtra::kable_styling(tb2)





