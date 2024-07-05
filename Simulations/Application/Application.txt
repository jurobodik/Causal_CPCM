# can also be downloaded from https://www.kaggle.com/datasets/grosvenpaul/family-income-and-expenditure

data <- read.csv("Family Income and Expenditure.csv")

data=data[data$`Total Household Income`>quantile(data$`Total Household Income`, 0.1), ] #above poverty line
data=data[data$`Total Number of Family members` == 1, ] #family of size one


X1 = data$`Total Household Income`
X2 = data$`Total Food Expenditure`
X3 = data$`Alcoholic Beverages Expenditure` #GAM estimation does not handle well too many ties. +5 is there only to ensure that it does not go into negative numbers
X = data.frame(X1,X2,X3)

plot(X)
hist(X1[X1<400000])
hist(X2)
hist(X3)

CPCM_graph_estimate(data.frame(X1, X2), family_of_distributions = 2)
CPCM_graph_estimate(data.frame(X1, X3), family_of_distributions = 2)
CPCM_graph_estimate(data.frame(X2, X3), family_of_distributions = 2)


#joint estimate of the graph. This takes a while (10 minutes?) since we only implemented brute force algorithm that goes through all DAGs. Also, independence test between d>2 variables is very slow
#result is the incidence matrix (X_i --> X_j if position (i,j) is equal to 1)
set.seed(1) 
CPCM_graph_estimate( X , family_of_distributions = 2)$result

#>result
#>     [,1] [,2] [,3]
#>[1,]    0    1    0
#>[2,]    0    0    0
#>[3,]    0    1    0








