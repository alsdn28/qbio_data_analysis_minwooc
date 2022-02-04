# exercise 1.1
data("attenu")
na_rows <- is.na(attenu$station)
attenu_cleaned <- attenu[na_rows==FALSE, ]
head(attenu_cleaned)
dim(attenu_cleaned)
  #166 x 5

# exercise 1.2
data("Theoph")
Theoph_2 <- Theoph
str(Theoph_2)
median(Theoph_2$Dose)
  # median = 4.53
Theoph_2$Dose_Class <- ifelse(Theoph_2$Dose>=4.53, "high", "low")
head(Theoph_2)
dim(Theoph_2)
  # 132 x 6

# exercise 1.3
library(readr)
library(here)
path <- here("week3_R", "week3_homework", "starbucks.csv")
starbucks <- read_csv(path)
View(starbucks)
is_row_empty <- is.na(starbucks)
# Need 6 na to contain no data
is_row_empty <- ifelse(rowSums(is_row_empty)==6, TRUE, FALSE)
nrow(starbucks)
length(is_row_empty)
# both 177 
starbucks_cleaned <- starbucks[is_row_empty==FALSE, ]
attach(starbucks_cleaned)
plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, xlab="Carbohydrates(g)", ylab="Calories")
max(starbucks_cleaned$Calories)
highcal <- starbucks_cleaned$Calories==430
which(highcal==TRUE)
starbucks_cleaned[65, "Drink"]
# Starbucks Signature Hot Chocolate
max(starbucks_cleaned$Fat)
starbucks_cleaned$is_highest_fat <- ifelse(starbucks_cleaned$Fat==26, TRUE, FALSE)
plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, xlab="Carbohydrates(g)", ylab="Calories", 
     col=as.factor(starbucks_cleaned$is_highest_fat))

# Exercise 1.4
library(readr)
path_2 <- here("week3_R", "week3_homework", "Batting.csv")
Batting <- read_csv(path_2)
View(Batting)
Batting$goodplayers <- ifelse(Batting$HR>=3, TRUE, FALSE)
sum(Batting$goodplayers)
# 26155 players with >3 homeruns
homeruns <- plot(Batting$yearID, Batting$HR, xlab="Year", ylab= "Homeruns")
angels <- ifelse(Batting$teamID=="LAA", TRUE, FALSE)
angels <- Batting[angels,]
angels$goodplayers <- ifelse(angels$HR>=3, TRUE, FALSE)
sum(angels$goodplayers)
# 259 players with >3 homeruns
teams <- ifelse(Batting$teamID=="ATL"|Batting$teamID=="PIT", TRUE, FALSE)
teams <- Batting[teams,]
na_homeruns <- is.na(teams$HR)
teams <- teams[na_homeruns==FALSE, ]
plot(teams$yearID, teams$HR, xlab="Players", ylab="Homeruns", col=ifelse(teams$teamID=="ATL", "red", "black"))

#Exercise 1.5
easy_plot <- function(x, y, color_data) {
  med <- median(color_data)
  levels <- ifelse(color_data>=med, "high", "low")
  levels <- factor(levels)
  print(med)
  plot(x=x, y=y, col=levels, pch=20)
  cor.test(x,y)
}
easy_plot(starbucks_cleaned$Fat, starbucks_cleaned$Carb, starbucks_cleaned$Calories)
easy_plot(Batting$G, Batting$AB, Batting$R)
  
# Exercise 2.1
data("iris")
# data set describes the iris flower with 150 data points and 5 variables
# Exercise 2.2
# sepal length, width and petal length and width are continuous variables
# species is categorical 
# Exercise 2.3
hist(x=iris$Sepal.Length)
hist(x=iris$Sepal.Width)
hist(x=iris$Petal.Length)
hist(x=iris$Petal.Width)
# Exercise 2.4
mean_sepal_width <- mean(iris$Sepal.Width)
iris_copy <- iris
iris$thickness <- ifelse(iris$Sepal.Width>=mean_sepal_width, TRUE, FALSE)
boxplot(iris$Sepal.Width~iris$thickness)
# Exercise 2.5
pairs(iris_copy[, 1:4], col=iris_copy$Species)
#black dataset most unique
