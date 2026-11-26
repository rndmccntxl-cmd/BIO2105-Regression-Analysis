# ==============================================================================
# BIO2105: Biostatistics — Chapter 9: Regression Analysis (Formative Assessment 12)
# Topic: Correlation, Simple Linear Regression, and Residual Diagnostics
# Datasets: trees (Part I) and women (Part II)
# ==============================================================================

# Reset any stuck graphics devices before plotting
graphics.off()

# ==============================================================================
# PART I: Black Cherry Trees Dataset (trees)
# ==============================================================================

# ------------------------------------------------------------------------------
# Task 1: Variable Roles & Data Exploration
# ------------------------------------------------------------------------------
data(trees)

head(trees)
str(trees)
summary(trees[, c("Girth", "Volume")])

# ------------------------------------------------------------------------------
# Task 2: Scatter Plot with Fitted Regression Line
# ------------------------------------------------------------------------------
tree_model <- lm(Volume ~ Girth, data = trees)

par(mar = c(4, 4, 2.5, 1))
plot(
  trees$Girth, trees$Volume,
  xlab = "Tree Girth (inches)",
  ylab = "Timber Volume (cubic feet)",
  main = "Black Cherry Trees: Girth vs. Timber Volume (n = 31)",
  pch = 19, col = "steelblue", cex = 1.2
)
grid(col = "gray80")
abline(tree_model, col = "darkorange", lwd = 2)
legend(
  "topleft",
  legend = c("Observed Trees", "OLS Regression Line"),
  pch = c(19, NA),
  lty = c(NA, 1),
  col = c("steelblue", "darkorange"),
  lwd = c(NA, 2)
)
par(mar = c(5.1, 4.1, 4.1, 2.1))

# ------------------------------------------------------------------------------
# Task 3: Pearson's Correlation Coefficient (r) — Built-in vs. Manual Formula
# ------------------------------------------------------------------------------
r_builtin <- cor(trees$Girth, trees$Volume)

n     <- nrow(trees)
x_bar <- mean(trees$Girth)
y_bar <- mean(trees$Volume)
s_x   <- sd(trees$Girth)
s_y   <- sd(trees$Volume)

r_manual <- (1 / (n - 1)) * sum(((trees$Girth - x_bar) / s_x) * ((trees$Volume - y_bar) / s_y))

cat(sprintf("Built-in correlation r = %.6f\n", r_builtin))
cat(sprintf("Manual correlation r   = %.6f\n", r_manual))
cor.test(trees$Girth, trees$Volume)

# ------------------------------------------------------------------------------
# Task 4: Invariance Under Linear Transformation
# ------------------------------------------------------------------------------
girth_cm <- 2.54 * trees$Girth
volume_L <- 28.3168 * trees$Volume
r_metric <- cor(girth_cm, volume_L)

cat(sprintf("Original r (in, ft^3): %.6f\n", r_builtin))
cat(sprintf("Metric r   (cm, L)   : %.6f\n", r_metric))
cat(sprintf("Absolute difference  : %.2e\n", abs(r_builtin - r_metric)))

# ------------------------------------------------------------------------------
# Task 5: Simple Linear Regression Estimation & Centroid Property
# ------------------------------------------------------------------------------
coeffs <- coef(tree_model)
b0 <- coeffs[1]
b1 <- coeffs[2]
print(coeffs)

b1_from_r <- r_builtin * (s_y / s_x)
cat(sprintf("Slope from lm()        : %.6f\n", b1))
cat(sprintf("Slope from r*(s_y/s_x) : %.6f\n", b1_from_r))

y_hat_at_xbar <- b0 + b1 * x_bar
cat(sprintf("Sample mean y_bar      : %.6f\n", y_bar))
cat(sprintf("Predicted y at x_bar   : %.6f\n", y_hat_at_xbar))

# ------------------------------------------------------------------------------
# Task 6: Coefficient of Determination (R^2)
# ------------------------------------------------------------------------------
r_squared <- summary(tree_model)$r.squared
cat(sprintf("R-squared from summary(): %.6f\n", r_squared))
cat(sprintf("r^2 from correlation    : %.6f\n", r_builtin^2))

# ------------------------------------------------------------------------------
# Task 7: Inference for Slope & Confidence vs. Prediction Intervals
# ------------------------------------------------------------------------------
summary(tree_model)
confint(tree_model, level = 0.95)

new_trees <- data.frame(Girth = c(10, 14, 18))

ci_mean <- predict(tree_model, newdata = new_trees, interval = "confidence", level = 0.95)
cat("\n95% Confidence Intervals (Mean Response):\n")
print(round(ci_mean, 3))

pi_indiv <- predict(tree_model, newdata = new_trees, interval = "prediction", level = 0.95)
cat("\n95% Prediction Intervals (Individual Response):\n")
print(round(pi_indiv, 3))


# ==============================================================================
# PART II: American Women Height & Weight Dataset (women)
# ==============================================================================

# ------------------------------------------------------------------------------
# Task 8: Model Fitting & Scatter Plot
# ------------------------------------------------------------------------------
data(women)
head(women)

women_model <- lm(weight ~ height, data = women)
summary(women_model)

par(mar = c(4, 4, 2.5, 1))
plot(
  women$height, women$weight,
  xlab = "Height (inches)",
  ylab = "Weight (lbs)",
  main = "American Women (Ages 30-39): Height vs. Weight (n = 15)",
  pch = 19, col = "steelblue", cex = 1.2
)
grid(col = "gray80")
abline(women_model, col = "darkorange", lwd = 2)
legend(
  "topleft",
  legend = c("Observed Data", "OLS Regression Line"),
  pch = c(19, NA),
  lty = c(NA, 1),
  col = c("steelblue", "darkorange"),
  lwd = c(NA, 2)
)
par(mar = c(5.1, 4.1, 4.1, 2.1))

# ------------------------------------------------------------------------------
# Task 9: ANOVA Table & Sum of Squares Decomposition
# ------------------------------------------------------------------------------
if (!exists("women_model")) {
  data(women)
  women_model <- lm(weight ~ height, data = women)
}

anova(women_model)

y_obs  <- women$weight
y_hat  <- fitted(women_model)
res    <- residuals(women_model)
y_mean <- mean(y_obs)

ss_total <- sum((y_obs - y_mean)^2)
ss_reg   <- sum((y_hat - y_mean)^2)
ss_error <- sum(res^2)

cat(sprintf("SS_Total             : %.4f\n", ss_total))
cat(sprintf("SS_Reg               : %.4f\n", ss_reg))
cat(sprintf("SS_Error             : %.4f\n", ss_error))
cat(sprintf("SS_Reg + SS_Error    : %.4f\n", ss_reg + ss_error))
cat(sprintf("R^2 = SS_Reg/SS_Total: %.6f\n", ss_reg / ss_total))

# ------------------------------------------------------------------------------
# Task 10: Residual Diagnostics & Assumption Checking
# ------------------------------------------------------------------------------
if (!exists("women_model")) {
  data(women)
  women_model <- lm(weight ~ height, data = women)
}

par(mfrow = c(2, 2), mar = c(3.8, 3.8, 2.2, 1))
plot(women_model)

par(mfrow = c(1, 1), mar = c(4, 4, 2.5, 1))

# Dedicated Residuals vs. Fitted Values plot with zero reference line
plot(
  fitted(women_model), residuals(women_model),
  xlab = "Fitted Values (lbs)",
  ylab = "Residuals (lbs)",
  main = "Residuals vs. Fitted Values (women_model)",
  pch = 19, col = "steelblue", cex = 1.2
)
grid(col = "gray80")
abline(h = 0, col = "red", lty = 2, lwd = 2)

# Dedicated Normal Q-Q plot of residuals
qqnorm(
  residuals(women_model),
  main = "Normal Q-Q Plot of Residuals (women_model)",
  pch = 19, col = "steelblue", cex = 1.2
)
qqline(residuals(women_model), col = "red", lty = 2, lwd = 2)
grid(col = "gray80")

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))