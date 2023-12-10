## -----------------------------------------------------------------------------
my.sample <- function(x, size, replace = TRUE, prob = rep(1/length(x), length(x))) {
  if (!replace && size > length(x)) {
    stop("Cannot take a sample larger than the population when 'replace = FALSE'")
  }
  
  cum_prob <- cumsum(prob)
  if (tail(cum_prob, 1) != 1) {
    cum_prob <- cum_prob / tail(cum_prob, 1)
  }
  
  results <- numeric(size)
  
  for (i in 1:size) {
    u <- runif(1)
    selected <- which.max(cum_prob > u)
    results[i] <- x[selected]
    if (!replace) {
      prob[selected] <- 0
      cum_prob <- cumsum(prob)
      cum_prob <- cum_prob / tail(cum_prob, 1)
    }
  }
  
  return(results)
}

# 测试
print(my.sample(1:10, 15))
print(my.sample(1:10, 9, replace=FALSE))
print(my.sample(1:10, 5, prob=c(0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05)))


## -----------------------------------------------------------------------------
generate_laplace <- function(n) {
  u <- runif(n)
  x <- ifelse(u < 0.5, log(2*u), -log(2*(1-u)))
  return(x)
}

samples <- generate_laplace(1000)

## -----------------------------------------------------------------------------
hist(samples, prob=TRUE, breaks=30, col="gray", main="Laplace Distribution", 
     xlab="x", xlim=c(-4, 4), ylim=c(0, 0.5), border="black")
curve(0.5 * exp(-abs(x)), add=TRUE, col="red", lwd=2)
legend("topright", legend=c("Generated Samples", "True Density"), 
       lty=1, col=c("gray", "red"))

## -----------------------------------------------------------------------------
generate_beta_ar <- function(a, b, n) {
  c <- 10
  samples <- numeric(n)
  count <- 1
  
  while(count <= n) {
    y <- runif(1) # proposal sample
    u <- runif(1)
    if(u <= dbeta(y, a, b) / (c * dunif(y))) {
      samples[count] <- y
      count <- count + 1
    }
  }
  
  return(samples)
}

samples <- generate_beta_ar(3, 2, 1000)


hist(samples, prob=TRUE, breaks=30, col="gray", main="Beta(3,2) Distribution", 
     xlab="x", xlim=c(0, 1), ylim=c(0, 2), border="black")
curve(dbeta(x, 3, 2), add=TRUE, col="red", lwd=2)
legend("topleft", legend=c("Generated Samples", "True Density"), 
       lty=1, col=c("gray", "red"))

## -----------------------------------------------------------------------------
generate_epanechnikov <- function(n) {
  samples <- numeric(n)
  for (i in 1:n) {
    U1 <- runif(1, -1, 1)
    U2 <- runif(1, -1, 1)
    U3 <- runif(1, -1, 1)
    
    if (abs(U3) >= abs(U2) && abs(U3) >= abs(U1)) {
      samples[i] <- U2
    } else {
      samples[i] <- U3
    }
  }
  return(samples)
}

# Simulate a large sample
samples <- generate_epanechnikov(10000)

# Plotting
hist(samples, prob=TRUE, breaks=30, col="gray", main="Rescaled Epanechnikov Kernel", 
     xlab="x", xlim=c(-1, 1), border="black")
curve((3/4)*(1-x^2) * (abs(x) <= 1), add=TRUE, col="red", lwd=2)
legend("topright", legend=c("Simulated Samples", "True Density"), 
       lty=1, col=c("gray", "red"), x.intersp=0.8, y.intersp=0.8, cex=0.8)

## -----------------------------------------------------------------------------
set.seed(12345)

simulate_pi_variance <- function(l) {
  d <- 1
  m <- 1e6
  pihat_values <- numeric(100)
  
  for (i in 1:100) {
    X <- runif(m, 0, d/2)
    Y <- runif(m, 0, pi/2)
    pihat <- 2*l/d/mean(l/2*sin(Y) > X)
    pihat_values[i] <- pihat
  }
  
  var_pihat <- var(pihat_values)
  
  list(l_value = l, variance = var_pihat)
}


result <- simulate_pi_variance(1)
print(paste("rho =", result$l_value, ", var =", result$variance))

result <- simulate_pi_variance(0.5)
print(paste("rho =", result$l_value, ", var =", result$variance))

result <- simulate_pi_variance(0.1)
print(paste("rho =", result$l_value, ", var =", result$variance))



## -----------------------------------------------------------------------------
compute_mean <- function(n, fun) {
  replicate(n, fun())
}

m <- 10000
n_rep <- 1000

# Direct Monte Carlo
mc_fun <- function() {
  mean(exp(runif(m)))
}
mc_results <- compute_mean(n_rep, mc_fun)

# Antithetic Variates
anti_fun <- function() {
  u <- runif(m / 2)
  v <- 1 - u
  mean((exp(u) + exp(v)) / 2)
}
anti_results <- compute_mean(n_rep, anti_fun)

v1 <- var(mc_results)
v2 <- var(anti_results)

list(
  mc_mean = mean(mc_results),
  anti_mean = mean(anti_results),
  mc_variance = v1,
  anti_variance = v2,
  variance_ratio = (v1 - v2) / v1
)


## -----------------------------------------------------------------------------
x <- seq(1, 10, 0.01)
y <- x^2 * exp(-x^2 / 2) / sqrt(2 * pi)
lambda <- 2  
p <- 3  
f1 <- lambda * exp(-lambda * (x - 1))
f2 <- dgamma(x-1, 3/2, 2)
plot(x, y, type="l", ylim=c(0, 1))
lines(x, f1, lty=2)
lines(x, f2, lty=3)
legend("topright", inset=0.02, legend=c("g(x)", "f_1", "f_2"), lty=1:3)

## -----------------------------------------------------------------------------
plot(x, y / f1, type="l", lty=2, ylab="ratio")
lines(x, y / f2, lty=3)
legend("topright", inset=0.02, legend=c("f_1", "f_2"), lty=2:3)

## -----------------------------------------------------------------------------
N <- 10000


I1 <- replicate(1000, expr = {
x <- rexp(N, rate = 2) + 1
f <- dexp(x-1, rate = 2)
g <- (x^2 / sqrt(2*pi))*exp(-x^2 / 2)
hx_gx <- g / f
mean(hx_gx)
})

I2 <- replicate(1000, expr = {
x <- rgamma(N, 3/2, 2) + 1  
f <- dgamma(x-1, 3/2, 2)
g <- (x^2 / sqrt(2*pi))*exp(-x^2 / 2)
hx_gx <- g / f
mean(hx_gx)
})


c(mean(I1), mean(I2))

c(var(I1), var(I2))

## -----------------------------------------------------------------------------
set.seed(123)
M <- 10000
k <- 5
m <- M / k
si <- numeric(k)
v <- numeric(k)
g <- function(x) exp(-x)/(1+x^2)
f <- function(x) (k/(1-exp(-1))) * exp(-x)

for (j in 1:k) {
    u <- runif(m, (j-1)/k, j/k)
    x <- -log(1-(1-exp(-1))*u)
    fg <- g(x)/f(x)
    si[j] <- mean(fg)
    v[j] <- var(fg)
}

sum(si)
mean(v)
sqrt(mean(v))

## -----------------------------------------------------------------------------
set.seed(123)
k <- 1
m <- M / k
si <- numeric(k)
v <- numeric(k)

for (j in 1:k) {
    u <- runif(m, (j-1)/k, j/k)
    x <- -log(1-(1-exp(-1))*u)
    fg <- g(x)/f(x)
    si[j] <- mean(fg)
    v[j] <- var(fg)
}

sum(si)
mean(v)
sqrt(mean(v))

## -----------------------------------------------------------------------------
set.seed(123)  

n <- 20
num_samples <- 10000
true_mean <- 2  
coverage_count <- 0

for(i in 1:num_samples) {
  sample_data <- rchisq(n, df=2)
  sample_mean <- mean(sample_data)
  sample_sd <- sd(sample_data)
  
  t_value <- qt(0.975, df=n-1)
  margin_error <- t_value * sample_sd / sqrt(n)
  
  lower_bound <- sample_mean - margin_error
  upper_bound <- sample_mean + margin_error
  
  if (lower_bound <= true_mean && upper_bound >= true_mean) {
    coverage_count <- coverage_count + 1
  }
}

coverage_prob <- coverage_count / num_samples
coverage_prob


## -----------------------------------------------------------------------------
set.seed(123)
n <- 20
alpha <- 0.05
num_samples <- 10000
reject_counts <- c(chi2=0, uniform=0, exponential=0)

for(i in 1:num_samples) {
  # Chi-squared(1) distribution
  sample_data <- rchisq(n, df=1)
  test_result <- t.test(sample_data, mu=1)  # Mean of chi-squared(1) is 1
  if(test_result$p.value < alpha) {
    reject_counts["chi2"] <- reject_counts["chi2"] + 1
  }
  
  # Uniform(0,2) distribution
  sample_data <- runif(n, 0, 2)
  test_result <- t.test(sample_data, mu=1)  # Mean of Uniform(0,2) is 1
  if(test_result$p.value < alpha) {
    reject_counts["uniform"] <- reject_counts["uniform"] + 1
  }
  
  # Exponential(1) distribution
  sample_data <- rexp(n, rate=1)
  test_result <- t.test(sample_data, mu=1)  # Mean of Exponential(1) is 1
  if(test_result$p.value < alpha) {
    reject_counts["exponential"] <- reject_counts["exponential"] + 1
  }
}

empirical_type1_errors <- reject_counts / num_samples
empirical_type1_errors

## -----------------------------------------------------------------------------
m <- 1000
n_simulations <- 1000
alpha <- 0.1

results <- data.frame(Method = character(),
                      FWER = numeric(),
                      FDR = numeric(),
                      TPR = numeric())

for (method in c("bonferroni", "BH")) {
  fwer_sum <- 0
  fdr_sum <- 0
  tpr_sum <- 0
  
  for (i in 1:n_simulations) {
    # 生成p值
    p_under_null <- runif(m * 0.95)
    p_under_alt <- rbeta(m * 0.05, 0.1, 1)
    p_values <- c(p_under_null, p_under_alt)
    
    # 应用校正
    adjusted_p_values <- p.adjust(p_values, method = method)
    
    # 计算错误率
    rejections <- adjusted_p_values < alpha
    false_positives <- sum(rejections[1:(m*0.95)])
    true_positives <- sum(rejections[(m*0.95+1):m])
    
    fwer_sum <- fwer_sum + (false_positives > 0)
    fdr_sum <- fdr_sum + false_positives/sum(rejections)
    tpr_sum <- tpr_sum + true_positives/(m*0.05)
  }
  
  results <- rbind(results, data.frame(Method = method,
                                       FWER = fwer_sum/n_simulations,
                                       FDR = fdr_sum/n_simulations,
                                       TPR = tpr_sum/n_simulations))
}

results

## -----------------------------------------------------------------------------
# 设定参数
lambda <- 2
sample_sizes <- c(5, 10, 20)
B <- 1000
m <- 1000

# 初始化结果矩阵
results <- matrix(0, length(sample_sizes), 4)
colnames(results) <- c("Theoretical Bias", "Bootstrap Mean Bias", "Theoretical SE", "Bootstrap Mean SE")

# 循环不同的样本大小
for (j in 1:length(sample_sizes)) {
  n <- sample_sizes[j]
  
  bootstrap_biases <- numeric(m)
  bootstrap_SEs <- numeric(m)
  
  for (i in 1:m) {
    # 从指数分布中抽取样本
    sample <- rexp(n, rate=lambda)
    sample_mean <- mean(sample)
    lambda_hat <- 1/sample_mean
    
    # 对lambda_hat进行bootstrap
    bootstrap_estimates <- numeric(B)
    for (b in 1:B) {
      resample <- sample(sample, replace=TRUE)
      bootstrap_estimates[b] <- 1/mean(resample)
    }
    
    bootstrap_bias <- mean(bootstrap_estimates) - lambda_hat
    bootstrap_SE <- sd(bootstrap_estimates)
    
    bootstrap_biases[i] <- bootstrap_bias
    bootstrap_SEs[i] <- bootstrap_SE
  }
  
  results[j, ] <- c(lambda/(n-1), mean(bootstrap_biases), lambda*n/((n-1)*sqrt(n-2)), mean(bootstrap_SEs))
}

rownames(results) <- paste("n=", sample_sizes)
results


## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)

data(law)
LSAT <- law$LSAT
GPA <- law$GPA

cor_function <- function(data, indices) {
  d <- data[indices, ]  
  return(cor(d$LSAT, d$GPA))
}

set.seed(12345)  
results <- boot(data=law, statistic=cor_function, R=10000)

conf_interval <- boot.ci(results, type="bca")

print(conf_interval)


## -----------------------------------------------------------------------------
library(boot)

x <- aircondit[,1]

meant <- function(x, i) {
  return(mean(as.matrix(x[i])))
}

b <- boot(x, statistic=meant, R=2000)

print(b)

ci_result <- boot.ci(b, type=c("norm", "basic", "perc", "bca"))

print(ci_result)

hist(b$t, prob=TRUE, main="")
points(b$t0, 0, cex=2, pch=16)



## -----------------------------------------------------------------------------
library(bootstrap)
attach(scor)

x <- as.matrix(scor)
n <- nrow(x)

theta.jack <- numeric(n)

lambda <- eigen(cov(x))$values
theta.hat <- max(lambda/sum(lambda))


for (i in 1:n) {
  y <- x[-i,]
  s <- cov(y)
  lambda <- eigen(s)$values
  theta.jack[i] <- max(lambda/sum(lambda))
}


bias.jack <- (n-1) * (mean(theta.jack) - theta.hat)
se.jack <- sqrt((n-1)/n * sum((theta.jack - mean(theta.jack))^2))

result <- list(est = theta.hat, bias = bias.jack, se = se.jack)
result

## -----------------------------------------------------------------------------
compareDistributions <- function(sampleA, sampleB, permCount = 199) {
  countA <- length(sampleA)
  countB <- length(sampleB)
  mergedSamples <- c(sampleA, sampleB)
  totalCount <- countA + countB
  ecdfA <- numeric(totalCount)
  ecdfB <- numeric(totalCount)
  
  # 构建经验分布函数
  for (index in 1:totalCount) {
    ecdfA[index] <- mean(mergedSamples[index] <= sampleA)
    ecdfB[index] <- mean(mergedSamples[index] <= sampleB)
  }
  
  # 计算Cramer-von Mises标准统计量
  cvmOriginal <- ((countA * countB) / totalCount^2) * sum((ecdfA - ecdfB)^2)
  
  # 生成经过排列的统计量分布
  cvmPermuted <- sapply(1:permCount, function(i) {
    permIndices <- sample(totalCount)
    permutedSamples <- mergedSamples[permIndices]
    samplePermutedA <- permutedSamples[1:countA]
    samplePermutedB <- permutedSamples[(countA + 1):totalCount]
    ecdfPermutedA <- sapply(1:totalCount, function(j) mean(permutedSamples[j] <= samplePermutedA))
    ecdfPermutedB <- sapply(1:totalCount, function(j) mean(permutedSamples[j] <= samplePermutedB))
    ((countA * countB) / totalCount^2) * sum((ecdfPermutedA - ecdfPermutedB)^2)
  })
  
  # 确定p值
  pVal <- mean(cvmPermuted >= cvmOriginal)
  
  return(list(cvmStat = cvmOriginal, pVal = pVal))
}

# 使用chickwts数据集
data(chickwts)
feedTypes <- split(chickwts$weight, chickwts$feed)

# 比较不同饲料类型
cvmComparison1 <- compareDistributions(feedTypes$soybean, feedTypes$linseed)
cvmComparison2 <- compareDistributions(feedTypes$sunflower, feedTypes$linseed)

# 结果输出
cvmComparison1
cvmComparison2


## -----------------------------------------------------------------------------
# 异常值边界计数函数
calcOutlierBounds <- function(data1, data2) {
  d1Adjusted <- data1 - mean(data1)
  d2Adjusted <- data2 - mean(data2)
  outliersD1 <- sum(d1Adjusted > max(d2Adjusted)) + sum(d1Adjusted < min(d2Adjusted))
  outliersD2 <- sum(d2Adjusted > max(d1Adjusted)) + sum(d2Adjusted < min(d1Adjusted))
  max(outliersD1, outliersD2)
}

# 基于最大异常值统计量的排列检验
permTestVarEquality <- function(sample1, sample2, numPermutations = 199) {
  combined <- c(sample1, sample2)
  size1 <- length(sample1)
  sizeTotal <- length(combined)
  permStats <- replicate(numPermutations, {
    permutedIndices <- sample(sizeTotal)
    calcOutlierBounds(combined[permutedIndices[1:size1]], combined[permutedIndices[(size1 + 1):sizeTotal]])
  })
  actualStat <- calcOutlierBounds(sample1, sample2)
  combinedStats <- c(permStats, actualStat)
  distribution <- table(combinedStats) / (numPermutations + 1)
  list(estimation = actualStat, pValue = mean(combinedStats >= actualStat), frequency = distribution, cumulative = cumsum(distribution))
}

# 测试等方差情形
set.seed(100)
group1Size <- 20
group2Size <- 40
avg <- 0
stdDeviation1 <- stdDeviation2 <- 1
group1 <- rnorm(group1Size, avg, stdDeviation1)
group2 <- rnorm(group2Size, avg, stdDeviation2)
permTestVarEquality(group1, group2)

# 测试不等方差情形
set.seed(100)
stdDeviation1 <- 1
stdDeviation2 <- 2
group1 <- rnorm(group1Size, avg, stdDeviation1)
group2 <- rnorm(group2Size, avg, stdDeviation2)
permTestVarEquality(group1, group2)


## -----------------------------------------------------------------------------
library(ggplot2)

# 设计函数
calculate_a <- function(N, b1, b2, b3, f0) {
    # 生成数据
    X1 <- rpois(N, 1)
    X2 <- rexp(N, 1)
    X3 <- rbinom(N, 1, 0.5)
    
    # 使用数值方法寻找 a 的值
    find_a <- function(a) {
        p <- exp(a + b1*X1 + b2*X2 + b3*X3) / (1 + exp(a + b1*X1 + b2*X2 + b3*X3))
        mean(p) - f0
    }
    
    # 使用优化函数来寻找合适的 a 值
    optimize(find_a, c(-10, 10))$minimum
}

# 调用函数并绘图
results <- data.frame()
f0_values <- c(0.1, 0.01, 0.001, 0.0001)

for (f0 in f0_values) {
    a_value <- calculate_a(10^6, 0, 1, -1, f0)
    results <- rbind(results, data.frame(f0 = f0, a = a_value))
}

results$log_neg_f0 <- -log(results$f0)

# 绘图
ggplot(results, aes(x = log_neg_f0, y = a)) +
    geom_point() +
    geom_line() +
    xlab("-log(f0)") +
    ylab("a")


## -----------------------------------------------------------------------------
# 设置样本大小和烧入样本大小
num_samples <- 5000
burn_in <- 1000

# 初始化样本矩阵
samples <- matrix(0, num_samples, 2)

# 设置相关系数和均值
correlation <- 0.9
mean1 <- 0
mean2 <- 0

# 设置标准差
std_dev1 <- 1
std_dev2 <- 1

# 计算条件标准差
cond_std1 <- sqrt(1 - correlation^2) * std_dev1
cond_std2 <- sqrt(1 - correlation^2) * std_dev2

# 初始化第一个样本
samples[1,] <- c(mean1, mean2)

# 吉布斯抽样
for (i in 2:num_samples) {
  prev_y <- samples[i - 1, 2]
  cond_mean1 <- mean1 + correlation * (prev_y - mean2) * std_dev1/std_dev2
  samples[i, 1] <- rnorm(1, cond_mean1, cond_std1)
  
  curr_x <- samples[i, 1]
  cond_mean2 <- mean2 + correlation * (curr_x - mean1) * std_dev2/std_dev1
  samples[i, 2] <- rnorm(1, cond_mean2, cond_std2)
}

# 去除烧入样本
final_samples <- samples[(burn_in + 1):num_samples, ]
X <- final_samples[, 1]
Y <- final_samples[, 2]

# 线性回归模型
linear_model <- lm(Y ~ X)

# 输出模型结果
summary(linear_model)

## -----------------------------------------------------------------------------
plot(X, Y, cex = 0.25)
abline(h = 0, v = 0)

## -----------------------------------------------------------------------------
# 残差图
plot(linear_model$fitted.values, linear_model$residuals, cex = 0.25)
abline(h = 0)

# QQ图
qqnorm(linear_model$residuals, cex = 0.25)
qqline(linear_model$residuals)

## -----------------------------------------------------------------------------
# 载入所需的包
library(coda)

# 定义Rayleigh分布的密度函数
Rayleigh_density <- function(z, beta) { 
  if (z < 0) 
    return(0) 
  stopifnot(beta > 0) 
  return((z/beta^2) * exp(-z^2/(2 * beta^2))) 
}

# 使用Metropolis-Hastings算法生成Rayleigh分布的马尔可夫链
Rayleigh_MH_chain <- function(beta, chain_length, initial_value) { 
  chain <- numeric(chain_length) 
  chain[1] <- initial_value 
  random_values <- runif(chain_length) 
  for (j in 2:chain_length) { 
    current <- chain[j - 1] 
    proposed <- rchisq(1, df = current) 
    acceptance_ratio <- Rayleigh_density(proposed, beta) * dchisq(current, df = proposed) /
                        Rayleigh_density(current, beta) * dchisq(proposed, df = current) 
    if (random_values[j] <= acceptance_ratio) 
      chain[j] <- proposed 
    else 
      chain[j] <- current
  }
  return(chain) 
}

# 设置参数
beta_value <- 4
initial_values <- c(1/beta_value^2, 1/beta_value, beta_value^2, beta_value^3)
chain_count <- 4
chain_length <- 2000

# 创建一个空矩阵来存储所有链
Markov_chains <- matrix(0, nrow = chain_count, ncol = chain_length)

# 生成多条马尔可夫链
for (j in 1:chain_count) 
  Markov_chains[j, ] <- Rayleigh_MH_chain(beta_value, chain_length, initial_values[j])

# 计算每一条链的累积均值
cumulative_means <- t(apply(Markov_chains, 1, cumsum))
for (j in 1:nrow(cumulative_means)) {
  cumulative_means[j, ] <- cumulative_means[j, ] / (1:ncol(cumulative_means))
}

# 将每条链转换为mcmc对象并进行Gelman-Rubin诊断
mcmc_objects <- lapply(1:nrow(cumulative_means), function(j) mcmc(cumulative_means[j, ]))
mcmc_chain_list <- mcmc.list(mcmc_objects)
gelman_rubin_diagnostic <- gelman.diag(mcmc_chain_list)
gelman_rubin_diagnostic



## -----------------------------------------------------------------------------
gelman.plot(mcmc_chain_list, col = c(1, 1))


## -----------------------------------------------------------------------------
library(stats4)

# 观测值区间
intervals <- matrix(c(11, 12, 8, 9, 27, 28, 13, 14, 16, 17, 0, 1, 23, 24, 10, 11, 24, 25, 2, 3), ncol = 2, byrow = TRUE)

# 定义负对数似然函数
neg_log_likelihood <- function(lambda, intervals) {
  likelihood <- sapply(1:nrow(intervals), function(i) {
    ui <- intervals[i, 1]
    vi <- intervals[i, 2]
    p <- integrate(function(x) lambda * exp(-lambda * x), lower = ui, upper = vi, rel.tol = 1e-8, abs.tol = 1e-10)$value
    max(p, .Machine$double.eps)  # 避免非有限值
  })
  -sum(log(likelihood))
}

# 使用 optim 来最大化似然函数
optimize_lambda <- function(intervals) {
  result <- optim(par = 1, fn = neg_log_likelihood, intervals = intervals, method = "L-BFGS-B", lower = 1e-7, upper = 100)
  return(result$par)
}

# 运行优化
lambda_est <- optimize_lambda(intervals)
print(lambda_est)



## -----------------------------------------------------------------------------
library(stats4)

# 观测值区间
intervals <- matrix(c(11, 12, 8, 9, 27, 28, 13, 14, 16, 17, 0, 1, 23, 24, 10, 11, 24, 25, 2, 3), ncol = 2, byrow = TRUE)

# E-step: 计算条件期望
E_step <- function(lambda, intervals) {
  sapply(1:nrow(intervals), function(i) {
    ui <- intervals[i, 1]
    vi <- intervals[i, 2]
    (integrate(function(x) x * lambda * exp(-lambda * x), lower = ui, upper = vi)$value) / 
    (integrate(function(x) lambda * exp(-lambda * x), lower = ui, upper = vi)$value)
  })
}

# M-step: 更新 lambda
M_step <- function(lambda, intervals) {
  expectations <- E_step(lambda, intervals)
  n <- nrow(intervals)
  lambda_new <- n / sum(expectations)
  return(lambda_new)
}

# EM 算法
EM_algorithm <- function(intervals, tolerance = 1e-6, max_iter = 1000) {
  lambda <- 1  # 初始 lambda
  for (i in 1:max_iter) {
    lambda_new <- M_step(lambda, intervals)
    if (abs(lambda_new - lambda) < tolerance) {
      break
    }
    lambda <- lambda_new
  }
  return(lambda)
}

# 运行 EM 算法
lambda_est <- EM_algorithm(intervals)
print(lambda_est)


## -----------------------------------------------------------------------------
library(lpSolve)

solve.game <- function(A) {
  A_min <- min(A)
  A_max <- max(A)
  A <- (A - A_min) / A_max

  m <- nrow(A)
  n <- ncol(A)

  objective <- c(rep(1, m), 0)

  constr <- rbind(cbind(-diag(m), rep(1, m)), 
                  cbind(A, rep(0, n)))
  rhs <- c(rep(0, m), rep(1, n))

  direction <- c(rep("<=", m + n))

  result <- lp("max", objective, constr, direction, rhs, 
               all.int = FALSE)

  # 提取结果
  strategy <- result$solution[1:m]
  value <- result$objval * A_max + A_min

  list(strategy = strategy, value = value)
}

# 定义原始游戏矩阵A
A <- matrix(c(0, -2, -2, 3, 0, 0, 4, 0, 0, 2, 0, 0, 0, 
              -3, -3, 4, 0, 0, 2, 0, 0, 3, 0, 0, 0, -4, -4, -3, 
              0, -3, 0, 4, 0, 0, 5, 0, 0, 3, 0, -4, 0, -4, 0, 5, 
              0, 0, 3, 0, 0, 4, 0, -5, 0, -5, -4, -4, 0, 0, 0, 
              5, 0, 0, 6, 0, 0, 4, -5, -5, 0, 0, 0, 6, 0, 0, 4, 
              0, 0, 5, -6, -6, 0), 9, 9, byrow = TRUE)

B <- A + 2

result_B <- solve.game(B)

result_B$value
result_B$strategy


## -----------------------------------------------------------------------------
# 创建一个有0行的数据框
df_zero_rows <- data.frame(x = numeric(0), y = character(0))

# 创建一个有0列的数据框
df_zero_columns <- data.frame(row.names = 1:10)


## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

df <- data.frame(
  numbers1 = c(1, 2, 3, 4, 5),
  numbers2 = c(10, 20, 30, 40, 50)
)

# 应用于数据框的每一列
df_scaled <- data.frame(lapply(df, scale01))
df_scaled

# 仅应用于数值型列
df_scaled_numeric <- data.frame(
  numbers1 = c(1, 2, 3, 4, 5),
  numbers2 = c(10, 20, 30, 40, 50),
  factors = factor(c("a", "b", "a", "b", "c")),
  characters = c("A", "B", "C", "D", "E")
)

numeric_cols <- sapply(df, is.numeric)  # 识别数值型列
df_scaled_numeric[numeric_cols] <- lapply(df[numeric_cols], scale01)
df_scaled_numeric

## -----------------------------------------------------------------------------
# 创建一个只包含数值列的数据框
df_numeric <- data.frame(
  column1 = rnorm(100),  # 生成100个正态分布的随机数
  column2 = runif(100),  # 生成100个均匀分布的随机数
  column3 = rpois(100, lambda = 2)  # 生成100个泊松分布的随机数
)

std_devs <- vapply(df_numeric, sd, numeric(1))
std_devs

## -----------------------------------------------------------------------------
# 创建一个包含数值型和非数值型列的数据框
df_mixed <- data.frame(
  numeric_column = rnorm(100),  # 数值列
  factor_column = factor(sample(letters[1:3], 100, replace = TRUE)),  # 因子列
  logical_column = sample(c(TRUE, FALSE), 100, replace = TRUE),  # 逻辑列
  character_column = sample(letters[1:3], 100, replace = TRUE)  # 字符列
)

# 识别数值型列
is_numeric_col <- vapply(df_mixed, is.numeric, logical(1))

# 对数值型列计算标准差
std_devs_numeric <- vapply(df_mixed[is_numeric_col], sd, numeric(1))
std_devs_numeric


## -----------------------------------------------------------------------------
gibbsSamplerR <- function(n, a, b, N) {
  x <- numeric(N)
  y <- numeric(N)

  x[1] <- round(n / 2)  # 初始值
  y[1] <- rbeta(1, x[1] + a, n - x[1] + b)

  for (i in 2:N) {
    x[i] <- rbinom(1, n, y[i - 1])
    y[i] <- rbeta(1, x[i] + a, n - x[i] + b)
  }

  return(data.frame(x, y))
}


## -----------------------------------------------------------------------------
library(Rcpp)

cppFunction('
NumericMatrix gibbsSamplerRcpp(int n, double a, double b, int N) {
  NumericMatrix samples(N, 2);
  double x = round(n / 2.0);
  double y = R::rbeta(x + a, n - x + b);

  samples(0, 0) = x;
  samples(0, 1) = y;

  for (int i = 1; i < N; ++i) {
    x = R::rbinom(n, y);
    y = R::rbeta(x + a, n - x + b);
    samples(i, 0) = x;
    samples(i, 1) = y;
  }

  return samples;
}')



## -----------------------------------------------------------------------------
# install.packages("microbenchmark")
library(microbenchmark)

# 设置参数
n <- 10
a <- 2
b <- 2
N <- 1000

# 比较
benchmark_result <- microbenchmark(
  R = gibbsSamplerR(n, a, b, N),
  Rcpp = gibbsSamplerRcpp(n, a, b, N),
  times = 10
)

print(benchmark_result)


