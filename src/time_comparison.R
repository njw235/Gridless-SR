library(microbenchmark)
i=1

compare_runtime<- function(M, n_alphas, runs=100){
  results <- data.frame(method = character(), time_ns = numeric(), value = numeric(),  run_id = integer())
  for (i in 1:runs) {
    
    set.seed(1)
    rho = runif(1)
    set.seed(12345+i)
    # Simulate an AR(1) chain
    chainParams = list(type="AR",  M = M, rho = rho)
    ch = generateChain(chainParams)
    x = ch$x
    r = autocov(x)
    delta = tune_delta(x,nSplits = 5)$delta*0.8
    alphaGrid = makeGrid(nX = n_alphas, upper_threshold = 1 - delta, cm = FALSE, scale = "equidist")
    s_alpha = sqrt((1 + alphaGrid^2)/(1 - alphaGrid^2))
    
    # execution time in nano-seconds
    t1 <- microbenchmark(m1 <- SR1(ch_mh$x,alphaGrid = alphaGrid, Xtr_approx = T), times = 1)$time
    t2 <- microbenchmark(m2 <- SR1(ch_mh$x,alphaGrid = alphaGrid, Xtr_approx = F), times = 1)$time
    
    # function values
    val1<-momentLS::L2diff_L2Moment(r = r,support = m1$support,weights = m1$weights)
    val2<-momentLS::L2diff_L2Moment(r = r,support = m2$support,weights = m2$weights)
    
    results <- rbind(results,
                     data.frame(method = "method1", time_ns = t1, value = val1,run_id=i),
                     data.frame(method = "method2", time_ns = t2, value = val2,run_id=i)
    )
  }
  return(results)
}

comp1 = compare_runtime(M = 1e5,n_alphas = 1000, runs = 100)

library(dplyr)
comp1 %>% data.frame() %>% 
  group_by(method) %>% 
  mutate(time_ms = time_ns / 10^6) %>% 
  select(method,time_ms,value) %>% 
  summarise(
    mean_time_ms = mean(time_ms),
    se_time_ms = sd(time_ms) / sqrt(n()),
    mean_value = mean(value),
    se_value = sd(value) / sqrt(n())
  )
