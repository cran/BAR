next_allocation_rate_BAR <- function(n, success_count, tot_num, power_c = "n/2N",
                     lower_bound = .05, control_arm = "", seed = 100){
  set.seed(seed)
  if (!control_arm %in% c("", "fixed")){
    stop("Please make sure the argument of control_arm is specified as default or fixed")
  }
  prior_a <- .5
  prior_b <- .5
  shapeBetas <- matrix(rep(c(prior_a, prior_b), length(success_count)), nrow = length(success_count), byrow = TRUE)
  for(i in 1:length(success_count)){
    shapeBetas[i, 1] <- shapeBetas[i, 1] + success_count[i]
    shapeBetas[i, 2] <- shapeBetas[i, 2] + n[i] - success_count[i]
  }
  for(i in 1:length(success_count)){
    assign(paste("a", i, sep = ""), rbeta(1000,shapeBetas[i,1],shapeBetas[i,2]))
  }
  pmaxfunction <- paste0("pmax(",paste(paste0("a", 1:length(success_count)), collapse =","),")")
  for(i in 1:length(success_count)){
    assign(paste("c", i, sep = ""), eval(parse(text=paste0("length(which(a",i,"==", pmaxfunction,"))/1000"))))
  }
  if (control_arm == "fixed"){
    previous_rem_sum <- 1-c1
    c1 <- 1/length(success_count)
    after_rem_sum <- 1-c1
    for(xx in 2:length(success_count)){
      assign(paste("c", xx, sep = ""), eval(parse(text = paste0("c", xx, "/previous_rem_sum*after_rem_sum"))))
    }
  }
  if(power_c == "n/2N"){
    lambda <- sum(n)/(2*tot_num)
  }
  else{
    lambda <- power_c
  }
  lambdafunction1 <- paste(paste0("c", 1:length(success_count),"^lambda"), collapse = "+")
  lambdafunction2 <- paste(paste0("c", 2:length(success_count),"^lambda"), collapse = "+")
  if(control_arm == "fixed"){
    rem_sum <- 1-c1
    r1 <- c1
    for(z in 2:length(success_count)){
      assign(paste("r", z, sep = ""), eval(parse(text=paste0("(c",z,"^lambda)/(",lambdafunction2,")*rem_sum"))))
    }

  }
  else{
    for(z in 1:length(success_count)){
      assign(paste("r", z, sep = ""), eval(parse(text=paste0("(c",z,"^lambda)/(",lambdafunction1,")"))))
    }
  }
  assign("r", eval(parse(text=paste("c(",paste(paste0("r", 1:length(success_count)), collapse = ","),")"))))
  all_i <- c(1:length(r))
  max_i <- which.max(r)
  min_i <- which.min(r)
  else_i <- all_i[! all_i %in% c(max_i, min_i)]
  min_threshold <- lower_bound
  max_threshold <- 1-(length(r)-1)*lower_bound
  if (max(r)>max_threshold & min(r)<min_threshold) {
    assign(paste0("r", max_i), eval(parse(text=paste0("min(",max_threshold, ",r",max_i,")"))))
    assign(paste0("r", min_i), eval(parse(text=paste0("max(",min_threshold, ",r",min_i,")"))))
    for(i in else_i){
      assign(paste0("r", i), eval(parse(text=paste0("max(",min_threshold, ",r",i,")"))))
      if (eval(paste0("r", i))>min_threshold){
        assign(paste0("r", i), min_threshold)
      }
    }
  }

  if (max(r)<=max_threshold & min(r)<min_threshold) {
    assign(paste0("r", min_i), eval(parse(text=paste0("max(",min_threshold, ",r",min_i,")"))))
    for(i in else_i){
      assign(paste0("r", i), eval(parse(text=paste0("max(",min_threshold, ",r",i,")"))))
    }
    assign(paste0("r", max_i), eval(parse(text=paste0("1-r",min_i,"-",paste(paste0("r", else_i, collapse = "-"))))))
  }
  assign("tempprob", eval(parse(text=paste("c(",paste(paste0("r", 1:length(success_count)), collapse = ","),")"))))
  Rho <- tempprob
  if(control_arm == "fixed"){
    phi.jk<-Rho*((Rho/(n/(sum(n))))^2)
    phi <- phi.jk
    ratio <- phi.jk/sum(phi.jk)
    previous_rem_sum <- 1-ratio[1]
    phi[1] <- 1/length(success_count)
    after_rem_sum <- 1-phi[1]
    for(zzz in 2:length(success_count)){
      phi[zzz] <- after_rem_sum*ratio[zzz]/previous_rem_sum
    }
  }

  else{
    phi.jk<-Rho*((Rho/(n/(sum(n))))^2)
    phi<-phi.jk/sum(phi.jk)
  }
  all_i <- c(1:length(phi))
  max_i <- which.max(phi)
  min_i <- which.min(phi)
  else_i <- all_i[! all_i %in% c(max_i, min_i)]

  min_threshold <- lower_bound
  max_threshold <- 1-(length(phi)-1)*lower_bound

  for (i in all_i){
    assign(paste0("phi",i),eval(parse(text=paste0("phi[",i,"]"))))
  }

  if (max(phi)>max_threshold & min(phi)<min_threshold) {
    assign(paste0("phi", max_i), eval(parse(text=paste0("min(",max_threshold, ",phi",max_i,")"))))
    assign(paste0("phi", min_i), eval(parse(text=paste0("max(",min_threshold, ",phi",min_i,")"))))
    #assign(paste0("phi", else_i), eval(parse(text=paste0("1-phi",max_i,"-phi",min_i,""))))
    for(i in else_i){
      assign(paste0("phi", i), eval(parse(text=paste0("max(",min_threshold, ",phi",i,")"))))
      if (eval(paste0("phi", i))>min_threshold){
        assign(paste0("phi", i), min_threshold)
      }
    }
  }

  if (max(phi)<=max_threshold & min(phi)<min_threshold) {
    assign(paste0("phi", min_i), eval(parse(text=paste0("max(",min_threshold, ",phi",min_i,")"))))
    for(i in else_i){
      assign(paste0("phi", i), eval(parse(text=paste0("max(",min_threshold, ",phi",i,")"))))
    }
    assign(paste0("phi", max_i), eval(parse(text=paste0("1-phi",min_i,"-",paste(paste0("phi", else_i,collapse = "-"))))))
  }
  assign("phi", eval(parse(text=paste("c(",paste(paste0("phi", 1:length(phi)), collapse = ","),")"))))
  x <- runif(1, 0, 1)
  ind <- findInterval(x, cumsum(phi), all.inside = FALSE, left.open = TRUE) + 1
  return(phi)

}
