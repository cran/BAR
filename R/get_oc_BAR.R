get_oc_BAR <- function(success_prob, n_burn_in, tot_num, block_size,
                       power_c = "n/2N", lower_bound = .05, reptime,
                       control_arm = "", output = "", seed = 100){
  set.seed(seed)
  if (!control_arm %in% c("", "fixed")){
    stop("Please make sure the argument of control_arm is specified as default or fixed")
  }
  if(n_burn_in * length(success_prob) > tot_num){
    stop("Burn-in sample size greater than total sample size")
  }
  if((tot_num - n_burn_in * length(success_prob)) %% block_size != 0){
    stop("Please make sure the sample size after burn-in period can be divided by block size")
  }
  assigning_prob_list <- list()
  alloc_prob_matrix <- matrix(NA, reptime, length(success_prob))
  pat_num_matrix <- matrix(NA, reptime, length(success_prob))
  for(m in 1:reptime){
    shapeBetas<-matrix(rep(0.5,2*length(success_prob)),nrow=length(success_prob),byrow=T)
    for(ii in 1:length(success_prob)){
      assign(paste("xx", ii, sep = ""), NULL)
    }
    for(ii in 1:length(success_prob)){
      assign(paste("xx", ii, sep = ""), c(eval(parse(text = paste("xx", ii, sep = ""))), rbinom(n_burn_in,1,success_prob[ii])))
    }
    N<-rep(n_burn_in,length(success_prob))
    for(ii in 1:length(success_prob)){
      shapeBetas[ii,1] <- shapeBetas[ii,1]+sum(eval(parse(text = paste("xx", ii, sep = ""))))
      shapeBetas[ii,2] <- shapeBetas[ii,2]+N[ii]-sum(eval(parse(text = paste("xx", ii, sep = ""))))
    }
    for(ii in 1:length(success_prob)){
      assign(paste("a", ii, sep = ""), rbeta(1000,shapeBetas[ii,1],shapeBetas[ii,2]))
    }
    #pmaxfunction <- paste0("pmax(",paste(paste0("a", 1:length(group.level)), collapse =","),")")
    pmaxfunction <- paste0("pmax(",paste(paste0("a", 1:length(success_prob)), collapse =","),")")
    for(ii in 1:length(success_prob)){
      #assign(paste("c", ii, sep = ""), eval(parse(text=paste0("length(which(a",ii,"==", paste0("pmax(a1,", "a", ii, ")"),"))/1000"))))
      #assign(paste("c", ii, sep = ""), eval(parse(text=paste0("length(which(a",ii,"==", success_prob))))
      assign(paste("c", ii, sep = ""), eval(parse(text=paste0("length(which(a",ii,"==", pmaxfunction,"))/1000"))))
    }
    if (control_arm == "fixed"){
      previous_rem_sum <- 1-c1
      c1 <- 1/length(success_prob)
      after_rem_sum <- 1-c1
      for(xx in 2:length(success_prob)){
        assign(paste("c", xx, sep = ""), eval(parse(text = paste0("c", xx, "/previous_rem_sum*after_rem_sum"))))
      }
    }
    prob_matrix <- matrix(NA, tot_num-n_burn_in*length(success_prob), length(success_prob))
    for (j in 1:((tot_num-n_burn_in*length(success_prob))/block_size)){
      if(power_c == "n/2N"){
        lambda <- (sum(N)+j-1)/(2*tot_num)
      }
      else{
        lambda <- power_c
      }
      lambdafunction1 <- paste(paste0("c", 1:length(success_prob),"^lambda"), collapse = "+")
      lambdafunction2 <- paste(paste0("c", 2:length(success_prob),"^lambda"), collapse = "+")
      if(control_arm == "fixed"){
        rem_sum <- 1-c1
        r1 <- c1
        for(z in 2:length(success_prob)){
          assign(paste("r", z, sep = ""), eval(parse(text=paste0("(c",z,"^lambda)/(",lambdafunction2,")*rem_sum"))))
        }

      }
      else{
        for(z in 1:length(success_prob)){
          assign(paste("r", z, sep = ""), eval(parse(text=paste0("(c",z,"^lambda)/(",lambdafunction1,")"))))
        }
      }
      assign("r", eval(parse(text=paste("c(",paste(paste0("r", 1:length(success_prob)), collapse = ","),")"))))

      #r1<-c1^lambda/(c1^lambda+c2^lambda+c3^lambda)
      #r2<-c2^lambda/(c1^lambda+c2^lambda+c3^lambda)
      #r3<-c3^lambda/(c1^lambda+c2^lambda+c3^lambda)
      all_i <- c(1:length(r))
      max_i <- which.max(r)
      min_i <- which.min(r)
      else_i <- all_i[! all_i %in% c(max_i, min_i)]

      min_threshold <- lower_bound
      max_threshold <- 1-(length(r)-1)*lower_bound

      #control arm 1/length(r) CANNOT less than min_threshold OR more than max_threshold
      #           if it less than min_threshold: 1/length(r)<0.05  ====> length(r)>20
      #           if it more than max_threshold: 1/length(r)>1-(length(r)-1)*0.05
      #                                   ====>  length(r)**2 -21length(r)+20>0
      #                                   ====>  (length(r)-20)(length(r)-1)>0
      #                                   ====>  length(r)>20 & length(r)>1
      #control arm 1/length(r) CANNOT be the maximum of the vector:
      #           for example: if it is the maximum of the vector
      #                     sum(   arm1         arm2         arm3         arm4         arm5 ...)
      #                     sum(   1/length(r) <1/length(r) <1/length(r) <1/length(r) <1/length(r) ...)
      #                     <1

      if (max(r)>max_threshold & min(r)<min_threshold) {
        assign(paste0("r", max_i), eval(parse(text=paste0("min(",max_threshold, ",r",max_i,")"))))
        assign(paste0("r", min_i), eval(parse(text=paste0("max(",min_threshold, ",r",min_i,")"))))
        #assign(paste0("r", else_i), eval(parse(text=paste0("1-r",max_i,"-r",min_i))))
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
      assign("tempprob", eval(parse(text=paste("c(",paste(paste0("r", 1:length(success_prob)), collapse = ","),")"))))
      #tempprob=threshold(r1,r2,r3)
      Rho=tempprob
      #gamma=2
      #allocation probability function
      if(control_arm == "fixed"){
        phi.jk<-Rho*((Rho/(N/(sum(N))))^2)
        phi <- phi.jk
        ratio <- phi.jk/sum(phi.jk)
        previous_rem_sum <- 1-ratio[1]
        phi[1] <- 1/length(success_prob)
        after_rem_sum <- 1-phi[1]
        for(zzz in 2:length(success_prob)){
          phi[zzz] <- after_rem_sum*ratio[zzz]/previous_rem_sum
        }
      }
      else{
        phi.jk<-Rho*((Rho/(N/(sum(N))))^2)
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
      #phi<-threshold(phi[1],phi[2],phi[3])
      prob_matrix[j, ] <- phi
      x<-rep(runif(1,0,1), block_size)
      ind <- findInterval(x, cumsum(phi), all.inside = FALSE, left.open = TRUE)+1
      for (i in all_i){
        assign(paste0("N",i),eval(parse(text=paste0("N[",i,"]"))))
      }
      for (i in all_i){
        assign(paste0("N",i), eval(parse(text = paste0("N", i, "+", "sum(ind==", i, ")"))))
      }
      assign("N", eval(parse(text=paste("c(",paste(paste0("N", all_i), collapse = ","),")"))))
      for (i in ind){
        assign(paste0("xx", i), eval(parse(text=paste0("c(xx",i, ",rbinom(1, 1, success_prob[",i,"]))"))))
      }
      ind_table <- table(factor(ind, levels = 1:length(success_prob)))
      which_ind_table <- which(ind_table != 0)
      for(ii in which_ind_table){
        for(i in 1:ind_table[ii]){
          if(eval(parse(text = paste0("xx", ii, "[length(xx", ii, ")+", -i+1, "]"))) == 1){
            shapeBetas[ii,1]<-shapeBetas[ii,1]+1
          }
          if(eval(parse(text = paste0("xx", ii, "[length(xx", ii, ")+", -i+1, "]"))) == 0){
            shapeBetas[ii,2]<-shapeBetas[ii,2]+1
          }
        }

      }


      #if (x>=0 & x<phi[1]) {
      #  N[1]<-N[1]+1
      #  xx1<-c(xx1,rbinom(1,1,success_prob[1]))
      #  if(xx1[length(xx1)] == 1){
      #    shapeBetas[1,1]<-shapeBetas[1,1]+1
      #  }
      #  if(xx1[length(xx1)] == 0){
      #    shapeBetas[1,2]<-shapeBetas[1,2]+1
      #  }
      #}
      #else if (x<=phi[1]+phi[2]){
      #  N[2]<-N[2]+1
      #  xx2<-c(xx2,rbinom(1,1,success_prob[2]))
      #  if(xx2[length(xx2)] == 1){
      #    shapeBetas[2,1]<-shapeBetas[2,1]+1
      #  }
      #  if(xx2[length(xx2)] == 0){
      #    shapeBetas[2,2]<-shapeBetas[2,2]+1
      #  }
      #}
      #else{
      #  N[3]<-N[3]+1
      #  xx3<-c(xx3,rbinom(1,1,success_prob[3]))
      #  if(xx3[length(xx3)] == 1){
      #    shapeBetas[3,1]<-shapeBetas[3,1]+1
      #  }
      #  if(xx3[length(xx3)] == 0){
      #    shapeBetas[3,2]<-shapeBetas[3,2]+1
      #  }
      #}

      for(ii in 1:length(success_prob)){
        assign(paste("a", ii, sep = ""), rbeta(1000,shapeBetas[ii,1],shapeBetas[ii,2]))
      }
      #pmaxfunction <- paste0("pmax(",paste(paste0("a", 1:length(group.level)), collapse =","),")")
      pmaxfunction <- paste0("pmax(",paste(paste0("a", 1:length(success_prob)), collapse =","),")")
      for(ii in 1:length(success_prob)){
        #assign(paste("c", ii, sep = ""), eval(parse(text=paste0("length(which(a",ii,"==", paste0("pmax(a1,", "a", ii, ")"),"))/1000"))))
        #assign(paste("c", ii, sep = ""), eval(parse(text=paste0("length(which(a",ii,"==", success_prob))))
        assign(paste("c", ii, sep = ""), eval(parse(text=paste0("length(which(a",ii,"==", pmaxfunction,"))/1000"))))
      }
      if (control_arm == "fixed"){
        previous_rem_sum <- 1-c1
        c1 <- 1/length(success_prob)
        after_rem_sum <- 1-c1
        for(xx in 2:length(success_prob)){
          assign(paste("c", xx, sep = ""), eval(parse(text = paste0("c", xx, "/previous_rem_sum*after_rem_sum"))))
        }
      }
      for(d in 1:length(success_prob)){
        Rho[d] <- eval(parse(text = paste0("c", d)))
      }
    }
    assigning_prob_list[[m]] <- prob_matrix
    alloc_prob_matrix[m, ] <- N/sum(N)
    pat_num_matrix[m, ] <- N
  }
  mean_alloc_prob <- apply(alloc_prob_matrix, 2, mean)
  mean_pat_num <- apply(pat_num_matrix, 2, mean)
  if(output == "raw"){
    return(assigning_prob_list)
  }
  return(list(avg_alloc_prob = mean_alloc_prob, avg_pat_num = mean_pat_num))
}
