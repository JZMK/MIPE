#Jin Zhang & Cong Chen, 2015

IPE_BW <- function(ipe.input, max.iter=200, param.conv=1e-4, distr="weibull",recensor=1){
  ##initialize the output
  result = list()

  ##initialize the survreg computation parameters
  survreg.control(maxiter=1000, rel.tolerance=1e-08)

  ##introduce latent event time
  result$latent.time = ipe.input$event.time

  ##Run Survreg for an initial param. estimate by comparing the groups as randomized
  survreg.out = survreg(Surv(event.time, censor.ind) ~ trt.ind, data =ipe.input, dist=distr)
  result$init <- survreg.out
  censor<-ipe.input$censor.ind
  ##Vector of estimates for algorithm checking
  result$trt.param = matrix(0,max.iter,1)

  for(i in c(1:max.iter)){
      result$trt.param[i,1] = as.numeric(survreg.out$coefficients[2])
	  if(i%%30==0) result$trt.param[i,1]<-result$trt.param[i,1] +.05
	  #result$trt.param[i,1] =ifelse(i==1,0.8,as.numeric(survreg.out$coefficients[2]))
      diff.param =abs(result$trt.param[i,1]-ifelse(i==1,0,result$trt.param[i-1,1]))
      result$conv = ifelse(diff.param < param.conv, "YES", "NO")
      if(result$conv=="YES"){
        result$survreg.out = survreg.out
        break
      }
      else {
	    if (recensor==0){
        tmp.index = which(ipe.input$trt.ind ==0 & ipe.input$switch.ind==1)
        result$latent.time[tmp.index] = ipe.input$switch.time[tmp.index] + exp(-result$trt.param[i,1])*(ipe.input$event.time[tmp.index] - ipe.input$switch.time[tmp.index])
        }

		if (recensor==1){
		phi<-result$trt.param[i,1]		
		result$latent.time = ipe.input$event.time
		tmp.index = which(ipe.input$trt.ind ==0 & ipe.input$switch.ind==1)
        result$latent.time[tmp.index] = ipe.input$switch.time[tmp.index] + exp(-phi)*(ipe.input$event.time[tmp.index] - ipe.input$switch.time[tmp.index])
		if(phi>=0){
		renloc<-which(censor==0 & ipe.input$trt.ind==0 )
		result$latent.time[renloc]<-exp(-phi)*(ipe.input$cutofftime[renloc])
		reeloc<-which(censor==1 & ipe.input$trt.ind==0 )
		comp<-(result$latent.time[reeloc]<=exp(-phi)*ipe.input$cutofftime[reeloc])*1
        result$latent.time[reeloc]<-(1-comp)*(exp(-phi)*ipe.input$cutofftime[reeloc])+comp*(result$latent.time[reeloc])
		}else{
		renloc<-which(censor==0 & ipe.input$trt.ind==0 )
		result$latent.time[renloc]<-(ipe.input$cutofftime[renloc])
		reeloc<-which(censor==1 & ipe.input$trt.ind==0 )
		comp<-(result$latent.time[reeloc]<=ipe.input$cutofftime[reeloc])*1
		result$latent.time[reeloc]<-(1-comp)*ipe.input$cutofftime[reeloc]+comp*(result$latent.time[reeloc])
		}
        ipe.input$censor.ind<-censor
		ipe.input$censor.ind[reeloc]<-comp
		}
		survreg.out = survreg(Surv(result$latent.time, ipe.input$censor.ind) ~ ipe.input$trt.ind, dist=distr)

	  }
  }
  return(result)
}



