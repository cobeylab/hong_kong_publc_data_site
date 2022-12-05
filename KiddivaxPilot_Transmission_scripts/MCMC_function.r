#
# R syntax to produce MCMC function used in:
#
# Klick B, Nishiura H, Ng S, et al.
# Transmissibility of seasonal and pandemic influenza in a cohort
# of households in Hong Kong in 2009
# Epidemiology. 2011 Nov;22(6):793-6.
#
# Last updated by Klick B, Fang VJ and Cowling BJ.
# January 3, 2012
#


### basic functions for the mcmc ###

### likelihood
log.like.alt8<-function(q_h.ch, q_h.ad, q_c.ch, q_c.ad, dframe) {
    w_0000<-1
    for (O1 in 0:5) {
        for (O2 in 0:O1) {
            for (O3 in 0:5) {
                for (O4 in 0:O3) {
                    if (O1==O2 && O3==O4 && O1+O3!=0) {
                        p.v<-0
                        assign(paste(paste(paste(paste("w_",O3,sep=""),O3,sep=""),O1,sep=""),O1,sep=""), 0)
                        for (i in 0:O3) {
                            for (j in 0:O1) {
                                p.v<- get(paste(paste(paste(paste("w_",i,sep=""),O3,sep=""),j,sep=""),O1,sep="")) + p.v
                            }
                        }
                        assign(paste(paste(paste(paste("w_",O3,sep=""),O3,sep=""),O1,sep=""),O1,sep=""), 1-p.v)}
                    else {
                        dum<-get( paste(paste(paste(paste("w_",O4,sep=""),O4,sep=""),O2,sep=""),O2,sep=""))
                        val<-dum*choose(O1,O2)*choose(O3,O4)*q_c.ch^(O3-O4)*q_h.ch^((O2+O4)*(O3-O4))*q_c.ad^(O1-O2)*q_h.ad^((O2+O4)*(O1-O2))
                        assign(paste(paste(paste(paste("w_",O4,sep=""),O3,sep=""),O2,sep=""),O1,sep=""), val)}
                }
            }
        }
    }
    llik.sum<-0
    for (i in 1:nrow(dframe)) {
        dum<-log(get(paste(paste(paste(paste("w_",dframe[i,1],sep=""),dframe[i,2],sep=""),dframe[i,3],sep=""),dframe[i,4],sep="")))
        llik.sum<-dum+llik.sum
    }
    llik.sum
}

### mcmc for pH1/sH3
mcmc.fun.1<-function (num,st.q_h.ch=3/4, st.q_h.ad=3/4, st.q_c.ch=7/8, st.q_c.ad=7/8, var.name, mcmc.data, inf.data) {
    cnum.m<-which(names(mcmc.data)==paste(var.name,"_miss",sep=""))
    cnum.i<-which(names(mcmc.data)==paste(var.name,"_inf",sep=""))

    q_h.ch.v<-c(st.q_h.ch,rep(NA,num-1))
    q_h.ad.v<-c(st.q_h.ad,rep(NA,num-1))
    q_c.ch.v<-c(st.q_c.ch,rep(NA,num-1))
    q_c.ad.v<-c(st.q_c.ad,rep(NA,num-1))

    for (i in 2:num) {
        ## home escape probability--children
        new.q_h.ch<-max(q_h.ch.v[i-1]+runif(1,-.05,.05),0.0001)
        new.q_h.ch<-min(0.9999,new.q_h.ch)

        old.like<-log.like.alt8(q_h.ch=q_h.ch.v[i-1], q_h.ad=q_h.ad.v[i-1], q_c.ch.v[i-1], q_c.ad.v[i-1], dframe=inf.data)
        cur.like<-log.like.alt8(q_h.ch=new.q_h.ch, q_h.ad=q_h.ad.v[i-1], q_c.ch.v[i-1],  q_c.ad.v[i-1], dframe=inf.data )
        old.val<-old.like+log(dbeta((1-q_h.ch.v[i-1]),1.5,6))
        cur.val<-cur.like+log(dbeta((1-new.q_h.ch),1.5,6))

        alpha <- cur.val-old.val
        alpha <- exp(alpha)
        #print(c(i,old.val,cur.val,q_h.ch.v[i-1],new.q_h.ch,alpha))
        u <-runif(1)
            ### Now accept or reject candidate
        if (u<alpha) {
            q_h.ch.v[i] <- new.q_h.ch
            }else {
            q_h.ch.v[i] <-q_h.ch.v[i-1]
            }

        ## home escape probability--adults
        new.q_h.ad<-max(q_h.ad.v[i-1]+runif(1,-.05,.05),0.0001)
        new.q_h.ad<-min(0.9999,new.q_h.ad)

        old.like<-log.like.alt8(q_h.ch=q_h.ch.v[i], q_h.ad=q_h.ad.v[i-1], q_c.ch.v[i-1], q_c.ad.v[i-1], dframe=inf.data)
        cur.like<-log.like.alt8(q_h.ch=q_h.ch.v[i], q_h.ad=new.q_h.ad, q_c.ch.v[i-1],  q_c.ad.v[i-1], dframe=inf.data)
        old.val<-old.like+log(dbeta((1-q_h.ad.v[i-1]),1.5,6))
        cur.val<-cur.like+log(dbeta((1-new.q_h.ad),1.5,6))

        alpha <- cur.val-old.val
        alpha <- exp(alpha)
        u <-runif(1)
            ### Now accept or reject candidate
        if (u<alpha) {
            q_h.ad.v[i] <- new.q_h.ad
            }else {
            q_h.ad.v[i] <-q_h.ad.v[i-1]
            }

        ## community escape probability--children
        new.q_c.ch<-max(q_c.ch.v[i-1]+runif(1,-.05,.05),0.0001)
        new.q_c.ch<-min(0.9999,new.q_c.ch)

        old.like<-log.like.alt8(q_h.ch=q_h.ch.v[i], q_h.ad=q_h.ad.v[i], q_c.ch.v[i-1], q_c.ad.v[i-1], dframe=inf.data)
        cur.like<-log.like.alt8(q_h.ch=q_h.ch.v[i], q_h.ad=q_h.ad.v[i], new.q_c.ch, q_c.ad.v[i-1], dframe=inf.data)
        old.val<-old.like+log(dbeta((1-q_c.ch.v[i-1]),1.2,6))
        cur.val<-cur.like+log(dbeta((1-new.q_c.ch),1.2,6))

        alpha <- cur.val-old.val
        alpha <- exp(alpha)
        u <-runif(1)
            ### Now accept or reject candidate
        if (u<alpha) {
            q_c.ch.v[i] <- new.q_c.ch
            }else {
            q_c.ch.v[i] <-q_c.ch.v[i-1]
            }

        ## community escape probability--adults
        new.q_c.ad<-max(q_c.ad.v[i-1]+runif(1,-.05,.05),0.0001)
        new.q_c.ad<-min(0.9999,new.q_c.ad)

        old.like<-log.like.alt8(q_h.ch=q_h.ch.v[i], q_h.ad=q_h.ad.v[i], q_c.ch.v[i], q_c.ad.v[i-1], dframe=inf.data)
        cur.like<-log.like.alt8(q_h.ch=q_h.ch.v[i], q_h.ad=q_h.ad.v[i], q_c.ch.v[i],  new.q_c.ad, dframe=inf.data)
        old.val<-old.like+log(dbeta((1-q_c.ad.v[i-1]),1.2,6))
        cur.val<-cur.like+log(dbeta((1-new.q_c.ad),1.2,6))

        alpha <- cur.val-old.val
        alpha <- exp(alpha)
        u <-runif(1)
            ### Now accept or reject candidate
        if (u<alpha) {
            q_c.ad.v[i] <- new.q_c.ad
            }else {
            q_c.ad.v[i] <-q_c.ad.v[i-1]
            }

        ch_prob<-1- q_c.ch.v[i]*(q_h.ch.v[i]^(pmax(mcmc.data$n_ch_inf-imp.inf+ mcmc.data$n_ad_inf,0)))
        ad_prob<-1- q_c.ad.v[i]*(q_h.ad.v[i]^(pmax(mcmc.data$n_ch_inf+ mcmc.data$n_ad_inf-imp.inf,0)))
        prob<-ifelse(mcmc.data$child==1,ch_prob,ad_prob)
        imp.inf<-rbinom(length(mcmc.data$d1_pmax),1,prob)

        mcmc.data$imput_inf<-ifelse(mcmc.data[,cnum.m],imp.inf,mcmc.data[,cnum.i])
        mcmc.data$imp_ch_inf<-mcmc.data$imput_inf*mcmc.data$child
        mcmc.data$imp_ad_inf<-mcmc.data$imput_inf*mcmc.data$adult

        t.ch<-as.vector(by(mcmc.data$child, mcmc.data$hhID, sum, na.rm=TRUE))
        inf.ch<-as.vector(by(mcmc.data$imp_ch_inf, mcmc.data$hhID, sum, na.rm=TRUE))
        t.ad<-as.vector(by(mcmc.data$adult, mcmc.data$hhID, sum, na.rm=TRUE))
        inf.ad<-as.vector(by(mcmc.data$imp_ad_inf, mcmc.data$hhID, sum, na.rm=TRUE))

        len_fam<-as.vector(by(mcmc.data$child, mcmc.data$hhID, length))

        mcmc.data$n_ch_inf<-rep(inf.ch, time=len_fam)
        mcmc.data$n_ad_inf<-rep(inf.ad, time=len_fam)

        inf.data<-data.frame(x1=inf.ch, x2=t.ch, x3=inf.ad, x4=t.ad)
        }

        retdata<-data.frame(q_h.ch.v=q_h.ch.v, q_h.ad.v=q_h.ad.v, q_c.ch.v=q_c.ch.v, q_c.ad.v=q_c.ad.v)
        return(retdata)
}

### mcmc for sH1
mcmc.fun.2<-function (num,st.q_h.ch=3/4, st.q_h.ad=3/4, st.q_c.ch=7/8, st.q_c.ad=7/8, var.name, mcmc.data, inf.data) {
    cnum.m<-which(names(mcmc.data)==paste(var.name,"_miss",sep=""))
    cnum.i<-which(names(mcmc.data)==paste(var.name,"_inf",sep=""))

    q_h.ch.v<-c(st.q_h.ch,rep(NA,num-1))
    q_h.ad.v<-c(st.q_h.ad,rep(NA,num-1))
    q_c.ch.v<-c(st.q_c.ch,rep(NA,num-1))
    q_c.ad.v<-c(st.q_c.ad,rep(NA,num-1))

    for (i in 2:num) {
        ## home escape probability--children
        new.q_h.ch<-max(q_h.ch.v[i-1]+runif(1,-.05,.05),0.0001)
        new.q_h.ch<-min(0.9999,new.q_h.ch)

        old.like<-log.like.alt8(q_h.ch=q_h.ch.v[i-1], q_h.ad=q_h.ad.v[i-1], q_c.ch.v[i-1], q_c.ad.v[i-1], dframe=inf.data)
        cur.like<-log.like.alt8(q_h.ch=new.q_h.ch, q_h.ad=q_h.ad.v[i-1], q_c.ch.v[i-1],  q_c.ad.v[i-1], dframe=inf.data )
        old.val<-old.like+log(dbeta((1-q_h.ch.v[i-1]),1.5,6))
        cur.val<-cur.like+log(dbeta((1-new.q_h.ch),1.5,6))

        alpha <- cur.val-old.val
        alpha <- exp(alpha)
        #print(c(i,old.val,cur.val,q_h.ch.v[i-1],new.q_h.ch,alpha))
        u <-runif(1)
            ### Now accept or reject candidate
        if (u<alpha) {
            q_h.ch.v[i] <- new.q_h.ch
            }else {
            q_h.ch.v[i] <-q_h.ch.v[i-1]
            }

        ## home escape probability--adults
        new.q_h.ad<-max(q_h.ad.v[i-1]+runif(1,-.05,.05),0.0001)
        new.q_h.ad<-min(0.9999,new.q_h.ad)

        old.like<-log.like.alt8(q_h.ch=q_h.ch.v[i], q_h.ad=q_h.ad.v[i-1], q_c.ch.v[i-1], q_c.ad.v[i-1], dframe=inf.data)
        cur.like<-log.like.alt8(q_h.ch=q_h.ch.v[i], q_h.ad=new.q_h.ad, q_c.ch.v[i-1],  q_c.ad.v[i-1], dframe=inf.data)
        old.val<-old.like+log(dbeta((1-q_h.ad.v[i-1]),1.5,6))
        cur.val<-cur.like+log(dbeta((1-new.q_h.ad),1.5,6))

        alpha <- cur.val-old.val
        alpha <- exp(alpha)
        u <-runif(1)
            ### Now accept or reject candidate
        if (u<alpha) {
            q_h.ad.v[i] <- new.q_h.ad
            }else {
            q_h.ad.v[i] <-q_h.ad.v[i-1]
            }

        ## community escape probability--children
        new.q_c.ch<-max(q_c.ch.v[i-1]+runif(1,-.05,.05),0.0001)
        new.q_c.ch<-min(0.9999,new.q_c.ch)

        old.like<-log.like.alt8(q_h.ch=q_h.ch.v[i], q_h.ad=q_h.ad.v[i], q_c.ch.v[i-1], q_c.ad.v[i-1], dframe=inf.data)
        cur.like<-log.like.alt8(q_h.ch=q_h.ch.v[i], q_h.ad=q_h.ad.v[i], new.q_c.ch, q_c.ad.v[i-1], dframe=inf.data)
        old.val<-old.like+log(dbeta((1-q_c.ch.v[i-1]),1.5,30))
        cur.val<-cur.like+log(dbeta((1-new.q_c.ch),1.5,30))

        alpha <- cur.val-old.val
        alpha <- exp(alpha)
        u <-runif(1)
            ### Now accept or reject candidate
        if (u<alpha) {
            q_c.ch.v[i] <- new.q_c.ch
            }else {
            q_c.ch.v[i] <-q_c.ch.v[i-1]
            }

        ## community escape probability--adults
        new.q_c.ad<-max(q_c.ad.v[i-1]+runif(1,-.05,.05),0.0001)
        new.q_c.ad<-min(0.9999,new.q_c.ad)

        old.like<-log.like.alt8(q_h.ch=q_h.ch.v[i], q_h.ad=q_h.ad.v[i], q_c.ch.v[i], q_c.ad.v[i-1], dframe=inf.data)
        cur.like<-log.like.alt8(q_h.ch=q_h.ch.v[i], q_h.ad=q_h.ad.v[i], q_c.ch.v[i],  new.q_c.ad, dframe=inf.data)
        old.val<-old.like+log(dbeta((1-q_c.ad.v[i-1]),1.5,30))
        cur.val<-cur.like+log(dbeta((1-new.q_c.ad),1.5,30))

        alpha <- cur.val-old.val
        alpha <- exp(alpha)
        u <-runif(1)
            ### Now accept or reject candidate
        if (u<alpha) {
            q_c.ad.v[i] <- new.q_c.ad
            }else {
            q_c.ad.v[i] <-q_c.ad.v[i-1]
            }

        ch_prob<-1- q_c.ch.v[i]*(q_h.ch.v[i]^(pmax(mcmc.data$n_ch_inf-imp.inf+ mcmc.data$n_ad_inf,0)))
        ad_prob<-1- q_c.ad.v[i]*(q_h.ad.v[i]^(pmax(mcmc.data$n_ch_inf+ mcmc.data$n_ad_inf-imp.inf,0)))
        prob<-ifelse(mcmc.data$child==1,ch_prob,ad_prob)
        imp.inf<-rbinom(length(mcmc.data$d1_pmax),1,prob)

        mcmc.data$imput_inf<-ifelse(mcmc.data[,cnum.m],imp.inf,mcmc.data[,cnum.i])
        mcmc.data$imp_ch_inf<-mcmc.data$imput_inf*mcmc.data$child
        mcmc.data$imp_ad_inf<-mcmc.data$imput_inf*mcmc.data$adult

        t.ch<-as.vector(by(mcmc.data$child, mcmc.data$hhID, sum, na.rm=TRUE))
        inf.ch<-as.vector(by(mcmc.data$imp_ch_inf, mcmc.data$hhID, sum, na.rm=TRUE))
        t.ad<-as.vector(by(mcmc.data$adult, mcmc.data$hhID, sum, na.rm=TRUE))
        inf.ad<-as.vector(by(mcmc.data$imp_ad_inf, mcmc.data$hhID, sum, na.rm=TRUE))

        len_fam<-as.vector(by(mcmc.data$child, mcmc.data$hhID, length))

        mcmc.data$n_ch_inf<-rep(inf.ch, time=len_fam)
        mcmc.data$n_ad_inf<-rep(inf.ad, time=len_fam)

        inf.data<-data.frame(x1=inf.ch, x2=t.ch, x3=inf.ad, x4=t.ad)
        }

        retdata<-data.frame(q_h.ch.v=q_h.ch.v, q_h.ad.v=q_h.ad.v, q_c.ch.v=q_c.ch.v, q_c.ad.v=q_c.ad.v)
        return(retdata)
}

# End of script.

