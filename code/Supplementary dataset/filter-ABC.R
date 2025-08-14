library('dplyr')
library("spatstat")

## sample filter
### per section
output_sample <- function(test, n){
    sec.no<- paste('C57BL6J-638850.',toString(n), sep='')
    test <-subset(test, brain_section_label == sec.no)
    sec <-subset(test,select=c(parcellation_substructure,x:y, average_correlation_score,class,subclass, supertype, cluster))
    res<-c(unique(sec$subclass))
    ## excitatory
    ie<-res[grepl("Glut", res)]
    sec_ie<- sec[sec$subclass %in% ie,]
    ex.cluster<- length(unique(sec_ie$cluster))
    ## inhibitory
    ii<-res[grepl("Gaba", res)]
    sec_ii<- sec[sec$subclass %in% ii,]
    in.cluster<- length(unique(sec_ii$cluster))
    
    sec_ie <-subset(sec_ie, average_correlation_score >=0.7)
    ex.no <- nrow(sec_ie)
    total<-count(sec_ie, cluster, sort=TRUE)
    etypes<- total[total$n>0.01*ex.no,]$cluster
    fil.ex<- length(etypes)
    
    sec_ii <-subset(sec_ii, average_correlation_score >=0.7)  
    in.no<- nrow(sec_ii)
    total<-count(sec_ii, cluster, sort=TRUE)
    itypes<- total[total$n>max(0.01*in.no,10),]$cluster  ## lower bar
    fil.in<- length(itypes)
    
    ex.data<- c(ex.cluster,ex.no, fil.ex)
    in.data<- c(in.cluster,in.no, fil.in)
    print(ex.data)
    print(in.data)
    newList <- list("ex" = ex.data, "et" = etypes,"edata" = sec_ie, "in" = in.data,"it" = itypes,"idata" = sec_ii)
    return(newList)
}



## create ppp 
### for all cells
create_pp<-function(types, sec_ic){
    pp<-list()
    flag<-FALSE
    for (t in types){
        if (grepl("L2/3", t)| grepl("L4/5", t)|grepl("L5",t)|grepl("L6",t)){
            flag<-TRUE
        }
        cells<- sec_ic[sec_ic$cluster %in% t,]
        ow<-init_window(t, sec_ic, flag)
        p<- ppp(cells$x, cells$y, ow)
        pp[[length(pp)+1]] = p
    }
    return(pp)
}


## initial window 
init_window<- function(t, sec_ic, flag){
    region<- sec_ic[sec_ic$cluster %in% t,]
    ow1<-ripras(region$x, region$y,shape='convex')
    
    if (flag==TRUE){
        p1<-ppp(region$x, region$y, ow1)
        b <- bw.ppl(p1)
        dp1<- density(p1, sigma=b)
        Zcut<- cut(dp1, breaks = c(0,0.001*max(dp1),max(dp1)), include.lowest=TRUE,label=1:2)
        V<-tess(image = Zcut)
        if (length(split(p1, V)) == 1){
            return(ow1)
        }
        p1m<- split(p1, V)[[2]]
        ow2<-Window(p1m)
        ow1<-simplify.owin(as.polygonal(ow2),0.05)
    }
    return(ow1)
}


## refined window for high density region
ref_window<- function(pp, i, pr){
    b<- bw.ppl(pp[[i]])
    dp1<- density(pp[[i]], sigma=b)
    Zcut<- cut(dp1, breaks = c(0,pr*max(dp1),max(dp1)), include.lowest=TRUE,label=1:2)  ## above-average high intensity 
    ## this controls window size
    V<-tess(image = Zcut)
    if (length(split(pp[[i]], V)) == 1){
        return(Window(pp[[i]]))
    }
    p1m<- split(pp[[i]], V)[[2]]
    ow2<-Window(p1m)
    ow1<-simplify.owin(as.polygonal(ow2),0.05)
    
    ## Take the largest window
    ow1.data<-as.data.frame(ow1)
    ids<-c(unique(ow1.data$id))
    if (length(ids)>1){  
        ids.m <- 1
        ids.area<-0
        for (i in ids){
            data<-subset(ow1.data,id==i)
            ## window with holes
            if (data$sign[1]==-1){
                next
            }
            bds<-split(data[,c('x','y')],data$i)
            a<-area(owin(poly=bds))
            if (a > ids.area){
                ids.area<-a
                ids.m <- i
            }
        }
        ow1.data<-subset(ow1.data,id==ids.m)
        bds<-split(ow1.data[,c('x','y')],ow1.data$id)
        ow1<-owin(poly=bds)
    }
    ## too small
    skip_to_next <- FALSE
    tryCatch(erosion.owin(ow1, r=0.5*b), error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { return(square(0.0001))}     
    ow<- erosion.owin(ow1, r=0.5*b)
    return(ow)
}


## calculate 70%/50% conf ratio
calc_conf<-function(pp, types, test, n){
    conf <- c()
    sec.no<- paste('C57BL6J-638850.',toString(n), sep='')
    test <-subset(test, brain_section_label == sec.no)
    for (i in 1:length(types)){
        cells<- test[test$cluster %in% types[[i]],]
        ow<-Window(pp[[i]])
        p<- ppp(cells$x, cells$y, ow)
        ratio <- npoints(pp[[i]])/npoints(p)
        conf<-append(conf,ratio)
    }
    return(conf)
}



## calculate HC R
calc_R<-function(pp, types, pr){
    nn<- c()
    for (i in 1:length(types)){
        m<- quantile(nndist(pp[[i]]), probs = pr)[[1]]
        nn<-append(nn,m)
    }
    return(nn)
}

### Intensity filter
intense_qualify<- function(pp, types, nn, conf){
    qq <-c()
    dist <- c()
    con <- c()
    for (i in 1:length(types)){
        print(types[i])
        ## filter dense region
        pr<-0.8
        while (area(ref_window(pp,i, pr))<(10*nn[i])^2){   ## area filter
            if (pr <= 0.05){
                pr <- 0.05
                break
            }
            else{
                pr<- pr-0.05
            }
        }
     
        ow1<- ref_window(pp,i, pr)
        flag <- TRUE
        if (area(ow1)<0.5*(10*nn[i])^2){
            flag <- FALSE
        }
        print(pr)
        Window(pp[[i]])<-ow1
        inten <- intensity(pp[[i]])
        if (inten>15*(0.1/nn[i])^2 & conf[i]>0.1& flag == TRUE){  
            qq<-append(qq,types[i])
            dist<-append(dist, nn[i])
            con<-append(con,conf[i])
        }
    }
    newList <- list("cluster" = qq, "dist" = dist, 'conf'=con)
    return(newList)
}