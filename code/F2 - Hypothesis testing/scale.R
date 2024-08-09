library("spatstat")
library('dplyr')
library('poolr')


## input dataset: sag-data
### contains all clusters in Isocortex region of 5 sections


### Area and intensity table by parcellations
calc_area<- function(data, n, cl){
    sec.no<- paste('Zhuang-ABCA-3.0',toString(n), sep='')
    sec<-subset(data, brain_section_label ==sec.no)
    c1<- sec[sec$cluster %in% cl,]
    total<-count(c1, parcellation_substructure)
    total <- rename(total,c('Counts'='n'))
    regions = c(total$parcellation_substructure)
    ar<-c()
    for (r in regions){
        region<- sec[sec$parcellation_substructure==r,]
        if (length(region$x)<=2){
            ar<-append(ar,1)
            next
        }
        ow<-ripras(region$x, region$y,shape='convex')
        p1 <- ppp(region$x, region$y, ow)
        dp1<- density(p1, sigma=bw.ppl(p1))
        Zcut<- cut(dp1, breaks = c(0,0.1*max(dp1),max(dp1)), include.lowest=TRUE,label=1:2)
        V<-tess(image = Zcut)
        if (length(split(p1, V)) == 1){
            ar<-append(ar,area(ow))
            next
        }
        p1m<- split(p1, V)[[2]]
        ow2<-Window(p1m)
        ow1<-simplify.owin(as.polygonal(ow2),0.05)
        ar<-append(ar,area(ow1))
    }
    total$area<-ar
    total$intensity <- total$'Counts'/ar
    df<- total[order(total$intensity,decreasing = TRUE), ]
    return(df)
}

### Find parcellation substructures
### This function is global across sections
find_parc <- function(data, nn, cl, flag){
    if (flag==TRUE){
        sub <- unique(data[data$cluster %in% cl,]$subclass)
        data<- data[data$subclass %in% sub,]  ## all cells in subclass
    }
    dd<-list()
    res<-list()
    c<-0
    a<-0
    for (n in nn){
        df<-calc_area(data, n, cl)
        dd[[length(dd)+1]] = df
        c=c+sum(df$Counts)
        a=a+sum(df$area)
    }
    ab <- c/a
   
    for (i in 1:length(nn)){
        df<-dd[[i]]
        med=ceiling(length(df$parcellation_substructure)/2)
        r1 = c(df$parcellation_substructure)[1:med]
        df <- df[df$intensity>ab ,]
        rr<- c(df$parcellation_substructure)
        rf<-union(rr,r1)    
        res[[as.character(nn[i])]] = rf
    }
    return(res)
}

## Take all parcellation substructures
all_parc<-function(data, nn, cl){
    res<-list()
    for (n in nn){  
        sec.no<- paste('Zhuang-ABCA-3.0',toString(n), sep='')
        sec<-subset(data, brain_section_label ==sec.no)
        c1<- sec[sec$cluster %in% cl,]
        total<-count(c1, parcellation_substructure)
        regions = c(total$parcellation_substructure)
        res[[as.character(n)]] = regions
    }
    return(res)
}



### Initial window defined by parcellation structures
### per section per cluster
init_window<-function(r,data, n, cl, flag){
    sec.no<- paste('Zhuang-ABCA-3.0',toString(n), sep='')
    sec<-subset(data, brain_section_label ==sec.no)
    region<- sec[sec$parcellation_substructure %in% r,]
    if (flag== TRUE){
        sub <- unique(region[region$cluster %in% cl,]$subclass)
        region<- region[region$subclass %in% sub,]  ## all cells in subclass
    }
    ow<-ripras(region$x, region$y,shape='convex')
    p1 <- ppp(region$x, region$y, ow)
    b <- bw.ppl(p1)
    dp1<- density(p1, sigma=b)
    Zcut<- cut(dp1, breaks = c(0,0.1*max(dp1),max(dp1)), include.lowest=TRUE,label=1:2)
    V<-tess(image = Zcut)
    c1<- region[region$cluster %in% cl,]
    owr<-ripras(c1$x, c1$y,shape='convex')
    if (length(split(p1, V)) == 1){
        ow2<-intersect.owin(ow,owr)
        newList <- list("data" = c1, "ow" = ow2)
        return(newList)
    }
    p1m<- split(p1, V)[[2]]
    ow2<-intersect.owin(Window(p1m),owr)
    ow1<-simplify.owin(as.polygonal(ow2),0.05)
    
    ow1.data<-as.data.frame(ow1)
    ids<-c(unique(ow1.data$id))  
    ## filter on erosions (small windows will be removed)
    ##if (length(ids)>1){ 
    ##    q = c()
    ##    for (i in ids){
    ##        data<-subset(ow1.data,id==i)
    ##        ## window with holes
    ##        if (data$sign[1]==-1){
    ##            next
    ##        }
    ##        bds<-split(data[,c('x','y')],data$id)
    ##        w<-owin(poly=bds)
    ##        skip_to_next <- FALSE
    ##        tryCatch(erosion.owin(w, r=0.5*b), error = function(e) { skip_to_next <<- TRUE})
    ##        if(skip_to_next) { 
    ##           next
    ##        }
    ##        q<-append(q,i)  
    ##    }
   
    ##    ow1.data<-subset(ow1.data,id %in% q)
    ##    bds<-split(ow1.data[,c('x','y')],ow1.data$id)
    ##    ow1<-owin(poly=bds)
    ##}
    ## ow1 <- erosion.owin(ow1, r=0.5*b)   
    newList <- list("data" = c1, "ow" = ow1)
    return(newList)
}


### Kernel smoother
### This function applies to all types and sections based on initial windows
final_window<-function(init, f){
    data<-init$data
    p1 <- ppp(data$x, data$y, init$ow, marks = data[,c('cluster_confidence_score')])
    b<-bw.ppl(p1)
    dp1<- density(p1, sigma=0.5*b)
    if (f == TRUE){  ## no overlap
        Zcut<- cut(dp1, breaks = c(0,mean(dp1),max(dp1)), include.lowest=TRUE,label=1:2)
    }
    else{  ## overlaps
        Zcut<- cut(dp1, breaks = c(0,mean(dp1),mean(dp1)+2*sd(dp1)), include.lowest=TRUE,label=1:2)
    }
    V<-tess(image = Zcut)
    p1m<- split(p1, V)[[2]]
    ow2<-Window(p1m)
    ow1<-simplify.owin(as.polygonal(ow2),0.05)
    
    p2 <- ppp(data$x, data$y, erosion.owin(ow1, r=0.25*b))
    if (is.empty(Window(p2)) | npoints(p2)< 10){  ## no erosion
        Window(p1)<- ow1 
    }
    else{
        ow1 <- erosion.owin(ow1, r=0.25*b)
        Window(p1)<- ow1
    }
        
    return(p1)
}

    
    
    
## Create ppp for types per section
### The uniform pipeline function
create_pp <- function(data, rr=NULL, nn, type){
    flag <- FALSE
    if (grepl("L2/3", type)| grepl("L4/5", type)|grepl("L5",type)|grepl("L6",type)){
        flag<-TRUE
    }
    pp <- list() 
    x <- unique(data[data$subclass %in% type,]$cluster)
    y <- unique(data[data$supertype %in% type,]$cluster)
    z <- unique(data[data$cluster %in% type,]$cluster)
    cl<- Reduce(union, list(x, y, z))
    if (is.null(rr)){
        rr <-find_parc(data,nn,cl, flag)
    }
    else{
        rr <- all_parc(data,nn, cl)
    }
    for (n in nn){           
        r <- rr[[as.character(n)]]   ## non-zero parcellations  
        init <- init_window(r, data, n, cl, flag)
        pf <- final_window(init, TRUE)
        ## comment for inhibitory types without projection error (return larger window)
        if (ecdf(nndist(pf))(0.01) > 0){
            pf <- final_window(init, FALSE)
        }
        if (!is.null(pf)){
            pp[[as.character(n)]] = pf
        }
    }
    return(pp)
}


## HT
### Use HC PP as benchmark if there is no overlap
### Use Poisson CSR as benchmark if there is overlap
get_pp_2d<-function(pp,rm, flag){
    rr<- seq(0.01,rm,by=0.001)
    ps<-matrix(0,length(rr), length(pp))  
    for (i in 1:length(pp)){
        npp <-npoints(pp[[i]])
        Wpp <-Window(pp[[i]])
        if (flag == TRUE){
            ## CSR   
            pp.base <-runifpoint(npp, Wpp, nsim=199)
        }else{
            ### HC PP
            pp.base <-rSSI(npp, r=0.01, Wpp, nsim=199)
        }
        for (r in 1:length(rr)){
            ps[r,i]<- dclf.test(pp[[i]], Lest, simulate=pp.base, nsim=199, rinterval=c(0,rr[r]), alternative='less',verbose=FALSE)$p.value
        }  
    }
    return(ps)
}

### Predomination table by parcellations
### req layer-wise input data
calc_dom <- function(data, n){
    sec.no<- paste('Zhuang-ABCA-3.0',toString(n), sep='')
    sec<-subset(data, brain_section_label ==sec.no)
    total<-count(sec, parcellation_substructure)
    total <- rename(total,c('total_count'='n'))
    df <- sec %>% group_by(parcellation_substructure,cluster) %>% 
      summarise(total_count=n(),.groups = 'drop') %>% as.data.frame()
    df<-df %>% group_by(parcellation_substructure) %>% top_n(1, total_count)
    df <- rename(df,c('cluster_dom'='cluster','cluster_count'='total_count'))
    total_table<-merge(total, df, by='parcellation_substructure')
    return(total_table)
}