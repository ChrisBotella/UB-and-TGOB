library('ggplot2')
library('glm2')
library('raster')

#####
# Functions
#####

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Simulate an IPP from intensity raster and number of points
simu_IPP = function(r_int,N){
  e = extent(r_int)
  Max = max(getValues(r_int),na.rm=T)
  df = matrix(0,0,1)
  while(dim(df)[1]<N){
    tmp = data.frame(x=runif(2*N,e[1],e[2]),y=runif(2*N,e[3],e[4]),z=runif(N,0,Max))
    v = extract(r_int,tmp[,c('x','y')])
    tmp = tmp[!is.na(v) & tmp$z<v,]
    if(dim(df)[1]==0){df = tmp[,c('x','y'),drop=F]}else{df=rbind(df,tmp[,c('x','y')])}
  }
  if(dim(df)[1]>N){df=df[1:N,]}
  colnames(df)=c('Longitude','Latitude')
  return(df)
}


# wraper uniform points
add.integration_pts.uniform=function(pts,D=NULL,n_0,Bias_killer=100,dom=NULL){
  
  if(is.null(dom)){dom=integration_pts.uniform(D,n_0)}
  dom$pseudo_y=0
  # columns of pts which are absent of dom are added to it, & filled with NAs
  other_col = colnames(pts)[!colnames(pts)%in%colnames(dom)]
  dom[,other_col]= NA
  
  n=NULL
  taxas = unique(as.character(pts$taxa))
  if(is.null(D)){AD = 1}else{AD = raster::area(D)}
  for(e in 1:length(taxas)){
    sp = pts$taxa==taxas[e]
    n = sum(sp)
    pts$w[sp] = AD/(Bias_killer*n)
    pts$pseudo_y[sp] = (Bias_killer*n)/AD
  }
  # columns of dom which are absent of pts are removed
  dom=dom[,colnames(dom)%in%colnames(pts)]
  
  # Merge species and background points
  # with a repetition of background points per species 
  for(e in 1:length(taxas)){dom$taxa=taxas[e];pts=rbind(pts,dom)}
  return(pts)
}

# A green data function
# Remove several useless heavy objects from the glm model object
# in order to save it afterwards
light_glm = function(GlmObject){
  GlmObject$data = 1
  GlmObject$model = 1
  GlmObject$qr$qr = GlmObject$qr$qr[1:2,,drop=F]
  GlmObject$family = "poisson(link = 'log')"
  vect = as.character(GlmObject$formula)
  GlmObject$formula = paste(vect[2],vect[1],vect[3])
  attr(GlmObject$terms,'.Environment') = NULL
  return(GlmObject)
}


tgb_yw = function(pts,Bias_killer=100,MainTaxa,D=NULL){
  # All lines whose taxa identifier is equal to MainTaxa
  # are considerations as observations points of the main
  # species
  n = sum(pts$taxa==MainTaxa)
  if(is.null(D)){A_D = 1}else{A_D = raster::area(D) / (1e4)}
  MainTaxaBool = (!is.na(pts$taxa) & pts$taxa==MainTaxa)
  # All lines that are not of the main species 
  # are considered as background points  
  n_0 = dim(pts)[1]-n
  pts$pseudo_y = MainTaxaBool * n * 100 / A_D
  pts$w = MainTaxaBool * A_D / (Bias_killer*n) +  !MainTaxaBool * (Bias_killer-1)*A_D/(Bias_killer*n_0) 
  return(pts)
}

# Personnal graphical colors 
mycolor=function(Names){
  lof_colors = c('magenta2','deeppink2','mediumorchid2','purple2','darkorchid4')
  k=0
  vec=NULL
  for(i in 1:length(Names)){
    if(Names[i]=="true species density"){
      vec[i] = "red2"
    }else if(Names[i]=="observation density"){
      vec[i] = "goldenrod4"
    }else if(Names[i]=="TG species density"){
      vec[i] = "darkorchid1"
    }else if(Names[i]=="observed points density"){
      vec[i] = "palegreen3"
    }else if(Names[i]%in%c("cla",'ppp')){
      vec[i] = 'blue3'
    }else if(Names[i]=="est_dens"){
      vec[i] = 'grey43'
    }else if(Names[i]=="true_dens"){
      vec[i] = 'black'
    }else if(Names[i]=="est_obs"){
      vec[i] = 'tan4'
    }else if(Names[i]=="p_x"){
      vec[i] = 'black'
    }else if(Names[i]=="true/TG density"){
      vec[i] = 'gray60'
    }else if(Names[i]=="tgb"){
      vec[i] = 'deepskyblue1'
    }else if(regexpr('lof_q',Names[i])>0){
      k=k+1
      #vec[i] = lof_colors[k]
      vec[i] ='black' 
    }
  }
  return(vec)
}



#####
# parameters
##### 

x_0 = 0
y_0 = 0
N_0 = 20000
n_is = c(20000)
nrep = 1
tested_methods = c('cla','tgb')
tested_simtypes = c('cst','lin','cut','hol','gs')

plotSize = 2.

dir = 'C:/Users/user/pCloud local/boulot/data/UB and TGOB/'
setwd(dir)

formu =  'I(axe3) + I(axe3^2)'

df_esp = read.csv('df_esp.csv',sep=";",header=T)
groups = readRDS('groups')
gr_names = names(groups)


#####
# dom & vals
#####

dom = data.frame(Longitude=runif(N_0,x_0,10),Latitude=runif(N_0,y_0,10))
dom$axe3 = dom$Longitude - 5
dom$w = 99 / (N_0*100)

vals = expand.grid(Longitude=seq(0,10,0.03),Latitude=seq(0,10,0.03),esp=NA)
vals$axe3 = vals$Longitude - 5
# CST : constant
vals$cst = 1
# LIN : triangle superior left
vals$lin = as.numeric(vals$Latitude>vals$Longitude)
# CUT: 1 on left, 0 on right
vals$cut = as.numeric(vals$Longitude<4)
# HOL 
vals$hol = exp((vals$axe3+1)^2/2)/(10+exp((vals$axe3+1)^2/2))
# GS : gaussian with mean zero 
vals$gs = exp(-vals$axe3^2)  / (sqrt(2*pi)*1)

#####
# RUN
#####

# variables matrix
X = matrix(NA,dim(vals)[1],2)
X[!is.na(vals$axe3),] = model.matrix( ~ I(axe3) +I(axe3^2) , vals )[,2:3]

loop = expand.grid(group=gr_names,simtype=tested_simtypes,n_i=n_is,rep=1:nrep)
res = expand.grid(meth=tested_methods,esp='main',group=gr_names,simtype=tested_simtypes,n_i=n_is,rep=1:nrep,time=NA,cor2=NA,cor2_obs=NA)

for(r in 1:dim(loop)[1]){
  print(loop[r,,drop=F])
  n_i = loop$n_i[r]
  rep = loop$rep[r]
  type = loop$simtype[r]
  gr = loop$group[r]
  especes = groups[gr][[1]] 
  
  for(i in 1:length(especes)){
    print(paste('esp',especes[i]))
    
    par = as.numeric(df_esp[df_esp$sp_id == especes[i],c('I.axe3','I.axe3.2.')])
    vals$esp = exp( X %*% par )
    vals$int = vals[,as.character(type)] * vals$esp
    
    r_int = rasterFromXYZ(vals[,c('Longitude','Latitude','int')])
    # df ~ IPP(r_int)
    df=simu_IPP(r_int,n_i)
    df$taxa=as.character(especes[i])
    if(i==1){DF = df}else{DF=rbind(DF,df)}
  }
  DF$axe3 = DF$Longitude - 5
  
  # CLA, on n'estime que la 4eme espèce
  spid = especes[length(especes)]
  print('UB')
  mod=list()
  deb=Sys.time()
  print(paste('inference on species',spid))
  df = add.integration_pts.uniform(DF[DF$taxa==as.character(spid),],D=NULL,n_0=20,dom=dom)

  
  #mod[[1]] = light_glm(cla_glm(df,formu)) 
  mod[[1]] = light_glm(
    glm2(data=df,
        formula=as.formula(paste('pseudo_y ~ ',formu)),
        weights = df$w,
        family=poisson(link="log"),
        control=list(maxit=100,trace=F)))
  
  
  print(Sys.time()-deb)
  setwd(dir)
  saveRDS(mod,paste('cla_rep',rep,'_n',n_i,'_simtype',type,'_grname',gr,sep=""))
  
  # tgob
  print('TGOB')
  mod=list()
  deb=Sys.time()
  print(paste('inference on species',spid))
  tab=tgb_yw(DF,Bias_killer=100,MainTaxa=as.character(spid))
  mod[[1]] = light_glm(
              glm(data=tab,
                  formula=as.formula(paste('pseudo_y ~', formu)),
                  weights = tab$w,
                  family=poisson(link="log"),
                  control=list(maxit=100,trace=F)))
  print(Sys.time()-deb)
  setwd(dir)
  saveRDS(mod,paste('tgb_rep',rep,'_n',n_i,'_simtype',type,'_grname',gr,sep=""))
}

#####
# pull results
#####

res = expand.grid(meth = tested_methods,simtype = tested_simtypes, gr = gr_names)
res$axe3 = NA
res$axe3.2 = NA

setwd(dir)

for(r in 1:dim(res)[1]){
  file = paste(res$meth[r],'_rep1_n',n_is,'_simtype',res$simtype[r],'_grname',res$gr[r],sep="")
  print(file)
  mod = readRDS(file)
  res$axe3[r] = mod[[1]]$coefficients['I(axe3)']
  res$axe3.2[r] = mod[[1]]$coefficients['I(axe3^2)']
}

#####
# MultiPlot for UB
#####

pas = 0.05
X=seq(-5,5,pas)

ploted_simtypes = tested_simtypes
ploted_methods = c('cla')
reord_grnames = gr_names[regexpr('BG_flat',gr_names)>0]

#dev.off()
plot_list = list()
k=1
for(gr in reord_grnames){
  for(type in ploted_simtypes){
    
    tab = expand.grid(x=X,curve=c('true species density','observation density','observed points density',ploted_methods),y=NA)
    # true species distribution
    especes = groups[gr][[1]]
    spid = especes[length(especes)]
    cd = tab$curve=='true species density'
    tmp = exp(df_esp$I.axe3[df_esp$sp_id==spid]*X + df_esp$I.axe3.2.[df_esp$sp_id==spid]*X^2)
    tab$y[cd] = tmp / sum(tmp*pas)
    # observation density
    cd = tab$curve=="observation density"
    if(type%in%c("cst",'std')){
      tab$y[cd]=1/length(X)
    }else if(type=="cut"){
      tab$y[cd]= as.numeric(X<(-1))
    }else if(type=="lin"){
      tab$y[cd]= (5-X)/sum(5-X)
    }else if(type=="ultracor"){
      tab$y[cd]= -min(X/25) + (X<(-2)|X>0) *(X/25) + (X>=-2 & X<=0)* -1/25
    }else if(type=="gs"){
      tab$y[cd]= exp(-X^2)/sqrt(2*pi)
      print(type)
    }else if(type=="cre"){
      tab$y[cd]= 1.*(X<(-1) | X>1)
    }else if(type=="hol"){
      tab$y[cd] = exp((X+1)^2/2)/(10+exp((X+1)^2/2))
    }
    tab$y[cd] = tab$y[cd]/sum(tab$y[cd]*pas)
    # observed points density
    cd = tab$curve == "observed points density"
    tab$y[cd] = tab$y[tab$curve=="observation density"] * tab$y[tab$curve=="true species density"]
    tab$y[cd] = tab$y[cd] / sum(tab$y[cd]*pas)
    
    # estimates
    for(m in ploted_methods){
      cd = res$meth==m & res$simtype==type & res$gr==gr
      tab$y[tab$curve==m] = exp(res$axe3[cd] * X + res$axe3.2[cd] * X^2)
      tab$y[tab$curve==m] = tab$y[tab$curve==m] / sum(tab$y[tab$curve==m]*pas)
    }
    
    p = ggplot(tab,aes(x=x,y=y,group=curve))
    p = p+geom_line(aes(colour=curve,size=curve,linetype=curve))
    p = p+scale_colour_manual(values=mycolor(levels(tab$curve)))
    p = p+scale_size_manual(values = plotSize+c(8,3,5,1))
    p = p+scale_linetype_manual(values= c('solid','dotdash','dotdash','dashed'))
    p = p+scale_x_continuous(breaks=seq(-5,5,2))
    p = p+ggtitle(element_blank())
    p = p+scale_y_continuous(limits=c(0,.7))
    p = p + ylab("Density (normalised  on [-5,5])")
    p=p+theme_bw()+theme(panel.grid.major = element_blank()#element_line(color="grey")
                         ,panel.grid.minor = element_blank(),
                         axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.x = element_text(size=40),
                         axis.text.y = element_text(size=40),
                         panel.border = element_blank(),
                         legend.position = "none",
                         axis.line.x = element_line(color="black", size = 2),
                         axis.line.y = element_line(color="black", size = 2))
    #print(p)
    plot_list[[k]]=p
    #plot_list[[k]]=ggplot(tab,aes(x=x,y=y,group=curve))+geom_line(aes(colour=curve))+scale_x_continuous(breaks=seq(-5,5,2))+ggtitle(paste('obs:',type,'; target:',gr))+theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ theme(legend.position = "none")
    k=k+1
  }
}


setwd(dir)
png(file="multiplot_UB_notext.png",width=5000,height=6000)
multiplot(plotlist = plot_list , cols=4)
dev.off()


#####
# MultiPlots for TG0B
#####


# a equal "flat", "thick", or "thin"
a = "thick"

pas = 0.05
X=seq(-5,5,pas)

ploted_simtypes = tested_simtypes
ploted_methods = c('tgb')
reord_grnames = gr_names[regexpr(paste('BG_',a,sep=""),gr_names)>0]

#dev.off()
plot_list = list()
k=1
for(gr in reord_grnames){
  for(type in ploted_simtypes){
    
    tab = expand.grid(x=X,curve=c('true species density',
                                  'observation density',
                                  'observed points density',
                                  'TG species density',
                                  "true/TG density",
                                  ploted_methods),y=NA)
    # true species distribution
    especes = groups[gr][[1]]
    spid = especes[length(especes)]
    cd = tab$curve=='true species density'
    tmp = exp(df_esp$I.axe3[df_esp$sp_id==spid]*X + df_esp$I.axe3.2.[df_esp$sp_id==spid]*X^2)
    tab$y[cd] = tmp / sum(tmp*pas)
    # observation density
    cd = tab$curve=="observation density"
    if(type%in%c("cst",'std')){
      tab$y[cd]=1/length(X)
    }else if(type=="cut"){
      tab$y[cd]= as.numeric(X<(-1))
    }else if(type=="lin"){
      tab$y[cd]= (5-X)/sum(5-X)
    }else if(type=="ultracor"){
      tab$y[cd]= -min(X/25) + (X<(-2)|X>0) *(X/25) + (X>=-2 & X<=0)* -1/25
    }else if(type=="gs"){
      tab$y[cd]= exp(-X^2)/sqrt(2*pi)
      print(type)
    }else if(type=="cre"){
      tab$y[cd]= 1.*(X<(-1) | X>1)
    }else if(type=="hol"){
      tab$y[cd] = exp((X+1)^2/2)/(10+exp((X+1)^2/2))
    }
    tab$y[cd] = tab$y[cd]/sum(tab$y[cd]*pas)
    # observed points density
    #cd = tab$curve == "observed points density"
    #tab$y[cd] = tab$y[tab$curve=="observation density"] * tab$y[tab$curve=="true species density"]
    #tab$y[cd] = tab$y[cd] / sum(tab$y[cd]*pas)
    
    # TG points density
    bg = especes[1]
    tmp = exp(df_esp$I.axe3[df_esp$sp_id==bg]*X+df_esp$I.axe3.2.[df_esp$sp_id==bg]*X^2)
    tab$y[tab$curve=="TG species density"] = tmp / sum(tmp*pas)
    
    # true/tg density 
    cd = tab$curve == "true/TG density"
    tab$y[cd] = tab$y[tab$curve=="true species density"] / tab$y[tab$curve=="TG species density"]
    tab$y[cd] = tab$y[cd] / sum(tab$y[cd]*pas)
    
    # estimates
    for(m in ploted_methods){
      cd = res$meth==m & res$simtype==type & res$gr==gr
      tab$y[tab$curve==m] = exp(res$axe3[cd] * X + res$axe3.2[cd] * X^2)
      tab$y[tab$curve==m] = tab$y[tab$curve==m] / sum(tab$y[tab$curve==m]*pas)
    }
    
    # plot
    p = ggplot(tab,aes(x=x,y=y,group=curve))
    p = p+geom_line(aes(colour=curve,size=curve,linetype=curve))
    p = p+scale_size_manual(values = plotSize+c(8,3,1,1,2,5))
    p = p+scale_colour_manual(values=mycolor(levels(tab$curve)))
    p = p+scale_linetype_manual(values= c('solid','dotdash','dotdash','solid','dashed','longdash'))
    p = p+scale_x_continuous(breaks=seq(-5,5,2))
    #p = p+ggtitle(paste('obs_',type,'; ',gr,sep=""))
    p = p+ggtitle(element_blank())
    p = p+scale_y_continuous(limits=c(0,.7))
    #p = p+theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position = "none",title=element_text(size=10))
    p = p + ylab("Density (normalised  on [-5,5])")
    p=p+theme_bw()+theme(panel.grid.major = element_blank()#element_line(color="grey")
                         ,panel.grid.minor = element_blank(),
                         axis.title.x=element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.x = element_text(size=60),
                         axis.text.y = element_text(size=60),
                         panel.border = element_blank(),
                         legend.position = "none",
                         axis.line.x = element_line(color="black", size = 2),
                         axis.line.y = element_line(color="black", size = 2))
    #print(p)
    plot_list[[k]]=p
    #plot_list[[k]]=ggplot(tab,aes(x=x,y=y,group=curve))+geom_line(aes(colour=curve))+scale_x_continuous(breaks=seq(-5,5,2))+ggtitle(paste('obs:',type,'; target:',gr))+theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ theme(legend.position = "none")
    k=k+1
  }
}


setwd(dir)
png(file=paste("multiplot_TGOB",a,"_notext.png",sep=""),width=5000,height=6000)
multiplot(plotlist = plot_list , cols=4)
dev.off()
