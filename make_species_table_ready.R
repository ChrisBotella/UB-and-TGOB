# created on 30/04/2018
# table of virtual species niches parameters


dir = 'C:/Users/user/pCloud local/boulot/data/UB and TGOB/'

dd = expand.grid( bg = c('1,','2,','3,') , esp = c('4','5','6','7'))
insidelist= paste(paste('c(',dd$bg,dd$esp,')'),collapse=',')
eval(parse(text=paste( 'groups = list(',insidelist,')')))

dd = expand.grid( bg = c('"BG_thin','"BG_thick','"BG_flat') , esp = c('SP_close_thin"','SP_close_thick"','SP_far_thin"','SP_far_thick"'))
insidelist= paste(paste(dd$bg,dd$esp),collapse=',')
eval(parse(text=paste( 'names(groups) = list(',insidelist,')')))

mus = c(0.,0.,0.,
        -1,
        -1.,
        -4,
        -4)
sigs = c(1.,2.,20.,
         .6,
         1.5,
         .6,
         1.5)
df_esp = data.frame(sp_id=1:length(mus),
                    I.axe3 = 0,
                    I.axe3.2. = 0)
for(i in 1:dim(df_esp)[1]){
  df_esp$I.axe3[i] = mus[i]/(sigs[i]^2)
  df_esp$I.axe3.2.[i] = -1/(2*sigs[i]^2)
}

setwd(dir)
write.table(df_esp,'df_esp.csv',sep=";",col.names=T,row.names=F)
saveRDS(groups,'groups')
