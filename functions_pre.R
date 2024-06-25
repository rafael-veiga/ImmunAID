library(tidyverse)
library(car)
library(fpp)
library(VIM)
library(boot)

library(mgcv)

interpolation_fun = function(x,y,x_ref){
  dref = data.frame(x_ref = x_ref)
  dref$y = NA
  dref$y[1] = 1
  dref$y[length(dref$y)] = 0
  dxy = data.frame(x=x,y=y)
  
  for(ix in x_ref){
    aux = dxy[dxy$x==ix,]
    fil = dref$x==ix
    fil[is.na(fil)]=FALSE
    dref$y[fil] = c(aux$y[1],aux$y[length(aux$y)])
  }
  #calculate new points
  x_ini = 0
  y_ini = 1
  x_fim = 1
  y_fim = 0
  
  for(i in 1:(length(dref$x_ref)-1)){
    if(!is.na(dref$y[i])){
      x_ini = dref$x[i]
      y_ini = dref$y[i]
    }else{
      j=i+1
      while(is.na(dref$y[j])){
        j = j+1
      }
      x_fim = dref$x_ref[j]
      y_fim = dref$y[j]
      if(y_fim==y_ini){
        dref$y[i] = y_ini 
      }else{
        a = (y_fim-y_ini)/(x_fim-x_ini)
        b=y_fim-a*x_fim
        dref$y[i] = a*dref$x[i]+b
      }
      
    }
  }
  return(dref$y)
}
# auc_mean_fun = function(sen,spe,numb){
#   #x= spe
#   #y=sen
#   
#   x_aux = sort(unique(spe))
#   x_ref = rep(NA,length(x_aux)*2)
#   con = 1
#   for(i in seq(1,length(x_ref),by=2)){
#     x_ref[i] = x_aux[con]
#     x_ref[i+1] = x_aux[con]
#     con = con+1
#   }
#    
#   cv = unique(numb)
#   y_list = list()
#   
#   for(i in cv){
#     fil = numb==i
#     y_list[[length(y_list)+1]] = interpolation_fun(x=spe[fil],y=sen[fil],x_ref = x_ref)
#   }
#   y_ref = rep(NA,length(x_ref))
#   y_ref_l = rep(NA,length(x_ref))
#   y_ref_h = rep(NA,length(x_ref))
#   for(i in 1:length(x_ref)){
#     y_aux = rep(NA,length(y_list))
#     for(j in 1:length(y_list)){
#       y_aux[j] = y_list[[j]][i]
#     }
#     y_ref[i] = mean(y_aux)
#     inter = 1.96*(sd(y_aux)/sqrt(length(y_aux)))
#     y_ref_l[i] = y_ref[i] - inter
#     y_ref_h[i] = y_ref[i] + inter
#   }
#  return(list(specificities = x_ref,sensitivities_mean = y_ref,sensitivities_l = y_ref_l,sensitivities_hn = y_ref_h)) 
# }




recur_reg <- function(df_,df_t,lis){
  if(length(lis)<3){
    return("erro")
  }
  tryCatch({
    df_ = df_ %>% select(all_of(c(lis,"out")))
    df_t = df_t %>% select(all_of(c(lis,"out")))
    model =  glm(out~., data= df_ ,family = "binomial")
    if(model$converged){
      return(predict(model,newdata = df_t))
    }else{
      model =  glm(out~., data= df_ ,family = "binomial")
      lis = names(sort(summary(model)$coefficients[,4],decreasing = TRUE))
      return(recur_reg(df_,df_t,lis[2:length(lis)])) 
    }
  }, warning = function(w) {
    model =  glm(out~., data= df_ ,family = "binomial")
    lis = names(sort(summary(model)$coefficients[,4],decreasing = TRUE))
    lis = gsub("`", "", lis)
    lis = lis[lis!="(Intercept)"]
    return(recur_reg(df_,df_t,lis[2:length(lis)]))
  }, error = function(w) {
    
    lis = lis[1:(length(lis)-1)]
    return(recur_reg(df_,df_t,lis[2:length(lis)]))
  })
}

open_file_p <- function(file){
  df = read_csv(paste0("./data/AutoGate_results_patients/",file))
  col = colnames(df)[2:length(colnames(df))]
  df <- df %>% mutate_at(col, as.double)
  df$batch =  unlist(strsplit(file,"_statistics.csv"))[1]
  return(df)
}


open_file_s <- function(file){
  df = read_csv(paste0("./data/AutoGate_results_IC/",file))
  col = colnames(df)[2:length(colnames(df))]
  df <- df %>% mutate_at(col, as.double)
  df$batch =  unlist(strsplit(file,"_IC_statistics.csv"))[1]
  return(df)
}

#fil_col
get_marks <- function(col){
  marks = list()
  pop = list()
  fil_other = (!grepl(" in ",col))
  fil_pop = (grepl(" in ",col) & (!grepl(" median ",col)))
  fil_mark = grepl(" median in ",col)
  other = col[fil_other]
  aux = col[fil_mark]
  aux = str_split(aux," median in ")
  for(i in aux){
    marks = append(marks,i[1])
    pop = append(pop,i[2])
  }
  aux = col[fil_pop]
  aux = str_split(aux," in ")
  for(i in aux){
    pop = append(pop,i[1])
    pop = append(pop,i[2])
  }
  
  type = c("other","pop_in_pop","mark_in_pop")
  res = list()
  res[["mark"]]= unique(unlist(marks))
  res[["pop"]] = unique(unlist(pop))
  res[["other"]] = other
  res[["type"]] = type
  
  return(res)
}

percent_missing <- function(column) {
  missing_values <- sum(is.na(column))
  percent_missing <- (missing_values / length(column)) * 100
  return(round(percent_missing, 2))
}

get_fil <- function(col,type="all",pop="all",mark="all"){
  fil_t = rep(FALSE,length(col))
  for(t in type){
    if(t=="all"){
      fil_t=rep(TRUE,length(col))
      break
    }
    if(t=="other"){
      fil_t = (fil_t) | (col=="Sample") | (col=="batch") | (col=="id") | (col=="type") | (col=="age") | (col=="sex") | (col=="disease") | (col=="fil_dis")
    }
    if(t=="pop_in_pop"){
      fil_t = (fil_t) | (grepl(" in ",col) & (!grepl(" median ",col)) )
    }
    if(t=="mark_in_pop"){
      fil_t = (fil_t) | grepl(" median in ",col)
    }
  }
  fil_p = rep(FALSE,length(col))
  for(p in pop){
    if(p=="all"){
      fil_p=rep(TRUE,length(col))
      break
    }
    f = as.logical(sapply(col,function(s){
      if(grepl(" in ",s)){
        return (unlist(str_split(s," in "))[2])
      }
      return(s)
    })==p)
    ini = as.logical(sapply(col,function(s){
      if(grepl(" in ",s) & (!grepl(" median ",s))){
        return (unlist(str_split(s," in "))[1])
      }
      return(s)
    })==p)
    f_ini = f | ini
    fil_p = fil_p | f_ini
  }
  fil_m = rep(FALSE,length(col))
  for(m in mark){
    if(m=="all"){
      fil_m=rep(TRUE,length(col))
      break
    }
    f = as.logical(sapply(col,function(s){return (unlist(str_split(s," median in "))[1])})==m)
    fil_m = fil_m | f
  }
  return(fil_m & fil_p & fil_t)
}

mfi_cell_pop.fit <- function(df,pop,tax=0.8){
marks = names(pop)
col = colnames(df)
batch = unique(df$batch)
lmark = list()
for(m in 1:length(marks)){
  lbatch = list()
  for(i in 1:length(batch)){
    lbatch[[i]] = list(l=NA,m=NA,h=NA)
  }
  names(lbatch) = batch
  lbatch$ref_l = NA
  lbatch$ref_m = NA
  lbatch$ref_h = NA
  lmark[[m]] = lbatch
}
names(lmark) = marks
for(m in 1:length(marks)){
  if(!is.na(pop[[m]])){
    a = df[,pop[[m]]][[1]]
    med = rep(NA,length(batch))
    l = rep(NA,length(batch))
    h = rep(NA,length(batch))
    for(b in batch){
      v=a[df$batch==b]
      med[b] = mean(v,na.rm=TRUE)
      lmark[[m]][[b]]$m = med[b]
      l[b] = quantile(v,probs=(1-tax),na.rm=TRUE)[[1]]
      lmark[[m]][[b]]$l = l[b]
      h[b] = quantile(v,probs=tax,na.rm=TRUE)[[1]]
      lmark[[m]][[b]]$h = h[b]
    }
    lmark[[m]]$ref_m = mean(med,na.rm=TRUE)
    lmark[[m]]$ref_l =mean(l,na.rm=TRUE)
    lmark[[m]]$ref_h =mean(h,na.rm=TRUE)
  }
}
return(lmark)
}



prop.fit <- function(df,var,tax=0.8){
  batch = unique(df_p$batch)
  med = list()
  low = list()
  hig = list()
  ref_l = rep(NA,length(var))
  ref_m = rep(NA,length(var))
  ref_h = rep(NA,length(var))
  for(v in 1:length(var)){
    a = df[[var[v]]]
    m = rep(NA,length(batch))
    l = rep(NA,length(batch))
    h = rep(NA,length(batch))
    for(b in 1:length(batch)){
      m[b] = mean(a[df$batch==batch[b]],na.rm = TRUE)
      l[b] = quantile(a[df$batch==batch[b]],probs=(1-tax),na.rm=TRUE)
      h[b] = quantile(a[df$batch==batch[b]],probs=tax,na.rm=TRUE)
    }
    med[[v]] = m
    low[[v]] = l
    hig[[v]] = h
    ref_l[v] = mean(l,na.rm=TRUE)
    ref_m[v] = mean(m,na.rm=TRUE)
    ref_h[v] = mean(h,na.rm=TRUE)
  }
  return(list(batch=batch,med = med,low=low,hig=hig,ref_h=ref_h,ref_m=ref_m,ref_l=ref_l,var=var))
}

prop.transform <- function(df,model){
  df_ = df
  lo1 = logit_fun(0)
  lo2 = logit_fun(100)
  for(v in 1:length(model$var)){
    #print(v)
    ref_l = logit_fun(model$ref_l[v])
    ref_m = logit_fun(model$ref_m[v])
    ref_h = logit_fun(model$ref_h[v])
    for(b in 1:length(model$batch)){
      #print(paste0("b ",b))
      low = logit_fun(model$low[[v]][b])
      med = logit_fun(model$med[[v]][b])
      hig = logit_fun(model$hig[[v]][b])
      a1 = (lo1 - ref_l)/(lo1 - low)
      a2 = (ref_l - ref_m)/(low - med)
      a3 = (ref_m - ref_h)/(med - hig)
      a4 = (ref_h - lo2)/(hig - lo2)
      b1 = ref_l - a1*low
      b2 = ref_l - a2*low
      b3 = ref_h - a3*hig
      b4 = ref_h - a4*hig
      x = df_[df_$batch==model$batch[b],colnames(df_)==model$var[v]]
      xl = logit_fun(as.vector(x)[[1]])
      fil1 = xl<=low
      fil1[is.na(fil1)] = FALSE
      fil2 = (xl>low) & (xl<=med)
      fil2[is.na(fil2)] = FALSE
      fil3 = (xl>med) & (xl<=hig)
      fil3[is.na(fil3)] = FALSE
      fil4 = xl>hig
      fil4[is.na(fil4)] = FALSE
      df_[df_$batch==model$batch[b],colnames(df_)==model$var[v]][fil1,1]=alogit_fun(a1*xl[fil1] + b1) 
      df_[df_$batch==model$batch[b],colnames(df_)==model$var[v]][fil2,1]=alogit_fun(a2*xl[fil2] + b2)
      df_[df_$batch==model$batch[b],colnames(df_)==model$var[v]][fil3,1]=alogit_fun(a3*xl[fil3] + b3)
      df_[df_$batch==model$batch[b],colnames(df_)==model$var[v]][fil4,1]=alogit_fun(a4*xl[fil4] + b4)
    }
  }
  return(df_)
}

mfi_cell_pop.transform <- function(df, model){
df_ = df
marks = names(model)
col = colnames(df_)
batch = unique(df_$batch)
for(m in marks){
  var = col[get_fil(col,mark=m)]
  ref_l = model[[m]]$ref_l
  ref_m = model[[m]]$ref_m
  ref_h = model[[m]]$ref_h
  if((!is.na(ref_l)) & (ref_l!=ref_m) & (ref_m != ref_h) ){
    for(b in batch){
      low = model[[m]][[b]]$l[[1]]
      med = model[[m]][[b]]$m[[1]]
      hig = model[[m]][[b]]$h[[1]]
      if((!is.na(low)) & (low!=med) & (med != hig) ){
        a1 = (ref_l)/(low)
        a2 = (ref_l - ref_m)/(low - med)
        a3 = (ref_m - ref_h)/(med - hig)
        b2 = ref_l - a2*low
        b3 = ref_h - a3*hig
        for(v in var){
          filcol = colnames(df_)==v
          filbath = df_$batch==b
          fil1 = filbath  & (df_[,filcol]<low)
          fil1[is.na(fil1)] = FALSE
          fil2 = filbath & (df_[,filcol]>=low) & (df_[,filcol]<med)
          fil2[is.na(fil2)] = FALSE
          fil3 = filbath & (df_[,filcol]>=med)
          fil3[is.na(fil3)] = FALSE
          df_[fil1,filcol] = a1*df_[fil1,filcol]
          df_[fil2,filcol] = a2*df_[fil2,filcol] + b2
          df_[fil3,filcol] = a3*df_[fil3,filcol] + b3
        }
      }else{
        print(paste0("mark ",m," batch ",b," not transform"))
      }
    }
  }else{
    print(paste0("mark ",m," not transform"))
  }
}
return(df_)
}


get_cell_pop <- function(df){
res = get_marks(colnames(df))
marks = res$mark
batch = unique(df$batch)
out = list()
for(m in marks){
  aux = df[,get_fil(colnames(df),mark = m)]
  col = colnames(aux)
  v = rep(NA,length(col)) 
  for (i in 1:length(col)){
    c = aux[,i]
    b = df$batch
    b = b[!is.na(c)]
    c = c[!is.na(c)]
    
    av = rep(NA,length(batch))
    for(bi in 1:length(batch)){
      if(batch[bi] %in% b){
        me =  mean(c[b==batch[bi]])
        if(m=="CD86" | m=="CD11c" | m=="CD57" | m=="CD40"){
          if(m>600){
            if(length(c[b==batch[bi]])>2){
              av[bi] = sd(c[b==batch[bi]])/me
            }
          }
        }else{
          if(m=="HLA-DR" | m=="CD27" | m=="CD38" | m=="CD94"){
            if(m>2000){
              if(length(c[b==batch[bi]])>2){
                av[bi] = sd(c[b==batch[bi]])/me
              }
            }
          }else{
            if(m>300){
              if(length(c[b==batch[bi]])>2){
                av[bi] = sd(c[b==batch[bi]])/me
              }
            }
          }
        }
      }
    }
    v[i] = mean(av)
  }
  if(length(v[!is.na(v)])>0){
    i = which.min(v)
    out[m] = col[i]
  }else{
    out[m] = NA
  }
}
return(out)
}
indirect_df_PROP <- function(df,type="no transform"){
  mix <- function(df,v){
    tryCatch(
      {
        model = Anova(aov(formula = as.formula(paste0(paste0("v",v)," ~ batch+sex+age")), data = df), type="III")
        #model = lmer(paste0(v,"~ 1+age+sex+(1|batch)"),data=df)
        #suppressWarnings(readLines(url))
      },
      warning = function(cond) {
        tryCatch({
          model = Anova(aov(formula = as.formula(paste0(paste0("v",v)," ~ batch+sex")), data = df), type="III")
        }, warning = function(cond) {
          tryCatch({
            model = Anova(aov(formula = as.formula(paste0(paste0("v",v)," ~ batch+age")), data = df), type="III")
          },warning = function(cond) {
            tryCatch({
              model = Anova(aov(formula = as.formula(paste0(paste0("v",v)," ~ batch")), data = df), type="III")
            },warning = function(cond) {model = 2
            })
          })
        })
      })
    return(model)
  }
  col = colnames(df)
  df = df[,get_fil(col,type=c("other","pop_in_pop"))]
  col = colnames(df)
  col_aux = get_fil(col,type="pop_in_pop")
  colnames(df)[col_aux] = sapply(1:length(colnames(df)[col_aux]), function(i){return(paste0("v",i))})
  var = col[col_aux]
  col_aux = colnames(df)
  pvalue = rep(NA,length(var))
    for(v in 1:length(var)){
      print(v)
      # v = "v141"
      aux = df[,c("batch","sex","age",paste0("v",v))]
      aux = aux[!is.na(aux[,4]),]
      aux2 = table(aux$batch)
      l = names(aux2[aux2<5])
      aux2 = aux
      for(b in l){
        #  b = l[1]
        if(!is.na(b)){
          aux2 = aux2[aux2$batch!=b,]
        }
      }
      
      if(nrow(aux2)>=10){
        if(length(unique(aux2$batch))>1){
        aux2 <- aux2 %>% mutate_at(c("age", paste0("v",v)), ~(scale(.) %>% as.vector))
        model = mix(aux2,v)
        if(!is.numeric(model)){
          pvalue[v] = model$`Pr(>F)`[row.names(model)=="batch"]
        }
        }
      }
    }
  aux = data.frame(vars=var,pvalue=pvalue)
  aux$type = type
  return(aux)
}

indirect_df_MFI <- function(df,type="no transform"){
  mix <- function(df,v){
    tryCatch(
      {
        model = Anova(aov(formula = as.formula(paste0(paste0("v",v)," ~ batch+sex+age")), data = df), type="III")
      },
      warning = function(cond) {
        tryCatch({
          model = Anova(aov(formula = as.formula(paste0(paste0("v",v)," ~ batch+sex")), data = df), type="III")
        }, warning = function(cond) {
          tryCatch({
            model = Anova(aov(formula = as.formula(paste0(paste0("v",v)," ~ batch+age")), data = df), type="III")
          },warning = function(cond) {
            tryCatch({
              model = Anova(aov(formula = as.formula(paste0(paste0("v",v)," ~ batch")), data = df), type="III")
            },warning = function(cond) {model = 2
            })
          })
        })
      })
    return(model)
  }
  col = colnames(df)
  df = df[,get_fil(col,type=c("other","mark_in_pop"))]
  col = colnames(df)
  col_aux = get_fil(col,type="mark_in_pop")
  colnames(df)[col_aux] = sapply(1:length(colnames(df)[col_aux]), function(i){return(paste0("v",i))})
  var = col[col_aux]
  col_aux = colnames(df)
  pvalue = rep(NA,length(var))
  for(v in 1:length(var)){
    aux = df[,c("batch","sex","age",paste0("v",v))]
    aux = aux[!is.na(aux[,4]),]
    aux2 = table(aux$batch)
    l = names(aux2[aux2<5])
    aux2 = aux
    for(b in l){
      if(!is.na(b)){
        aux2 = aux2[aux2$batch!=b,]
      }
    }
    if(nrow(aux2)>=10){
      if(length(unique(aux2$batch))>1){
        aux2 <- aux2 %>% mutate_at(c("age", paste0("v",v)), ~(scale(.) %>% as.vector))
        model = mix(aux2,v)
        if(!is.numeric(model)){
          pvalue[v] = model$`Pr(>F)`[row.names(model)=="batch"]
        }
      }
    }
  }
  aux = data.frame(vars=var,pvalue=pvalue)
  aux$type = type
  return(aux)
}



correct_by_POP <- function(df_fit,df_trans,prop=FALSE){
  type = "transform by pop"
  df = df_fit[,get_fil(colnames(df_fit),type = c("other","mark_in_pop"))]
  df_t = df_trans
  res = get_marks(colnames(df_fit))
  mark = res$mark
  batch = unique(df_fit$batch)
  for(m in mark){
    # m = mark[1]
    print(m)
    aux = df[,get_fil(colnames(df),mark = m)]
    col = colnames(aux)
    v = rep(NA,length(col))
    for (i in 1:length(col)){
      # i = 1
      c = aux[,i]
      b = df$batch
      b = b[!is.na(c)]
      c = c[!is.na(c)]
      
      av = rep(NA,length(batch))
      for(bi in 1:length(batch)){
        if(batch[bi] %in% b){
          me =  mean(c[b==batch[bi]])
          if(m>300){
            if(length(c[b==batch[bi]])>2){
            av[bi] = sd(c[b==batch[bi]])/me
            }
          }
        }
      }
      v[i] = mean(av)
    }
    if(length(v[!is.na(v)])>0){
      i = which.min(v)
      c = aux[,i]
      b = df$batch
      b = b[!is.na(c)]
      c = c[!is.na(c)]
      av = rep(NA,length(batch))
      for(bi in 1:length(batch)){
        av[bi] =  mean(c[b==batch[bi]])
      }
      me = mean(av,na.rm=TRUE)
      for(b in 1:length(batch)){
        tax = me/av[b]
        df_t[df_t$batch==batch[b],get_fil(colnames(df_t),mark = m)] = df_t[df_t$batch==batch[b],get_fil(colnames(df_t),mark = m)]*tax
      }
      print(col[i])
    }else{
      print("no found pop")
    }
    
  }
  return(df_t)
  }


remove_batch <- function(df,size=3){
  df_t = df
  t = as.data.frame(table(df_t$batch,df_t$disease))
  batch = unique(df_t$batch)
  t = t[t$Var2=="Negative_control",]
  for(b in batch){
    if(t[t$Var1==b,]$Freq<size){
      df_t = df_t[df_t$batch!=b,]
    }
  }
  return(df_t)
}


IC.fit <- function(df){
  col = colnames(df)
  res = get_marks(col)
  fil_mfi = get_fil(col,type = "mark_in_pop")
  batch = unique(df$batch)
  aux = matrix(NA,nrow = length(batch),ncol = length(res$mark))
  colnames(aux) = res$mark
  row.names(aux) = batch
  mfi_cor = aux
  for(m in 1:length(res$mark)){
    # m = res$mark[1]
    fil = get_fil(col,mark =res$mark[m])
    aux = df[,fil & fil_mfi]
    aux = apply(aux,1, function(l){return(median(l,na.rm=TRUE))})
    aux = data.frame(value=aux,batch=df$batch)
    batch = unique(aux$batch)
    ref = median(aux$value,na.rm=TRUE)
    for(b in 1:length(batch)){
      #  b = batch[1]
      aux2 = aux[aux$batch==batch[b],]
      tax = ref/median(aux2$value,na.rm=TRUE)
      mfi_cor[b,m] = tax
    }
  }
  return(mfi_cor)
}

IC.transform <- function(df,tax){
  col = colnames(df)
  marks = colnames(tax)
  batch = row.names(tax)
  for(m in 1:length(marks)){
    # m = 1
    fil = get_fil(col,mark =marks[m])
    for(b in 1:length(batch)){
      #  b = batch[1]
      df[df$batch==batch[b],fil] = df[df$batch==batch[b],fil]*tax[b,m]
    }
  }
  return(df)
}

direct_df <-function(df,lab="no transform",pro=FALSE){
  res = get_marks(colnames(df))
  batch = unique(df$batch)
  mark = res$mark
  data = data.frame()
  if(!pro){
    for(i in 1:length(mark)){
      for(b in 1:length(batch)){
        aux2 = df[df$batch==batch[b],get_fil(colnames(df), mark = mark[i])]
        data = bind_rows(data,data_frame(batch=batch[b],mark=mark[i],value=mean(colMeans(aux2,na.rm = TRUE),na.rm = TRUE)))
      }
    }
    data$mark = as.factor(data$mark)
  }else{
    for(i in 1:length(res$pop)){
      print(paste0(i," de ",length(res$pop)))
      for(b in 1:length(batch)){
        aux2 = df[df$batch==batch[b],get_fil(colnames(df), pop = res$pop[i],type = "pop_in_pop")]
        data = bind_rows(data,data_frame(batch=batch[b],pop = res$pop[i],value=mean(colMeans(aux2,na.rm = TRUE),na.rm = TRUE)))
      }
  }
  }
  
  return(data)
}

remove_missing <- function(df,per){
  df1 = df
  #df1$type = NULL
  vari = colnames(df1)[get_fil(colnames(df1),type=c("mark_in_pop","pop_in_pop"))]
  for(v in vari){
    tab = as.data.frame(table(is.na(df1[,colnames(df1)==v])))
    if("TRUE" %in% tab$Var1){
      perc = tab$Freq[tab$Var1=="TRUE"]/ sum(tab$Freq)
      
    }else{
      perc = 0
    }
    
    if(perc>per){
      #print(colnames(df1))
      df1 = df1[,colnames(df1)!=v]
    }
  }
  return(df1)
}



box_cox_trans <- function(df,fil){
  df1 = df[,fil]
  df_min = rep(NA,length(colnames(df1)))
  lampda = rep(NA,length(colnames(df1)))
  std_mean =  rep(NA,length(colnames(df1)))
  std_sd =  rep(NA,length(colnames(df1)))
  for(i in 1:length(colnames(df1))){
    df_min[i] = min(df1[,i],na.rm = TRUE)-1
    df1[,i] = df1[,i] - df_min[i]
    lampda[i] = BoxCox.lambda(df1[,i][!is.na(df1[,i])], method = c("loglik"), lower = -5, upper = 5)
    if(lampda[i]==0){
      df1[,i] = log(df1[,i])
    }else{
      df1[,i] = (sign(df1[,i])*((abs(df1[,i]))^lampda[i])-1)/lampda[i]
    }
    #standat
    std_mean[i] = colMeans(df1[,i],na.rm = TRUE)[[1]]
    std_sd[i] = sd(unlist(df1[!is.na(df1[,i]),i]))
    df1[,i] = (df1[,i]-std_mean[i])/std_sd[i]
  }
  df[,fil] = df1
  trans = data.frame(lampda=lampda,std_mean=std_mean,std_sd = std_sd,df_min = df_min)
  return(list(df=df,trans=trans))
}

box_cox_trans_add <- function(df,trans,fil){
  df1 = df[,fil]
  lampda = trans$lampda
  std_mean =  trans$std_mean
  std_sd =  trans$std_sd
  df_min = trans$df_min
  for(i in 1:length(colnames(df1))){
    #add min
    df1[,i] = df1[,i] - df_min[i]
    #inv boxcox
    if (lampda[i] == 0) {
      df1[,i] = log(df1[,i])
    }else{
      df1[,i] =  (sign(df1[,i])*((abs(df1[,i]))^lampda[i])-1)/lampda[i]
    }
    #inv standart
    df1[,i] = (df1[,i]-std_mean[i])/std_sd[i]
    
    
  }
  df[,fil] = df1
  return(df)
}

inv_box_cox_trans <- function(df,trans,fil){
  df1 = df[,fil]
  lambda = trans$lampda
  std_mean =  trans$std_mean
  std_sd =  trans$std_sd
  df_min = trans$df_min
  for(i in 1:length(colnames(df1))){
    #inv standart
    df1[,i] = df1[,i] * std_sd[i]+std_mean[i]
    #inv boxcox
    if (lambda[i] == 0) {
      df1[,i] = exp(df1[,i])
    }else{
      df1[,i] =  (df1[,i] * lambda[i] + 1)^(1/lambda[i])
    }
    #add min
    df1[,i] = df1[,i] + df_min[i]
  }
  df[,fil] = df1
  return(df)
}

imputation <- function(df,fil,k=3){
  df1 = df[,fil]
  df1 = kNN(df1,k = 3,imp_var = FALSE)
  df[,fil] = df1
  return(df)
}
logit_fun <- function(p){
  a = p/100
  a[a<0.00001] = 0.00001
  a[a>0.99999] = 0.99999
  a = boot::logit(a)
  return(a)
}

alogit_fun <- function(a){
  a_ = boot::inv.logit(a)
  a_[a_<=0.00001] = 0
  a_[a_>=0.99999] = 1
  a_ = a_*100
  return(a_)
}




gam_all_fit <- function(df,var,type){
  disease = unique(df$disease)
  batch = unique(df$batch)
  # ref by disease
  ref = list()
  for(i in 1:length(disease)){
    fil_dis = df$disease==disease[i]
    if(type=="mfi"){
      lo1 = 0
      #lo1 = median(apply(df[fil_dis,] %>% select(all_of(var)),2,function(con){return(quantile(con,probs = 0.05,na.rm = TRUE))}),na.rm=TRUE)
      l = median(apply(df[fil_dis,] %>% select(all_of(var)),2,function(con){return(quantile(con,probs = 0.25,na.rm = TRUE))}),na.rm=TRUE)
      m = median(apply(df[fil_dis,] %>% select(all_of(var)),2,function(con){return(median(con,na.rm = TRUE))}),na.rm=TRUE)
      h = median(apply(df[fil_dis,] %>% select(all_of(var)),2,function(con){return(quantile(con,probs = 0.75,na.rm = TRUE))}),na.rm=TRUE)
      lo2 = median(apply(df[fil_dis,] %>% select(all_of(var)),2,function(con){return(quantile(con,probs = 0.95,na.rm = TRUE))}),na.rm=TRUE)
    }else{
      lo1 = logit_fun(0)
      aux = df[fil_dis,] %>% select(all_of(var))
      aux = logit_fun(aux[[1]])
      #lo1 = quantile(aux,probs = 0.25,na.rm = TRUE)[[1]]
      l = quantile(aux,probs = 0.25,na.rm = TRUE)[[1]]
      m = quantile(aux,probs = 0.50,na.rm = TRUE)[[1]]
      h = quantile(aux,probs = 0.75,na.rm = TRUE)[[1]]
      lo2 = logit_fun(1)
    }
    if(is.na(l) | is.na(m) | is.na(h)| lo1==l | l==m | m==h | h==lo2){
      ref[[i]] = list(lo1=NA,l=NA,m=NA,h=NA,lo2=NA)
    }else{
      ref[[i]] = list(lo1=lo1,l=l,m=m,h=h,lo2=lo2)
    }
    ref[[i]] = list(lo1=lo1,l=l,m=m,h=h,lo2=lo2)
  }
  names(ref) = disease
  # weight by batch
  aux = df %>% select(all_of(var))
  model = list()
  for(iw in 1:length(batch)){
    fil_batch = df$batch==batch[iw]
    x = list()
    y = list()
    weit= list()
    for(i in 1:length(disease)){
      #print(paste0("iw ",iw, " i",i))
      va = aux[fil_batch & (df$disease == disease[i]),] %>% select(all_of(var))
      va = apply(va,1,function(con){return(median(con,na.rm = TRUE))})
      if(length(va[!is.na(va)])==0){
        #print(paste0("iw ",iw, " i",i))
        next
      }
      #print(i)
      w = length(va[!is.na(va)])
      if(type=="mfi"){
        log1 = 0
        #lo1 = quantile(va,probs = 0.05,na.rm = TRUE)[[1]]
        l = quantile(va,probs = 0.25,na.rm = TRUE)[[1]]
        m = quantile(va,probs = 0.50,na.rm = TRUE)[[1]]
        h = quantile(va,probs = 0.75,na.rm = TRUE)[[1]]
        lo2 = quantile(va,probs = 0.95,na.rm = TRUE)[[1]]
        
      }else{
        lo1 = ref[[i]]$lo1
        lo2 = ref[[i]]$lo2
        va = aux[fil_batch & (df$disease == disease[i]),] %>% select(all_of(var))
        va = apply(va,2,logit_fun)
        #lo1 = quantile(va,probs = 0.05,na.rm = TRUE)[[1]]
        l = quantile(va,probs = 0.25,na.rm = TRUE)[[1]]
        m = quantile(va,probs = 0.50,na.rm = TRUE)[[1]]
        h = quantile(va,probs = 0.75,na.rm = TRUE)[[1]]
        #lo2 = quantile(va,probs = 0.95,na.rm = TRUE)[[1]]
      }
      if((!is.na(l)) & (!is.na(m)) & (!is.na(h)) & lo1!=l & l!=m & m!=h & h!=lo2){
        x = append(x,ref[[i]]$lo1)
        y = append(y,ref[[i]]$lo1)
        weit = append(weit,w)
        x = append(x,l)
        y = append(y,ref[[i]]$l)
        weit = append(weit,w)
        x = append(x,m)
        y = append(y,ref[[i]]$m)
        weit = append(weit,w)
        x = append(x,h)
        y = append(y,ref[[i]]$h)
        weit = append(weit,w)
        x = append(x,lo2)
        y = append(y,ref[[i]]$lo2)
        weit = append(weit,w)
      }
    }
    if(length(x)>0){
      x = unlist(x)
      y = unlist(y)
      weit = unlist(weit)
      #print(length(x))
      data = data.frame(x=x,y=y,weit=weit)
      
      k = 4
      mn = gam(y~s(x,k=k,bs="cr"),data = data,weights = weit)
      ms = smoothCon(s(x,k=k,bs="cr"),data,knots=NULL)[[1]]
      F = mono.con(ms$xp)
      G <- list(X=ms$X,C=matrix(0,0,0),sp=mn$sp,p=ms$xp,y=y,w= weit)
      G$Ain <- F$A
      G$bin <- F$b
      G$S <- ms$S         
      G$off <- 0          
      p <- pcls(G)
     model[[iw]] <-list(p=p,ms=ms) 
     
    }else{
      model[[iw]] = NA
    }
   
  }
  names(model) = batch
  return(list(model=model,type=type))
}



gam_all_transform = function(df,var,model){
  df_ = df
  batch = names(model$model)
  for(v in var){
    fil_var = colnames(df_)==v
    for(b in 1:length(batch)){
      if(!is.na(model$model[[b]])[1]){
        fil = df_$batch==batch[b] & (!is.na(df[,fil_var]))
        x = df_[fil,fil_var]
        if(nrow(x)>0){
        if(model$type=="mfi"){
          df_[fil,fil_var] = Predict.matrix(model$model[[b]]$ms, data.frame(x = x[[1]])) %*% model$model[[b]]$p
        }else{
          df_[fil,fil_var] = alogit_fun(Predict.matrix(model$model[[b]]$ms, data.frame(x = logit_fun(x[[1]]))) %*% model$model[[b]]$p)
        }
        }
      }
    }
  }
  return(df_)
}

auc_mean_fun_dis = function(dat){
  interpolation_fun = function(x,y,x_ref,y_list){
    ix = 2
    dref = data.frame(x_ref = x_ref)
    dref$y = NA
    dref$y[1] = 1
    dref$y[length(dref$y)] = 0
    #dref$y[length(dref$y)] = 0
    dxy = data.frame(x=x,y=y)
    
    x_uni = unique(dxy$x)
    for(ix in x_uni){
      aux_x = dxy[dxy$x==ix,]
      fil = dref$x==ix
      dref$y[fil] = c(aux_x$y[1],aux_x$y[length(aux_x$y)])
    }
    #calculate new points
    x_ini = 0
    y_ini = 1
    x_fim = 1
    y_fim = 0
    
    for(i in 1:(length(dref$x_ref)-1)){
      if(!is.na(dref$y[i])){
        x_ini = dref$x[i]
        y_ini = dref$y[i]
      }else{
        j=i+1
        while(is.na(dref$y[j])){
          j = j+1
        }
        x_fim = dref$x_ref[j]
        y_fim = dref$y[j]
        if(y_fim==y_ini){
          dref$y[i] = y_ini 
        }else{
          a = (y_fim-y_ini)/(x_fim-x_ini)
          b=y_fim-a*x_fim
          dref$y[i] = a*dref$x[i]+b
        }
        
      }
    }
    y_list[[length(y_list)+1]] = dref$y
    return(y_list)
  }
  
  
  tipe = unique(dat$type)
  
  dis = unique(dat$d1)
  dis = c("AOSD",dis)
  cv = unique(dat$cv)
  data = tibble()
  for(ti in tipe){
    for(d1 in 1:length(dis)){
      for(d2 in 1:length(dis)){
        if(d1>d2){
          print(paste0("mean ",ti," ",dis[d1]," ",dis[d2]))
          aux_b = dat[dat$type==ti & dat$d1==dis[d1] & dat$d2==dis[d2],]
          #create x_ref
          x = sort(unique(aux_b$sensitivities))
          x_ref = list()
          x_ref = append(x_ref,0)
          x_ref = append(x_ref,0)
          for(i in 2:(length(x)-1)){
            x_ref = append(x_ref,x[i])
            x_ref = append(x_ref,x[i])
          }
          x_ref = append(x_ref,1)
          x_ref = append(x_ref,1)
          x_ref = unlist(x_ref)
          y_mean = rep(NA,length(x_ref))
          y_cl = y_mean
          y_ch = y_mean
          y_list = list()
          auc_list = list()
          cv = unique(aux_b$cv)
          for(icv in cv){
            aux0 = aux_b[aux_b$cv==icv,]
            auc_list = append(auc_list,aux0$auc[1])
            aux0 = aux0[order(aux0$specificities,decreasing = TRUE),]
            aux0 = aux0[order(aux0$sensitivities),]
            x = aux0$sensitivities
            y = aux0$specificities
            y_list = interpolation_fun(x,y,x_ref,y_list)
          }
          auc_list = unlist(auc_list)
          auc_list = auc_list[!is.na(auc_list)]
          auc_mean = mean(auc_list)
          print(auc_mean)
          #x  +/-  z*(s/√n)
          inter = 1.96*(sd(auc_list)/sqrt(length(auc_list)))
          auc_cl = auc_mean - inter
          auc_ch = auc_mean + inter
          for(i in 1:length(y_mean)){
            resul = rep(NA,length(cv))
            
            for(icv in 1:length(y_list)){
              resul[icv] = y_list[[icv]][i]
            }
            y_mean[i] = mean(resul)
            inter = 1.96*(sd(resul)/sqrt(length(resul)))
            y_cl[i] = y_mean[i] - inter
            y_ch[i] = y_mean[i] + inter
            
          }
          x_ref2 = sort(unique(x_ref))
          y_cl2 = rep(NA,length(x_ref2))
          y_ch2 = rep(NA,length(x_ref2))
          y_mean2 = rep(NA,length(x_ref2))
          for(v in 1:length(x_ref2)){
            y_mean2[v] = mean(y_mean[x_ref==x_ref2[v]])
            y_cl2[v] = mean(y_cl[x_ref==x_ref2[v]])
            y_ch2[v] = mean(y_ch[x_ref==x_ref2[v]])
          }
          y_mean2[1] = 1
          y_cl2[1] = 1
          y_ch2[1] = 1
          tam = length(y_mean2)
          y_mean2[tam] = 0
          y_cl2[tam] = 0
          y_ch2[tam] = 0
          aux_t =data.frame(sensitivities = x_ref2,specificities = y_mean2, cl = y_cl2,ch=y_ch2)
          aux_t$ch[aux_t$ch>1] = 1
          aux_t$ch[aux_t$cl<0] = 0
          aux_t$auc_mean = auc_mean
          aux_t$auc_cl=auc_cl
          aux_t$auc_ch = auc_ch
          aux_t$type = ti
          aux_t$d1 = dis[d1]
          aux_t$d2 = dis[d2]
          data = rbind(data,aux_t)
          aux_t$d2 = dis[d1]
          aux_t$d1 = dis[d2]
          data = rbind(data,aux_t)
        }
      }
    }
    
  }
  return(data)
}

