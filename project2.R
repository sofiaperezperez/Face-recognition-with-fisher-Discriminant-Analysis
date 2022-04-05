####### PART A############

FDA=function(Images=list.files(getwd())){
  
  #Implementation PCA for dimensionality reduction
  
  library(OpenImageR)
  vector=matrix(0, nrow=108000, ncol=length(Images))
  
  #108000 is the number of pixels
  
  #creation of the labels
  
  people=rep(0,length(Images))
  for(i in 1:length(Images)){
    if (nchar(Images[i])==8){
      people[i]=substr(Images[i],0,2)
    }
    else{
      people [i]=substr(Images[i],0,1)
      
    }
    
  }
  
  
  
  #Transformation of the iamges into vectors, and then
  
  #all of them contatenated crating a matrix
  
  for (i in 1:length(Images)){
    Im=readImage(Images[i])
    red=Im[,,1]
    green=Im[,,2]
    blue=Im[,,3]
    
    
    data=NULL
    data=cbind(data,as.vector(red)) 
    data=cbind(data,as.vector(green))
    data=cbind(data,as.vector(blue))
    vector[,i]=c(data[,1],data[,2],data[,3]) #Creation of an image into a vector
    
  }
  
  
  
  data_original=t(as.data.frame(vector))#people as rows
  mean_train=apply(data_original,2,mean)
  #we standarize our data substracting the mean
  
  for (i in 1:nrow(data_original)){
    data_original[i,]=data_original[i,]-mean_train
  }
  
  data_transposed=t(data_original) #people as cols
  
  
  #Since the number of observations is much bigger than the number of variables:
  G=data_original
  small=(G%*%t(G))/(nrow(data_original-1))
  
  
  # Calculation of the eigenvectors and eigenvalues of the small matrix.
  Eigen = eigen(small)
  Eigenvalues = Eigen$values
  
  #the variance explained by each of the principal components.
  Cumulative.Var = cumsum(Eigenvalues)/sum(Eigenvalues)
  Cumulative.Var=round(Cumulative.Var,4)
  Eigenvectors = Eigen$vectors  #our eigenvectors are our new principal axis
  
  
  #P is going to be the projection matrix, that will multiply the original
  
  #data in order to obtain the projected one.
  P=as.matrix(t(G))%*%as.matrix(Eigenvectors[,1:150]) #matrix to multiply the people
  
  newdata=as.matrix(data_original)%*%as.matrix(P)  #people projected in the new variables
  
  #Selection of 24 first principal components, retaining 95% of the original variance
  
  newdata=data.frame(newdata=newdata[,1:24], labels=as.integer(people))
  
  colnames(newdata)=c("PC.1","PC.2","PC.3","PC.4","PC.5","PC.6","PC.7","PC.8","PC.9","PC.10",
                      
                      "PC.11","PC.12","PC.13","PC.14","PC.15","PC.16","PC.17","PC.18","PC.19","PC.20",
                      
                      "PC.21","PC.22","PC.23","PC.24","labels")
  
  
  #NOW WE START WITH FISHER DISCRIMINANT ANALYSIS
  #computation mean vectors
  m = colMeans(newdata[,-25])
  mean_vectors=matrix(0,ncol=ncol(newdata)-1,nrow= 25)
  
  for (i in 1:24){
    mean_vectors[i,]=colMeans(newdata[newdata$labels==as.character(i),1:24])
    
  }
  
  
  #Compute the scatter matrices
  SW=0
  for (i in 1:25){
    SW=SW+cov(newdata[which(newdata$labels==as.character(i)),-25])*(table(newdata$labels)[i]-1)
    
  }
  
  #Compute the distance between classes
  S.B=0
  for ( i in 1:25){
    S.B=S.B+ (table(newdata$labels)[i])*(mean_vectors[i,]-m)%*%t(mean_vectors[i,]-m)
  }
  
  #There are a lot of projections, but we want to obtain the "best" one
  
  #compute the eigenvectors
  eig=eigen((solve(SW)%*%S.B))
  eigenvalues=eig$values
  eigenvectors=eig$vectors  #if we choose 10, 95% of the variance retained
  prop.var=eigenvalues/sum(eigenvalues)
  cummulative.var=cumsum(eig$values/sum(eig$values))
  
  
  #As we can see, if we keep adding more eigenvectors, the %of variance explained is bigger.
  
  #However, there is a point (we have considered 10), in which is not really significant to keep adding
  
  library(ggplot2)
  ggplot(data.frame(cummulative.var))+aes(x=1:24,y=cummulative.var)+geom_point()+geom_line()
  
  
  
  #cumputation new dataset
  data.new=as.matrix(newdata[,1:24])%*%as.matrix(eigenvectors[,1:10])
  data.new=data.frame(data.new)
  colnames(data.new)=c("PC.1","PC.2","PC.3","PC.4","PC.5","PC.6","PC.7","PC.8","PC.9","PC.10")
  data.new=data.frame(data.new,labels=people)
  
  
  
  library(ggplot2)
  ggplot()+aes(data.new$PC.1,color=data.new$labels,fill=data.new$labels)+
    geom_histogram(alpha=0.5,position = "identity",bins=300)+
    theme(legend.position = "bottom")
  
  #As we can see, most the data is projected into one line (best projection to minimice
  
  #the inside variance and to maximice the distance between points of different classes)
  
  return(list(eigenvectors_fisher=eigenvectors,mean=mean_train,D=eigenvalues,P=P,Var=Cumulative.Var,
              train=data.new,raw_data=newdata,people=people))
  
}


partA=FDA()
mean_train=partA$mean
data.new=partA$train
eigenvectors_fisher=partA$eigenvectors_fisher
data_without_fisher=partA$raw_data
P=partA$P
people=partA$people
########DISTANCE MATRIX

#matrix that stores the distance of one photo with the rest.It is used to compute the threshold.

matrix_dist=matrix(0, nrow=nrow(data.new), ncol=nrow(data.new))

for (j in 1:nrow(data.new)){
  
  for( i in 1:nrow(data.new)){
    
    matrix_dist[i,j]=sqrt(sum((data.new[i,-11] - data.new[j,-11])^2))
    
  }
  
}



ggplot()+aes(x=matrix_dist,color="red",fill="red")+geom_histogram(bins = 25)+
  geom_vline(xintercept = 1600,color="black")+
  xlab("distances")+labs("distances matrix")



#THERSHOLD=1700

#####KNN

#function to predict

KNN=function(data_tr,data_tst,train_labels,k,distance){
  
  predictions=rep(0,nrow(data_tst))
  
  for (j in 1:nrow(data_tst)){ #aqui predices cada row del test (25 predicciones)
    dmatrix=dist(rbind(data_tst[j,],data_tr), method = distance, diag = TRUE, upper = TRUE)
    dmatrix=as.matrix(dmatrix)
    dmatrix=dmatrix[1,2:(nrow(data_tr)+1)]
    ordenados=sort(dmatrix,index.return=TRUE,decreasing=FALSE) #ordea las distancias de menor a mayor
    labels_sel=as.character(train_labels[ordenados$ix[1:k]])
    uniqv <- unique(labels_sel)
    predictions[j]=uniqv[which.max(tabulate(match(labels_sel, uniqv)))]
    }
  
  return(predictions)
}



#TRAINING OF THE ALGORITHM

nfolds=5
n=nrow(data.new)
folds=sample(cut(1:n,breaks=nfolds,labels=FALSE))

best_accs=NULL #aqui vas a almacenar la mejor acc de cada fold



k=c(3,5,7)
fisher=c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
distances=c("euclidean","manhattan","minkowski")
params=expand.grid(k=k,dist=distances,fisher=fisher)
result=data.frame(params,acc=0)

length(fisher)

for (i in 1:nfolds){
  
  for (j in 2:length(fisher)){
    data.new=as.matrix(data_without_fisher[,1:24])%*%as.matrix(eigenvectors_fisher[,1:j])
    data.new=data.frame(data.new)
    data.new=data.frame(data.new,labels=people)
    
    data_tr=data.new[which(folds!=i),-ncol(data.new)]
    data_tst=data.new[which(folds==i),-ncol(data.new)]
    train_labels=data.new[folds!=i,ncol(data.new)]
    test_labels=data.new[folds==i,ncol(data.new)]
    
    result$acc=0 
    
    for (p in 1:nrow(result)){ #15 combinaciones
      
      final=KNN(data_tr,data_tst,train_labels,result$k[p],result$dist[p])
      result$acc[p]= mean(final==test_labels) #esto te da la accuracy jeje
      
    }
  }
  
  best=result[which.max(result$acc),] #por cada divison del data una accuracy
  best_accs=rbind(best_accs,best)
  
}



opt_params=best_accs[which.max(best_accs$acc),]  #the combination that provides the best accuracy

#optimal combination of parameters

k=opt_params$k
distance=opt_params$dist
n_cols=opt_params$fisher

best_data=as.matrix(data_without_fisher[,1:24])%*%as.matrix(eigenvectors_fisher[,1:n_cols])
best_data=data.frame(best_data)
best_data=data.frame(best_data,labels=people)

save(best_data,P,eigenvectors_fisher,n_cols=n_cols,mean_train,k,distance,KNN,file ="data_files_P2.RData" ) #esto luego cuando tengamos la funcion en la parte a)



#######PART B########

classifier=function(Image){
  load("data_files_P2.RData")
  
  library(stringr)
  library(OpenImageR)
  
  #transformation of the image into a vector
  Image=readImage(Image)
  red=Image[,,1]
  green=Image[,,2]
  blue=Image[,,3]
  
  data_image=NULL
  data_image=cbind(data_image,as.vector(red)) 
  data_image=cbind(data_image,as.vector(green))
  data_image=cbind(data_image,as.vector(blue))
  
  new_image=c(data_image[,1],data_image[,2],data_image[,3])
  
  
  #projection of the image and standarization
  
  new_image=new_image-mean_train
  new_vector=as.vector(new_image)%*%as.matrix(P[,1:24])  #projection into the new axis
  new_vector=new_vector%*%as.matrix(eigenvectors_fisher[,1:n_cols]) #projection into fisher axis
  colnames(new_vector)=colnames(best_data[,-ncol(best_data)])
  
  #now we check whether it is a person from our dataset
  prediction=KNN(data_tr=best_data[,-ncol(best_data)],data_tst=new_vector,train_labels = as.character(best_data$labels),k=k,distance = distance)
  
  #now we check the threshold
  
  min_distance=as.matrix(dist(rbind(new_vector,best_data[which(prediction==(best_data$labels)),-ncol(best_data)]),method="euclidean"))
  if (min_distance[1,1]>=1700){
    
    print(-1)
  }
  
  else{
    return(str_c("The person belongs to the data set, corresponds to person number ",prediction))
    
  }
  
}

classifier("3DT.jpg")
classifier("11AT.jpg")
