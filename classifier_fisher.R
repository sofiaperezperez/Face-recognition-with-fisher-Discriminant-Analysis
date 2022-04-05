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

classifier("1AT.jpg")
