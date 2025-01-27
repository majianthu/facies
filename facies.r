library(copent)
library(randomForest)

is.adjfacies<-function(a,b){
  r = 0
  if(a==1){
    if(b==1|b==2){r=1}
  }
  if(a==2){
    if(b==1|b==2|b==3){r=1}
  }
  if(a==3){
    if(b==2|b==3){r=1}
  }
  if(a==4){
    if(b==4|b==5){r=1}
  }
  if(a==5){
    if(b==4|b==5|b==6){r=1}
  }
  if(a==6){
    if(b==5|b==6|b==7){r=1}
  }
  if(a==7){
    if(b==6|b==7|b==8){r=1}
  }
  if(a==8){
    if(b==6|b==7|b==8|b==9){r=1}
  }
  if(a==9){
    if(b==7|b==8|b==9){r=1}
  }
  r
}

nadj<-function(x1,x2){
  n = 0
  for(i in 1:length(x1)){
    n = n + is.adjfacies(x1[i],x2[i])
  }
  n
}

tr1 = read.csv("~/Rworks/facies/facies.csv",header = T)

x11();
pairs(tr1[,c(1,4:11)], col = tr1[,1], pch = ".")

x11();
c1 = hist(tr1[,1],breaks = 0:9)$counts
names(c1) = c("SS","CSiS","FSiS","SiSh","MS","WS","D","PS","BS")
at0 = barplot(c1,col = 1:9, ylab = "frequency", ylim = c(0,max(c1)+30))
text(x = at0, y = c1+15, labels = c1, col = "red")

x11(height = 10, width = 10);
par(mfrow = c(1,5))
idx1 = which(tr1$Well.Name=="SHRIMPLIN")
plot(tr1[idx1,5],tr1[idx1,4], col = tr1[idx1,1], ylim = rev(range(tr1[idx1,4])), 
     xlab = "GR", ylab = "Depth", type = "b", pch = 20)
legend(x = 200, y = min(tr1[idx1,4]), names(c1), col = 1:9, pch = 20)
plot(tr1[idx1,6],tr1[idx1,4], col = tr1[idx1,1], ylim = rev(range(tr1[idx1,4])), 
     xlab = "ILD_log10", ylab = "", type = "b", pch = 20)
plot(tr1[idx1,7],tr1[idx1,4], col = tr1[idx1,1], ylim = rev(range(tr1[idx1,4])), 
     xlab = "DeltaPHI", ylab = "", type = "b", pch = 20)
plot(tr1[idx1,8],tr1[idx1,4], col = tr1[idx1,1], ylim = rev(range(tr1[idx1,4])), 
     xlab = "PHIND", ylab = "", type = "b", pch = 20)
plot(tr1[idx1,9],tr1[idx1,4], col = tr1[idx1,1], ylim = rev(range(tr1[idx1,4])), 
     xlab = "PE", ylab = "", type = "b", pch = 20)

ce1 = rep(0,8)
for(i in 1:8){
  ce1[i] = copent(tr1[,c(1,i+3)])
}
x11(); 
at1 = barplot(ce1, xaxt = "n", ylab = "CE")
text(x = at1, y = -0.01, srt = 60, adj = 1, xpd = TRUE, labels = names(tr1)[4:11], col = 'red')

k1 = 30
mae1 = mae2 = mae3 = mae4 = mae5 = mae6 = mae7 = rep(0,k1)
acc1 = acc2 = acc3 = acc4 = acc5 = acc6 = acc7 = rep(0,k1)
adj1 = adj2 = adj3 = adj4 = adj5 = adj6 = adj7 = rep(0,k1)
n = dim(tr1)[1]; n1 = round(n*0.8)
for(k in 1:k1){
  r1 = rank(runif(n))
  y = tr1[r1[1:n1],1]
  x1 = tr1[r1[1:n1],4:11]
  x2 = tr1[r1[1:n1],4:10]
  x3 = tr1[r1[1:n1],c(4:6,8:10)]
  x4 = tr1[r1[1:n1],c(4:6,9,10)]
  x5 = tr1[r1[1:n1],c(4,5,9,10)]
  x6 = tr1[r1[1:n1],c(4,9,10)]
  x7 = tr1[r1[1:n1],c(4,10)]
  
  yt = tr1[r1[(n1+1):n],1]
  x1t = tr1[r1[(n1+1):n],4:11]
  x2t = tr1[r1[(n1+1):n],4:10]
  x3t = tr1[r1[(n1+1):n],c(4:6,8:10)]
  x4t = tr1[r1[(n1+1):n],c(4:6,9,10)]
  x5t = tr1[r1[(n1+1):n],c(4,5,9,10)]
  x6t = tr1[r1[(n1+1):n],c(4,9,10)]
  x7t = tr1[r1[(n1+1):n],c(4,10)]

  rf1 = randomForest(x1,y)
  rf2 = randomForest(x2,y)
  rf3 = randomForest(x3,y)
  rf4 = randomForest(x4,y)
  rf5 = randomForest(x5,y)
  rf6 = randomForest(x6,y)
  rf7 = randomForest(x7,y)
  
  pred1 = predict(rf1,x1t)
  pred2 = predict(rf2,x2t)
  pred3 = predict(rf3,x3t)
  pred4 = predict(rf4,x4t)
  pred5 = predict(rf5,x5t)
  pred6 = predict(rf6,x6t)
  pred7 = predict(rf7,x7t)
  
  mae1[k] = mean(abs(yt-round(pred1)))
  mae2[k] = mean(abs(yt-round(pred2)))
  mae3[k] = mean(abs(yt-round(pred3)))
  mae4[k] = mean(abs(yt-round(pred4)))
  mae5[k] = mean(abs(yt-round(pred5)))
  mae6[k] = mean(abs(yt-round(pred6)))
  mae7[k] = mean(abs(yt-round(pred7)))
  
  acc1[k] = length(which(yt==round(pred1)))
  acc2[k] = length(which(yt==round(pred2)))
  acc3[k] = length(which(yt==round(pred3)))
  acc4[k] = length(which(yt==round(pred4)))
  acc5[k] = length(which(yt==round(pred5)))
  acc6[k] = length(which(yt==round(pred6)))
  acc7[k] = length(which(yt==round(pred7)))
  
  adj1[k] = nadj(yt,round(pred1))
  adj2[k] = nadj(yt,round(pred2))
  adj3[k] = nadj(yt,round(pred3))
  adj4[k] = nadj(yt,round(pred4))
  adj5[k] = nadj(yt,round(pred5))
  adj6[k] = nadj(yt,round(pred6))
  adj7[k] = nadj(yt,round(pred7))
}
mae = c(mean(mae1),mean(mae2),mean(mae3),mean(mae4),mean(mae5),mean(mae6),mean(mae7))
names(mae) = 8:2
acc = c(mean(acc1),mean(acc2),mean(acc3),mean(acc4),mean(acc5),mean(acc6),mean(acc7)) / (n-n1)
names(acc) = 8:2
adj = c(mean(adj1),mean(adj2),mean(adj3),mean(adj4),mean(adj5),mean(adj6),mean(adj7)) / (n-n1)
names(adj) = 8:2

x11(); 
at2 = barplot(mae, ylim = c(0,max(mae)+0.2), ylab = "MAE")
text(x = at2, y = mae + 0.02, labels = round(mae,4), col = 'red')

x11(); 
at3 = barplot(acc, ylim = c(0,max(acc)+0.05), ylab = "Accuracy")
text(x = at3, y = acc + 0.015, labels = round(acc,4), col = 'red')

x11(); 
at4 = barplot(adj, ylim = c(0,max(adj)+0.05), ylab = "Adjusted Accuracy")
text(x = at4, y = adj + 0.015, labels = round(adj,4), col = 'red')
