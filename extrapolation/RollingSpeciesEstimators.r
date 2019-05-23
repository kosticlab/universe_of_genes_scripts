## install the latest version from github
library(devtools)

## import packages
library(iNEXT)
library(ggplot2)
library(PopED)

library(SPECIES)
library(optimbase)

Root = './'
File = 'OccurenceMatrix.csv'

NumSamples = 1473
AA = read.csv(paste(Root,File,sep = ''), header = FALSE)
MaxTons = size(AA)[1] #The maximum number of times the gene has been found in the sampling

Chao2 = zeros(NumSamples,4)
ChaoBungeData = zeros(NumSamples,4)
ChaoLeeACE = zeros(NumSamples,4)
ChaoLeeACE2 = zeros(NumSamples,4)
JackknifeData = zeros(NumSamples,4)


for(i in 1:NumSamples){
  Result = chao1984(cbind(1:MaxTons, AA[i]))
  Chao2[i,] = c(Result$Nhat, Result$SE, Result$CI[1], Result$CI[2])}
write.csv(x = Chao2, file = paste(Root,'95Chao2.csv',sep=''), sep = ',')

for(i in 1:NumSamples){
  Result = ChaoLee1992(cbind(1:MaxTons,AA[i]))
  ChaoLeeACE[i,] = c(Result$Nhat[1], Result$SE[1], Result$CI[1,1], Result$CI[1,2])
  ChaoLeeACE2[i,] = c(Result$Nhat[2], Result$SE[2], Result$CI[2,1], Result$CI[2,2])
}
write.csv(x = ChaoLeeACE, file = paste(Root,'95ChaoLeeACE.csv',sep=''), sep = ',')
write.csv(x = ChaoLeeACE2, file = paste(Root,'95ChaoLeeACE2.csv', sep=''), sep = ',')

for(i in 1:NumSamples){
  Result = ChaoBunge(cbind(1:MaxTons, AA[i]), t = 9)
  ChaoBungeData[i,] = c(Result$Nhat, Result$SE, Result$CI[1], Result$CI[2])}
write.csv(x = ChaoBungeData, file = paste(Root,'95ChaoBunge9.csv',sep=''), sep = ',')

for(i in 1:NumSamples){
  Result = ChaoBunge(cbind(1:MaxTons, AA[i]), t = 6)
  ChaoBungeData[i,] = c(Result$Nhat, Result$SE, Result$CI[1], Result$CI[2])}
write.csv(x = ChaoBungeData, file = paste(Root,'95ChaoBunge6.csv',sep=''), sep = ',')

for(i in 1:NumSamples){
  Result = ChaoBunge(cbind(1:MaxTons, AA[i]), t = 3)
  ChaoBungeData[i,] = c(Result$Nhat, Result$SE, Result$CI[1], Result$CI[2])}
write.csv(x = ChaoBungeData, file = paste(Root,'95ChaoBunge3.csv',sep=''), sep = ',')

for(i in 2:NumSamples){
  Result = jackknife(cbind(1:MaxTons, AA[i]),k=6)
  JackknifeData[i,] = c(Result$Nhat, Result$SE, Result$CI[1], Result$CI[2])}
write.csv(x = JackknifeData, file = paste(Root,'95Jackknife6.csv', sep=''), sep = ',')

for(i in 2:NumSamples){
  Result = jackknife(cbind(1:MaxTons, AA[i]),k=3)
  JackknifeData[i,] = c(Result$Nhat, Result$SE, Result$CI[1], Result$CI[2])}
write.csv(x = JackknifeData, file = paste(Root,'95Jackknife3.csv', sep=''), sep = ',')

for(i in 2:NumSamples){
  Result = jackknife(cbind(1:MaxTons, AA[i]),k=9)
  JackknifeData[i,] = c(Result$Nhat, Result$SE, Result$CI[1], Result$CI[2])}
write.csv(x = JackknifeData, file = paste(Root,'95Jackknife9.csv', sep = ''), sep = ',')
