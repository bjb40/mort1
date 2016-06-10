#dev R 3.2.1 "World-Famous Astronaut"
#Bryce Bartlett
#uses data downloaded from https://github.com/bjb40/cause_death

#generals
indir = 'H:/projects/cause_death/output/'
outdir = 'H:/Academic Projects/Mortality I/output/'


#input and clean

icd9=read.csv(paste0(indir,'icd9.csv'))
icd10=read.csv(paste0(indir,'icd10.csv'))


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Output of Tables and Histograms
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

sink(paste0(outdir,'codes_compare.txt'))


#table and histogram
#all types
print(table(icd9['type']))
print(table(icd10[icd10['mc_only']==0,'type']))


#3 code types
icd9_3 = table(icd9[,'div_no'],icd9[,'type'])
icd10_3 = table(icd10[icd10['mc_only']==0,'code3'],icd10[icd10['mc_only']==0,'type'])

cat('\n\nthree code types: External, Residual, acute, chronic\n')

cat('\nICD9\n')
cat(sum(icd9_3[,'External']>0),
sum(icd9_3[,'Residual']>0),
sum(icd9_3[,'Acute']>0),
sum(icd9_3[,'Chronic']>0))

cat('\nICD10\n')
cat(sum(icd10_3[,'External']>0),
sum(icd10_3[,'Residual']>0),
sum(icd10_3[,'Acute']>0),
sum(icd10_3[,'Chronic']>0))

#divisions (ICD10 only)
icd10_div = table(icd10[icd10['mc_only']==0,'div_nos'],icd10[icd10['mc_only']==0,'type'])

cat('\n\nICD10 divisions only\n')
cat(sum(icd10_div[,'External']>0),
sum(icd10_div[,'Residual']>0),
sum(icd10_div[,'Acute']>0),
sum(icd10_div[,'Chronic']>0))


#nchs113 tables
cat('\n\nNCHS 113 tables')
cat('\nICD9')
print(table(icd9[,'nchs113'],icd9[,'type']))

cat('\n\nICD10')
print(table(icd10[icd10['mc_only']==0,'nchs113'],icd10[icd10['mc_only']==0,'type']),width=150)



#chapter tables
cat('\n\nChapter tables')
cat('\nICD9')
print(table(icd9[,'chapter'],icd9[,'type']))

cat('\n\nICD10')
print(table(icd10[icd10['mc_only']==0,'chapter'],icd10[icd10['mc_only']==0,'type']),width=150)

sink()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#final table (revised from output above)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

nm = c('Acute','Chronic','External','Residual')
icd9_four = c(1355,1026,624+30,1517) #motor vehicle accidents
icd10_four = c(2537,1580+18,1296,2612)
icd9_three = c(202,172,94,217)
icd10_three = c(502,312+2,359,460)

#standardized proportions--THIS ONE -- NEED TO ADD 30 EXTERNAL FOR MV ACCIDENTS
icd9_four/sum(icd9_four)
icd10_four/sum(icd10_four)

icd9_four
icd10_four

icd9_three
icd10_three


icd9_three/sum(icd9_three)
icd10_three/sum(icd10_three)


#standardized by hypothetically available
#icd9_four/11000
#icd10_four/26000

#for (i in 1:4){
#  print(i)
#  print(prop.test(c(icd9_four[i],icd10_four[i]),c(11000,26000)))
#}

#icd9_three/1100
#icd10_three/2600

#for (i in 1:4){
#  print(i)
#  print(prop.test(c(icd9_three[i],icd10_three[i]),c(1100,2600)))
#}


#@@@@@@@@@@@@@@@@@@@@@
#Make table
#@@@@@@@@@@@@@@@@@@@@@

sink(paste0(outdir,'table2.md'))

cat('\nTable__. Summary of Causes of Death by category for ICD9 and ICD10.\n')

cat('\n|\t|',paste(nm,'|'),'Total|')
cat('\n|:----|-----:|-----:|-----:|-----:|------:|')
cat('\n|Three Character Codes|\t|\t|\t|\t|\t|')
cat('\n| Icd9 | ',paste(icd9_three,'|'),sum(icd9_three),'|')
cat('\n|\t|',paste0('(',round(prop.table(icd9_three),3),') |'),'(1.000) |')
cat('\n| icd10 |',paste(icd10_three,'|'),sum(icd10_three),'|')
cat('\n|\t|',paste0('(',round(prop.table(icd10_three),3),') |'),'(1.000) |')

cat('\n|Four Character Codes|\t|\t|\t|\t|\t|')
cat('\n| Icd9 | ',paste(icd9_four,'|'),sum(icd9_four),'|')
cat('\n|\t|',paste0('(',round(prop.table(icd9_four),3),') |'),'(1.000) |')
cat('\n| icd10 |',paste(icd10_four,'|'),sum(icd10_four),'|')
cat('\n|\t|',paste0('(',round(prop.table(icd10_four),3),') |'),'(1.000) |')

cat('\n\nProportions in parenthesis.')

cat('\n\nSignificance Tests:\n\n')

cat('\n\nThree code ICD10 versus ICD9:\n\n')
for(i in 1:4){
  cat('\n- ',nm[i],":\n")
  t = prop.test(c(icd9_three[i],icd10_three[i]),c(sum(icd9_three),sum(icd10_three)))
  cat('Est (9,10): ', round(t$estimate,3), ' pval: ', round(t$p.value,3))
}

cat('\n\nFour code ICD10 versus ICD9:\n\n')
for(i in 1:4){
  cat('\n- ',nm[i],":\n")
  t = prop.test(c(icd9_four[i],icd10_four[i]),c(sum(icd9_four),sum(icd10_four)))
  cat('Est (9,10): ', round(t$estimate,3), ' pval: ', round(t$p.value,3))
}


cat('\n\nFour code ICD9 versus three code icd9:\n\n')
for(i in 1:4){
  cat('\n- ',nm[i],":\n")
  t = prop.test(c(icd9_four[i],icd9_three[i]),c(sum(icd9_four),sum(icd9_three)))
  cat('Est (three,four): ', round(t$estimate,3), ' pval: ', round(t$p.value,3))
}

cat('\n\nFour code ICD10 versus three code icd10:\n\n')
for(i in 1:4){
  cat('\n- ',nm[i],":\n")
  t = prop.test(c(icd10_four[i],icd10_three[i]),c(sum(icd10_four),sum(icd10_three)))
  cat('Est (three,four): ', round(t$estimate,3), ' pval: ', round(t$p.value,3))
}

sink()

#@@@@@@@@@@@@@@@@@@@@@
#Make supplement table
#@@@@@@@@@@@@@@@@@@@@@

#group, 113, 113 ti, icd9 numbers, icd10 numbers

#View(icd10$code[icd10$nchsti == 'All other diseases (Residual)'])

lim = (icd10$nchs113 == 113)

as.character(unique(icd10$code3ti[lim]))



