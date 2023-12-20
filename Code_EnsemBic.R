# Check the path whether the input files are in same directory or not.
getwd()

#### Step 1 and Step 2 #####

# This program needs two input files "Total_Biclusters.csv" and "Ensemble Bicluster.csv"
# Initially create two empty .csv files, viz. "Total_Biclusters.csv" and "Ensemble Bicluster.csv"

# To generate the contents of "Total_Biclusters.csv" fetch the dataset into BiclustGui 3.0 tool.
# The parameter settings of BiclustGui is provided in Table 3 of the manuscript.

# Execute different biclustering algorithms using same dataset.
# Store the identified biclusters in  "Total_Biclusters.csv"


#### Step 3 ####

# Compute the p-value for all the biclusters one by one using the web tool FuncAssociate 3.0
# Parameter settings of FuncAssociate 3.0 are “Species: Homo sapiens” and “Namespace: hgnc_symbol"

# While checking the p-value, manually delete the the biclusters from "Total_Biclusters.csv" whose p>0.05
# For the remaining biclusters keep the p-value in a file "p-value.csv" in the same order the biclusters are stored in "Total_Biclusters.csv"

# Now, both the input files are ready to fetch and execute. 

Tot_Bic = read.csv("Total_Biclusters.csv")
p-value = read.csv("p-value.csv")

r = ncol(Tot_Bic)
# Each column of Tot_Bic depicts individual bicluster, i.e., the list of genes of that bicluster. 

# Declare one data frame to keep the ensemble bicluster
Ens_Bic = data.frame()

p_value = data.frame()
p_value[1,1] = 7
p_value[1,2] = 4
p_value[1,3] = 3
p_value[1,4] = 6

if(r > 0)
{
  a<-which.min(p_value)
  
  for(i in 1:nrow(Tot_Bic))
  {
    Ens_Bic[i,1] = Tot_Bic[i,a]
    Tot_Bic[i,a] = ""
  }
  colnames(Ens_Bic) <- c("Ensemble_Bicluster")
  Tot_Bic = Tot_Bic[,-c(a)]
}

# Next, update the "Total_Biclusters.csv" file.
# Delete the genes present in "Ensemble Bicluster.csv"

New_Tot_Bic = Tot_Bic

nc_Tot = ncol(Tot_Bic)
nr_Tot = nrow(Tot_Bic)
nc_Ens = ncol(Ens_Bic)
nr_Ens = nrow(Ens_Bic)

for(j in 1:nc_Ens)
{
  for(i in 1:nr_Ens)
  {
    for(n in 1:nc_Tot)
    {
      for(m in 1:nr_Tot)
      {
        if(Tot_Bic[m,n] == "")
        {
          break
        }
        
        if(Ens_Bic[i,j] == Tot_Bic[m,n])
        {
          New_Tot_Bic[m,n] = NA
          break
        }
        #        print(paste("m = ",m, "n = ",n))
      }
    }
    print(paste("i = ", i, "J = ", j))
  }
}

df = data.frame()
nr_Tot = nrow(New_Tot_Bic)
nc_Tot = ncol(New_Tot_Bic)

for(q in 1:nc_Tot)
{
  count = 1
  for(p in 1:nr_Tot)
  {
    if(is.na(New_Tot_Bic[p,q]) == FALSE)
    {
      if(New_Tot_Bic[p,q] == "")
      {
        break
      }
      df[count,q] = New_Tot_Bic[p,q]
      count = count + 1
    }
    print(paste("p=",p, "q = ",q))
  }
  colnames(df)[q] = colnames(New_Tot_Bic)[q]
}

df1 = df
df1[is.na(df1)] = ""
write.csv(df1, "")
