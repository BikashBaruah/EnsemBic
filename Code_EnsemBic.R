############### start #################

# Check the path whether the input files are in same directory or not.
getwd()

#### Step 1 and Step 2 #####

# This program needs two input files "Total_Biclusters.csv" and "p-value.csv"
# Initially create two empty .csv files (manually), viz. "Total_Biclusters.csv" and "p-value.csv"

# To generate the contents of "Total_Biclusters.csv" fetch the dataset into "BiclustGui 3.0" tool.
# The parameter settings of "BiclustGui 3.0" is provided in Table 3 of the manuscript entitled 
# "EnsemBic: An Effective Ensemble of Biclustering to Identify Potential Biomarkers of Esophageal Squamous Cell Carcinoma".

# Execute different biclustering algorithms using same dataset.
# Store the identified biclusters in  "Total_Biclusters.csv"


#### Step 3 ####

# Compute the p-value for all the biclusters one by one using the web tool "FuncAssociate 3.0"
# Parameter settings of "FuncAssociate 3.0" are “Species: Homo sapiens” and “Namespace: hgnc_symbol"
# "FuncAssociate 3.0" parameter settings are also available in the above mentioned manuscript.
 
# While checking the p-value, manually delete the the biclusters from "Total_Biclusters.csv" whose p>0.05
# For the remaining biclusters keep the p-value in a file "p-value.csv" in the same order the biclusters are stored in "Total_Biclusters.csv"

# Now, both the input files are ready to fetch and execute. 

Tot_Bic = read.csv("Total_Biclusters.csv")
p-value = read.csv("p-value.csv")

r = ncol(Tot_Bic)
# Each column of "Tot_Bic" depicts a list of genes comprising a bicluster. 

# Reject the biclusters whose p>0.005
for(i in r:1)
{
  if(p_value[1,i]>0.05)
  {
    Tot_Bic = Tot_Bic[,-c(i)]
    p_value = p_value[,-c(i)]
  }
  print(i)
}

# Update the value of r
r = ncol(Tot_Bic)

# Declare one data frame to keep the ensemble bicluster
Ens_Bic = data.frame()

# Find ensemble bicluster and keep it in "Ens_Bic"
if(r > 0)
{
  a<-which.min(p_value)
  
  for(i in 1:nrow(Tot_Bic))
  {
    Ens_Bic[i,1] = Tot_Bic[i,a]
    Tot_Bic[i,a] = ""
  }
  colnames(Ens_Bic) <- c("Ensemble_Bicluster")
  
  # Remove the ensemble bicluster from the "Tot_Bic"
  Tot_Bic = Tot_Bic[,-c(a)]
  
  # Decrement the value of "r" by 1
  r = r-1
}else
{
  stop("Since r value is less than 0; Execution is stopeed here")
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
        print(paste("m = ",m, "n = ",n))
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
Tot_Bic = df1

# write the updated "Total_Bicluster.csv" file 
write.csv(Tot_Bic, "Total_Biclusters.csv")

# Check again, whether "Tot_Bic" is empty or not
r = ncol(Tot_Bic)

# Check if r > 0 then again repeat the "Step 3" to get the next ensemble Bicluster
# For that, generate new "p_value.csv" using updated "Total_Bicluster.csv" file and "Func_Associate 3.0" 
# Otherwise no need of further execution, i.e., all the ensemble biclusters are already genereted.

######### end #########
