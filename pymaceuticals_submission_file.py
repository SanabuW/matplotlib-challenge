#!/usr/bin/env python
# coding: utf-8

# In[140]:


# Import dependencies
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import numpy as np

# Set data set directory strings to variables
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read dataset csv's into dataframe
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Begin data review
# Display mouse data
mouse_metadata


# In[141]:


# Display study data
study_results


# In[142]:


# Merge data to attribute Drug Regimen to Tumor Volume
mergeData = pd.merge(study_results, mouse_metadata, on = "Mouse ID", how = "left")
mergeData
#mergeData.sort_values(["Mouse ID", "Timepoint"])


# In[143]:


# Use .duplicated to retrieve Mouse IDs in rows where Mouse IDs and timepoints are a duplicate
mergeDataDupe = mergeData[mergeData.duplicated(subset = ["Mouse ID","Timepoint"], keep = False)]
# Get list of any Mouse ID that has been duplicated
dupeMice = mergeDataDupe["Mouse ID"].unique()


# In[144]:


# Use for loop to regenerate database without each Mouse ID found to be a duplicate
# *All data on the mouse will be removed, as instructed by the assignment instructions "...check the data 
# for any mouse ID with duplicate time points and remove any data associated with that mouse ID."
for x in dupeMice:
    mergeDataClean = mergeData[mergeData["Mouse ID"] != x]

# Recheck for duplicates
mergeDataCleanChk = mergeDataClean[mergeDataClean.duplicated(subset = ["Mouse ID","Timepoint"], keep = False)]
print(mergeDataCleanChk)


# In[145]:


# Review final cleaned dataframe
print(mergeDataClean)
print(f'-------------------------')
print(mergeDataClean.dtypes)
print(f'-------------------------')
print(mergeDataClean.count())
print(f'-------------------------')


# In[146]:


#Final Review
    # 1880 rows × 8 columns
    # No Null values
    # Data types OK
    # Values OK
    # No duplicate Mouse IDs with duplicate timepoints in mergeDataClean
# OK to proceed with calculations
# Use mergeDataClean as the dataset


# In[147]:


# Generate a summary statistics table consisting of the following stats on the tumor 
# volume for each drug regimen

# Separate the Tumor Volume and Drug Regimen
tumorVoldf = mergeDataClean.loc[:,["Tumor Volume (mm3)","Drug Regimen"]]

# Group by the drug regimen
drugTumorGroups = tumorVoldf.groupby("Drug Regimen")
# Find mean
drugMeanSer = drugTumorGroups.mean()
# Find median
drugMedSer = drugTumorGroups.median()
# Find variance
drugVarSer = drugTumorGroups.var(ddof = 0)
# Find standard deviation
drugSTDSer = drugTumorGroups.std(ddof = 0)
# Find SEM
drugSEMSer = drugTumorGroups.sem(ddof=0)

# Assemble summary table
summaryStatsdf = pd.DataFrame({"Mean":drugMeanSer.iloc[:,0], 
                             "Median":drugMedSer.iloc[:,0],
                             "Variance":drugVarSer.iloc[:,0],
                             "St Dev":drugSTDSer.iloc[:,0],
                             "SEM":drugSEMSer.iloc[:,0]
                            })
# Format summary table
summaryStatsdf["Mean"] = summaryStatsdf["Mean"].map("{:,.3f}".format)
summaryStatsdf["Median"] = summaryStatsdf["Median"].map("{:,.3f}".format)
summaryStatsdf["Variance"] = summaryStatsdf["Variance"].map("{:,.3f}".format)
summaryStatsdf["St Dev"] = summaryStatsdf["St Dev"].map("{:,.3f}".format)
summaryStatsdf["SEM"] = summaryStatsdf["SEM"].map("{:,.3f}".format)
summaryStatsdf


# In[148]:


# Create and display bar plot for total number of measurements for each regimen as instructed by the
# homework instructions.

# DataFrame.plot() version

# Get count values
drugGroups = mergeDataClean.groupby("Drug Regimen")
countValues = drugGroups["Drug Regimen"].count()

# Plot bar chart
drugMeasureChart = countValues.plot(kind = "bar", width = 0.8, title = "No. of Measurements for Each Drug Regimen")

# Format chart
drugMeasureChart.set_ylabel("No. of Measurements")
plt.xticks(rotation = 90)
plt.xlim(-1.25, len(np.arange(len(drugGroups))) + .25)
plt.ylim(0, max(countValues) + 40)
plt.legend()
plt.show()


# In[149]:


#pyplot version
# Get count values
xValues = countValues.index
yValues = countValues

# Plot bar chart
plt.bar(xValues, yValues, width=0.8)

# Format bar chart
plt.xticks(rotation = 90)
plt.title("No. of Measurements for Each Drug Regimen")
plt.xlabel("Mouse ID")
plt.ylabel("No. of Measurements")
plt.xlim(-1.25, len(np.arange(len(yValues.index))) + .25)
plt.ylim(0, max(yValues) + 40)
plt.legend([yValues.index.name], loc = "best")
plt.show()


# In[150]:


# Create and display pie plot for total number of mice per gender
# DataFrame.plot() version

# Get count values
genderGroups = mergeDataClean.groupby("Sex")
genderCountValues = genderGroups["Mouse ID"].nunique()
femaleCount = genderCountValues["Female"]
maleCount = genderCountValues["Male"]

# Plot pie chart
genderCountValues.plot(kind = "pie", autopct = "%1.1f%%", shadow = False,startangle = 90, 
                                     title = "Gender Distribution of Mice", colors =["lightblue", "lightskyblue"])
# Format pie chart
plt.axis("equal")
plt.xlabel("")
plt.ylabel("")
plt.show()


# In[151]:


# Create and display pie plot for total number of mice per gender
# pyplot version
# Plot pie chart
# Create labels for pyplot input
labels = list(genderCountValues.index.values) 
#labels
plt.pie(genderCountValues, labels = labels, autopct="%1.1f%%", shadow=False,startangle=90, 
        colors =["lightblue", "lightskyblue"])
# Format pie chart
plt.title("Gender Distribution of Mice")
plt.xlabel("")
plt.ylabel("")
plt.axis("equal")
plt.show()


# In[152]:


# Calculate the final tumor volume of each mouse across four of the most promising treatment regimens: 
# Capomulin, Ramicane, Infubinol, and Ceftamin. 
# Calculate the quartiles and IQR and quantitatively determine if there are any potential 
# outliers across all four treatment regimens.
topTreatmentsdf = mergeDataClean.loc[(mergeDataClean["Drug Regimen"] == "Capomulin") | 
                                      (mergeDataClean["Drug Regimen"] == "Ramicane") | 
                                      (mergeDataClean["Drug Regimen"] == "Infubinol") | 
                                      (mergeDataClean["Drug Regimen"] == "Ceftamin")
                                     ,:]
# Find the final timepoint for all mice
topTreatmentGroup = topTreatmentsdf.groupby("Mouse ID")
maxTime = topTreatmentGroup["Timepoint"].max()
# Attach the Tumor Volumes to the found Mice and Timepoints by inner joining on both series to "filter out"
# all other entries
finalTumordf = pd.merge(maxTime, mergeDataClean, on=["Mouse ID", "Timepoint"], how="inner")
# Rename tumor volume series to reflect intent of the data
finalTumordf = finalTumordf.rename(columns = {"Tumor Volume (mm3)":"Final Tumor Volume"})

# Review merged dataframe
finalTumordf


# In[164]:


drugsList = ["Capomulin","Ramicane","Infubinol","Ceftamin"]
tempVol = []

print(f'Quartile Data on Final Tumor Volumes')
print()
for x in drugsList:
    tempVolValue = finalTumordf.loc[finalTumordf["Drug Regimen"] == x, "Final Tumor Volume"]
    tempVol.append(tempVolValue)
    print(f'{x}')
    print(f'--------------------------------')
    qrtl = tempVolValue.quantile([.25,.5,.75])
    lowerQrtl = qrtl[0.25]
    print(f'Lower Quartile: {"%.3f" % lowerQrtl}')
    upperQrtl = qrtl[0.75]
    print(f'Upper Quartile: {"%.3f" % upperQrtl}')
    iqrQrtl = upperQrtl - lowerQrtl
    print(f'IQR: {"%.3f" % iqrQrtl}')
    lowerBoundQrtl = lowerQrtl - (1.5 * iqrQrtl)
    print(f'Lower Bound: {"%.3f" % lowerBoundQrtl}')
    upperBoundQrtl = upperQrtl + (1.5 * iqrQrtl)
    print(f'Upper Bound: {"%.3f" % upperBoundQrtl}')
    outliersQrtl = tempVolValue.loc[(tempVolValue < lowerBoundQrtl) | (tempVolValue > upperBoundQrtl)]
    outliersQrtl
    if outliersQrtl.empty:
        print(f'No outliers for {x}')
    else:
        print(f'Outliers present for {x}')
        print(finalTumordf.loc[outliersQrtl.index,["Mouse ID", "Final Tumor Volume"]])
    print()


# In[154]:


# Plot boxplots
# Refit data from dataframes into 2D array
obsArray = [capdf["Final Tumor Volume"],ramdf["Final Tumor Volume"],infdf["Final Tumor Volume"],cefdf["Final Tumor Volume"]]

# Draw boxplot chart
fig2, ax2 = plt.subplots()

# Format boxplot chart
ax2.set_title('Final Tumor Size Per Drug Regimen')
ax2.set_ylabel('Final Tumor Volume')
red_diamond = dict(markerfacecolor='r', marker='D')
plt.xticks([], ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]) 
ax2.boxplot(obsArray, flierprops = red_diamond)
plt.show()


# In[161]:


# Select a mouse that was treated with Capomulin and generate a line plot of tumor volume vs. time point 
# for that mouse.
# Get a random mouse to sample under Capomulin regimen
sampleMousedf = mergeDataClean[mergeDataClean["Drug Regimen"] == "Capomulin"].sample(1)
sampleMouseID = sampleMousedf.iloc[0]["Mouse ID"]
capSingledf = mergeDataClean.loc[mergeDataClean["Mouse ID"] == sampleMouseID,:]
print(capSingledf[["Mouse ID", "Timepoint", "Tumor Volume (mm3)", "Drug Regimen", "Weight (g)"]])

# Plot line graph for sample mouse
capSingleline, = plt.plot(capSingledf["Timepoint"], capSingledf["Tumor Volume (mm3)"], 
                          label=f'Mouse ID: {sampleMouseID}')
# Format line graph
plt.title("Tumor Volume Over Time in Specimen Under Capomulin")
plt.xlabel("Timepoint")
plt.ylabel("Tumor Volume (mm3)")
plt.xlim(0, max(capSingledf["Timepoint"]) + 2)
plt.ylim(min(capSingledf["Tumor Volume (mm3)"]) - 0.5 , max(capSingledf["Tumor Volume (mm3)"]) + 1)
plt.legend(handles=[capSingleline], loc="best")
plt.show()


# In[160]:


# Generate a scatter plot of mouse weight versus average tumor volume for the Capomulin treatment regimen.
capAlldf = mergeDataClean.loc[mergeDataClean["Drug Regimen"] == "Capomulin",:]
# Find average tumor volume per mouse
capAllGrouped = capAlldf.groupby("Mouse ID")
capAvgVol = capAllGrouped["Tumor Volume (mm3)"].mean()
capWeight = capAllGrouped["Weight (g)"].unique().astype(int)
# Plot scatter chart
plt.scatter(capWeight,capAvgVol)
# Format scatter chart
plt.title("Mouse Weight Versus Average Tumor Volume for Capomulin")
plt.xlabel("Weight")
plt.ylabel("Average Tumor Size")
plt.xlim(min(capWeight) - 0.8, max(capWeight) + 2)
plt.ylim(min(capAvgVol) - 0.5 , max(capAvgVol) + 2)
plt.show()


# In[159]:


# Calculate the correlation coefficient and linear regression model between mouse weight and 
# average tumor volume for the Capomulin treatment. Plot the linear regression model on top of 
# the previous scatter plot.
# Calculating correlation coefficient
correlation = st.pearsonr(capWeight,capAvgVol)
print(f'The correlation between Mouse Weight and Average Tumor Volume for Capomulin is: {round(correlation[0],2)}.')

# Prepare objects for drawing
(slope, intercept, rvalue, pvalue, stderr) = st.linregress(capWeight, capAvgVol)
regress_values = capWeight * slope + intercept
line_eq = "y = " + str(round(slope,2)) + "x + " + str(round(intercept,2))
# Plot scatter chart
plt.scatter(capWeight,capAvgVol)
# Draw regression line and equation notation
plt.plot(capWeight,regress_values,"r-")
plt.annotate(line_eq,(21.5,40),fontsize=12,color="red")
# Format scatter chart
plt.title("Mouse Weight Versus Average Tumor Volume for Capomulin")
plt.xlabel("Weight")
plt.ylabel("Average Tumor Size")
plt.xlim(min(capWeight) - 0.8, max(capWeight) + 2)
plt.ylim(min(capAvgVol) - 0.5 , max(capAvgVol) + 2)
plt.show()


# In[ ]:





# In[ ]:




