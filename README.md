# Bioinfo Project 2023
## UPDATE 4 APRIL 2023
> Add more explanation in STEP 5.

> Add code in STEP 5 for T-test ALS vs Control samples.

> Add STEP 6 (Elastic-Net) and a new dataset

## UPDATE 29 March 2023 Bis
> Add PCA normalization explanation and code.

## UPDATE 29 March 2023
> Reorder the "descriptive analysis" step and add explanation in the "sample description" section.

## UPDATE 22 March 2023
>Adding examples in Step 1 "make your first functions".

## UPDATE 16 March 2023
>Adding code example for Step 1. 

This repository contains materials for the bioinfomatics project (2023 L3 INFO Paris-Saclay).
The goal of this project is to identify biomarkers associated with ALS disease. 
To do this, you have access to RNA-Seq sequencing data from post-mortem brain cortex biopsies from individuals who had ALS and others who did not.

The data for this project comes from the study "Postmortem Cortex Samples Identify Distinct Molecular Subtypes of ALS: Retrotransposon Activation, Oxidative Stress, and Activated Glia" by Tam et al. The full study can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6866666/).

## Introduction

This README will guide you throughout the project and your analyses. All steps below must be completed. Unless explicitly stated otherwise in the document (E.g., the mandaotry instructions), you have full freedom regarding the different options available to accomplish all tasks, especially concerning the coding.
This document will be updated at each session. 
Refer to the update date above.

You will work in pairs, from the begining. 
You are allowed to discuss with other pairs, but no direct code sharing (copy/paste etc). 
During all the project, note precisely your personnal contribution and the contribution of your pair.
You be asked for this during the last session (~15 min oral presentation). 

During each of the first sessions, you will be given lessons/introductions on notions you may not have yet. 
Once this lessons start, stop all of your activities and follow the lessons, even if you are familar with the concept.

At the end of your project, you must send by email (phd.rinaudo@gmail.com):
- All your code,
- A report containing the result of your analyses,
- A text file (.csv, .txt...) contaning the ordered list of your top 100 genes (first = better).

The email must have the subject: "Projet_bioinfo_2023_nom1_prenom1_nom2_prenom2" and contains an archive named: "nom1_prenom1_nom2_prenom2" with all the listed requirements.

Note: It is possible to mix code and report in a notebook, and to mix notebook and independent code. 
However, the report must be contained in a single document.
The criteria for evaluation will be (in order):
- Compliance with instructions,
- Clarity of code,
- Clarity of the report,
- Relevance of the code,
- Results found.

## Mandatory general instructions:
- The programming language used must be Python (and only Python),
- No dependencies other than Python modules should be used (Jupiter notebooks accepted),
- Unzipping your archive will produce a folder name "nom_prenom", this is "your folder",
- The code should run once your folder is placed in this repository (i.e. use relative paths only!!),
- The code should not produce any file outside of your folder,
- The code must be commented,
- The code must follow PEP 8 as much as possible,
- The code must contain unit tests,
- The code must be object-oriented.

Note: the objective is to perform an analysis. Therefore, you should not (necessarily) aimed to produce a bioinformatic pipeline robust to changes in the data format for example.
However, it is expected that you ensure that all analyses are correct, in a way or in another.

## Final note:
You will find in the following several steps to help you for the analysis of the data. 
Some steps will tell you exactly what to do, and some will be very open, requiring an investigation and a reflexion. 
Before coding an analysis you could have imagined, present it to me so that I can advise you. 
Also do not rush the avalaible steps, follow the cadence of the sessions. 
For difficult steps, I will initiate and conduct the reflexion. 

# Step 1 - Data Preprocessing
You have access to the raw data of the study in the folder "Data".
This folder contains two information:
- The RNA counts of each sample,
- An annotation file, giving information on the experiment, and (more interestingly) the samples.

Download the data and start to preprocess them.

## Gather RNA counts
To analyze the samples, you will need to merge them into a single Python object.
One standard way to do this is to build a dataframe (or any "table like" struture) such as each row is "sample" and each column is a "gene". Make sure to test your dataset, so that if you change something later on, errors can be catch easily.

Don't forget that the code should be object-oriented.

You can fit your code on the exact state of the data, i.e., design your code to work on thise data on not necessarily on other. 
Here is an example on how you can load the data
```python
import pandas as pd
import glob
import re

path = "./Data" # the path of the data

pdList = [] # variable to temporary store all dataframes (one for each txt file)
# For all txt file
for fname in glob.glob(path+"/*.txt"):
    df = pd.read_table(fname) # put the file in a dataframe
    sample_name = re.search("GSM\d+", fname).group() # search the name of the sample in the file name
    df.rename(index= df["gene/TE"], inplace=True) # rename the index (=rows) using the column containing the gene name
    df.drop(columns=df.columns[0], axis=1, inplace=True) # drop the first column containing the gene name, no more need
    df.rename(columns={ df.columns[0]: sample_name }, inplace = True) # rename the column (there is only one column at this step) using the sample name
    pdList.append(df) # add the current dataframe in the list
data_matrix = pd.concat(pdList, 1) # concat all dataframe in 1 dataframe
data_matrix = data_matrix.transpose() # transpose the dataframe to get a more standard shape (samples x variables)
```
This code should be included in an oriented object code. 
Also, this code assummes that each file is correct, contains no error etc... 
This is something to check to be rigorous. 

## Gather sample annotations
The sample annotations are all placed in a unique "xml" file. First, open the file with any text editor, and try to understand its architecture. Then, identify the information that could be relevant for your analysis. 

Finally, create a dataframe (or any "table like" structure) such as each row is a samples and each column an annotation. Make sur to test your dataset (gene counts+annotations) so that you can catch any error in the next steps.

To parse an xml, you can use the library "xml.etree.ElementTree".
You need to explore the file manually (using a text editor) to catch the structure and the name of all blocks etc...
The samples are contained in blocks named "Sample", and other information are in other blocks that you need to identify.
Here is an example to make a dataframe containing only on column corresponding to the "Cns_subregion".

```python
data_annotation = pd.DataFrame(columns = ['Sample_id', 'Cns_subregion']) # initialisation of the dataframe
xtree = et.parse('./Data/GSE124439_family.xml') # create a variable containing the xml in a tree shape
xroot = xtree.getroot() # get the root of the tree to start the exploration of the tree/xml
# for each element named "sample" that can be found from the root
for child in xroot.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Sample"):
    temp_sample_id = child.attrib['iid'] # the attribut of this node contains the sample id ()
    # for each element named "Characteristics" that can be found from the current sample
    for child2 in child.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Characteristics"):
        if(child2.attrib["tag"] == "cns subregion"):
            temp_cns_subregion = child2.text.replace('\n', '')
    temp_df = pd.DataFrame({'Sample_id': [temp_sample_id], 'Cns_subregion': [temp_cns_subregion]})
    data_annotation = pd.concat([data_annotation, temp_df])
```

## Make first preprocessing functions
By now, you should already have at least one class in your code, with associated getters and setters.
By getters and setters I mean methods that should be used to access to the attributes of instances. 
Good pratice is to prevent the direct access of the instance attributes, so that you can control the way they can be update. 
To do this in python, you can use the prefixe "__" before the attribute name. 
```python
class ALS_RNAseq:
    def __init__(self):
        self.__data_matrix = []
```
This way, the attribute "__data_matrix" can not be change or even read outside of the class. 
A getter (or setter) method need to be developped for this purpose. 
```python
class ALS_RNAseq:
    def __init__(self):
        self.__data_matrix = []

    def get_data_matrix(self):
        return self.__data_matrix
```
In this example, the getter method have no real interest, but anyway it is good pratice. 
Think about other preprocessing functions that could usefull for the next steps. 
For example, functions that can check if you have all needed annotations, or function that can subset your data based on some annotation criteria (i.e., get the subdataframe the "control" samples only). Or functions that modified one attribute (e.g., the data matrix) and update accordingly the other attributes (e.g., the annotations).

At last, you may want a convient way to "check" or look at your objects. 
Therefore a "print" function can be defined. 
By simply defining a "__str__" method in your class you can change the behevior of the "print" function when used with your class objects. 
```python
class ALS_RNAseq:
    def __init__(self):
        self.__data_matrix = []

    def get_data_matrix(self):
        return self.__data_matrix
    
    def __str__(self):
        return "put here what you want to print"
```

# Step 2 - Descriptive analysis
The descriptive analysis covers all kind of analyses that are direct description of the data, such as computing mean, standard deviation, histogram, boxplot...

## Samples description:
For each sample, compute the mean (across all genes), the median, the standard deviation. Find a way to efficiently report all those data.

Use the annotation data to describe your whole dataset. How many "disease groups", how many sample "sources", how many samples per individual etc... 
This description should guide you for the next step, in order to either correctly group your data when comparing subsets of samples or to avoid potential bias.

You should output all information relative to the samples concisely, so that we have a full (but summarize) view of our sample.

## RNA counts description:
For each gene, compute the mean (across all samples), the median and the standard deviation. Find a way to efficiently report all those data, and make your first interpretation. *Spoiler* : you may have to use graphs and transform/manipulate the data.

Samples correspond to different individuals, to ALS or control individuals etc... Think about what kind of subsets you could analyse and why (and do the descriptive analysis for those subsets).

## Start your report
At this step, you should have already begin your report. 
Before going futher, clean your code and refine your report.

# Step 3 - PCA
As you may have observed, the number of genes is far too high to compare all samples using all genes with simple analyses. 
The PCA is a classical first step analysis in those cases, and offers (among other things) a good way to visualize your data.
To understand what a PCA is, let's check at my favorite youtube channel:
[StatQuest: PCA Step-by-Step.](https://www.youtube.com/watch?v=FgakZw6K1QQ) 
We will review the video together, wait for me please.

To implement a PCA in python, a simple way is to use the PCA function in the sklearn.decomposition package. 
Scikit-learn is a wonderfull Python library, and contains a lot of "must-have" features needed for a data-scientist. 
Take some time to visite the official [website.](https://scikit-learn.org/stable/)
For a pratical python PCA tutorial, let's check again a [Josh Starmer's video](https://www.youtube.com/watch?v=Lsue2gEM9D0&ab_channel=StatQuestwithJoshStarmer).

Before doing a PCA on your data, it is mandatory to "center" your data (so that each gene has a mean equal to zero) and it is advise to "scale" your data (so that each gene has a standard deviation of 1).
This way, you ensure that the PCA is accuratly done (with the centering) and all genes are considered base on their relative variability (with the scaling) and not their absolute value. It is not advise to consider multiple variables with different scale all together.

You can use the following code before doing a PCA (X will be used in the PCA after this code)

```python
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X = scaler.fit_transform(my_data) # my_data being your dataframe containing your genes in column and samples in row
```

Now, perform a PCA and plot the samples using their coordinates in the first PCs. 
TIPs: to select the good number of PCs, compute the percenatage of variance their capture.
Use the annotations to color your plots, and see if you can already observe some kind of signal.

(Bonus) PCA is also good way to find outliers. 
Outliers are samples that are greatly different from the other samples. 
The difference should be "huge", so that only experimental errors could explain it.
Using the PCA and visualization, look at possible outliers.

# Step 4 - tSNE and UMAP (optional)
Another (more recent) good vizualization tool for high dimensional data is the [t-SNE](https://www.youtube.com/watch?v=NEaUSP4YerM&ab_channel=StatQuestwithJoshStarmer), and its little brother, the [UMAP](https://www.youtube.com/watch?v=eN0wFzBA4Sc&t=482s&ab_channel=StatQuestwithJoshStarmer). 
The advantage of this two methods is that they can reduce the dimension of your data using a desired number of components (2 most of the time), not (too much) leaving alway a part of your data variability (in theory). 
On the other hand, they do not preserve large distance interpretation, so that only "local similarities" must be interpreted (e.g., outliers are much more difficult to spot). 
UMAP tends to preserve much better large distances, but still not reach the PCA in this topic.

Try to implement a t-SNE and/or a UMAP. 
UMAP can be implemented using the "umap" module, whereas t-SNE has a scikit-learn implementation. 

Compare this visualition vs the PCA one.

# Step 5 - Univariate analysis
We have started to explore our data by computing basic statistics, making visualizations etc... 
Now it's time to perform more advanced analyses, making use of more advanced statistical background. 

In modern datascience, we are used to manipulate complexe algorithms, like regression, ensemble methods or even more complex methods (like deep learning algorithms). 
However, results found using those algorithms can be hard to interpret.
Therefore, standard univariate analyses are always a good idea to start an analysis. 
From all univariate analyses, the student-test, or t-test, (and in less extend the wilcoxon test) is one of the simplest and powerful method when applicable (it is in our case). 
This test rests on [hypothesis testing](https://www.youtube.com/watch?v=0oc49DyA3hU), making interpretation extremely straightforward (thanks to the sacrosanct [p-value](https://www.youtube.com/watch?v=vemZtEM63GY))

### T-test implementation
The t-test compare the means of a selected variable (here one gene) of two different groups (e.g., control vs ALS persons). 
Therefore, to investigate the data using a t-test, we need to perform one t-test per gene, resulting in as many results as genes.
To implement a t-test and other frequestist tests, [scipy.stats](https://docs.scipy.org/doc/scipy/reference/stats.html) contains a good sets of functions.

Let's say we want to compare the mean count of RNA of ALS subjects against the mean count of RNA of control subjects for the first gene in our dataset. 
We first need to keep the first column only:
```python
first_gene = data_matrix.iloc[:,0] # with data_matrix a dataframe containing the RNA count, with samples in row and genes in column 
```
Then, split the column into two series, one for ALS subjects, and one for control subjects:
```python
als_index = data_annotation['Group'] == "ALS Spectrum MND" # with data_annotation the dataframe containing the annotation
ctrl_index = data_annotation['Group'] == "Non-Neurological Control"
first_gene_als = first_gene[als_index.values]
first_gene_ctrl = first_gene[ctrl_index.values]
```
Finally we compare the two series using a t-test and we get the p-value:
```python
from scipy.stats import ttest_ind
pvalue_first_gene = ttest_ind(first_gene_als, first_gene_ctrl).pvalue
```

### T-test for each gene
Obviously we need to put that code in a loop or any code that will perform the t-test iteratively on each column/gene. 
One way to do this is to use a for loop.
Another way is to use the ".apply" method of pandas dataframe (if your data_matrix is a pandas dataframe):
```python
pvalues = data_matrix.apply(lambda x: ttest_ind(x[als_index.values], x[ctrl_index.values]).pvalue)
```

### Multiple testing
The p-values give the probability to obtain the observed the difference in mean of the RNA count if there is actually no difference in reality. 
So low p-values can be interpreted as: "unexpected observations". 
However, the more observations we make, the more likely we are to have surprising observations.
So we should account for the "multiple testing" we do.
In general, when exploring data (and not looking for a confirmation), controling the false discovery rate (fdr) is the way to go. 
The fdr "limits" the number of false positive discoveries. 
For example, a fdr of 5% tries to ensure that only 5% of our results will be false positives. 
To correct (raw) p-values using fdr, we can use the "statsmodels.stats.multitest.multipletests" function with the "method" parameters equals to "fdr_bh".

### Graphical representation
P-values are computed using 3 factors:
- the difference in mean
- the standard deviation(s)
- the number of samples in each groups

The first factor (difference in mean) is what we can call the "effect". 
The two other factors will determine some kind of confidence for the first factor (how much the effect is trustable). 
So we can have huge effects with relatively high p-values, and small effect with low p-values. In biology, we always look at big and significant effect first (i.e., big difference, low p-value).
A good way to report results is to use the [volcano plot](https://en.wikipedia.org/wiki/Volcano_plot_(statistics)). 

## STEP 6 - Multivariate Analysis - Elastic-Net
The t-test is an excellent tool because it offers a way to estimate the relevance of each gene (to distinguish ALS vs control samples for example). 
With this analysis, we are able to rank the genes from the most relevant to the less usefull for our question.
However, this analysis looks at the genes one by one, and doesn't try to take advantage of gene combinations. 

We can extend the t-test (using multiple genes simultaneously) with the help of logistic regression. 
The idea behind a logistic regression is to use the RNA counts of a sample as variables (for all genes) and try to predict the sample group (ALS or control for example).
To do so, we try to find a linear combination of the variables that compute a probability for each possible classes (ALS or control).

Regularized logistic regressions are extentions of standard logistic regressions that penalized variables with no clear interest.
In another words, they try to attribute a coefficient of 0 to variables with low interest (no useful signal).
Therefore this methods are very good candidates for biomarker selections. 

Among these regularized regressions, [elastic-net](https://www.youtube.com/watch?v=1dKRdX9bfIo&ab_channel=StatQuestwithJoshStarmer) is one of the most versatile. 
It can handle "wide" dataset (more variables than samples), correlation (do not have to choose between very similar variables) and are easily configurable (parameter tuning).

Warning: like the PCA, the elastic-net algorithm takes multiple variables into account. 
Therefore, variables with high values and high standard deviation can introduce bias (being artifically more important). 
So, like PCA, I advise you to normalize the data before using elastic-net.

(Optional) Also, you can also use the output of the PCA and perform the elastic-net on this transformed data.
Interpretations will be more difficult (you will have to look at loadings of your components at the end of the analysis) but potentially more powerfull.

### Elastic-Net implementation
You can use the method "LogisticRegression" from scikitlearn (sklearn.linear_model) to implement an elastic-net analysis:
```python
from sklearn.linear_model import LogisticRegression
elasticNet = LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, C=0.5).fit(x, y)
```
Where x contains the data and y the groups (2 groups only).
Two parameters can be fine tune: 
- l1_ratio: the ration between l1 and l2 regularizations
- C: inverse of regularization strength

### Fine tuning of parameters
To fine tune parameters while avoiding overfitting, we can perform a [cross validation](https://www.youtube.com/watch?v=fSytzGwwBVw&ab_channel=StatQuestwithJoshStarmer): 
```python
from sklearn.linear_model import LogisticRegression
elasticNet = LogisticRegressionCV(penalty='elasticnet', cv= 3, solver='saga', l1_ratios=[0.25,0.5,0.75], Cs=[0.1,0.5], scoring= 'accuracy').fit(x, y)
```
You have to specify the number of folds (here cv=3). 
Choose a number of fold adapted to the data. 
By default, the folds are stratified (they split the dataset taking the class balance into account, something that we absolutely want here).

You also have to specify the scoring function used to evaluate and compare the models (accuracy in this example). 
Check the [documentation](https://scikit-learn.org/stable/modules/classes.html#module-sklearn.metrics) to get a list of scoring functions and choose one adapted to the our data.
Have a special look to [balanced accuracy](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.balanced_accuracy_score.html#sklearn.metrics.balanced_accuracy_score), [f1 score](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html#sklearn.metrics.f1_score), [ROC AUC](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html#sklearn.metrics.roc_auc_score), and [the matthews correlation coefficient](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.matthews_corrcoef.html#sklearn.metrics.matthews_corrcoef).

Once you have your elastic-net model (with the best parameters), you can look at its overall performance.
```python
model_prediction = elasticNet.predict(x)
accuracy_train_dataset = sklearn.metrics.accuracy_score(y, model_prediction)
```
Here I used the accuracy metric, but remember to use (at least) the metric you used during the tuning.

### Testing your model
We built and evaluated our model on the same dataset. 
Even if we used a cross-validation to choose the parameters, the final model used all the avalaible data. 
Therefore we have no idea how the model performs on new data. 
For this purpose, we need a "test set". 

Download the "data set", process the data (warning: there might be some differences in the xml file) and use this new dataset to evaluate the model.
I advise you to use different scores (begining with the one used to choose the parameters).

### Get the best genes
The elastic-net model contains (by definition) all the coefficients of the regression, i.e., the importance of each gene to differentiate the ALS samples from the control samples.
```python
elasticNet.coef_
```
You can rank these coefficients (using their absolute values) to obtain a good list of candidates.





