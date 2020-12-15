# Flexible Cells Algorithm
Consider a set of geo-tagged events like location of COVID-19 cases and location of car accidents. One important thing about these events is to tell where is the center of this events? By center we mean the location with maximum number of events. For example in case of COVID-19, in which country, city, or even smaller areas it starts to spread to other places? Another important thing about these events is that how they spread from the center to other locations?
So, in this project we try to find centers of a set of geo-tagged events and their spead/dispersion.

## FCA
This folder contains our algorithm's implementation on finding multiple centers. It is much like BFS algorithm plus early termination of bigger cells. 
It has three main phases
1) Searching Phase: finding candidate cells using BFS and early termination
2) Checking Phase: removing false positives from results of searching phase using brute-force checking of all possible cases
3) Parameter Tuning Phase: using Backstrom's probabilistic framework

## SDG
This folder contains our SyntheticDataGenerator algorithm which tries to genertae some dataset with a specific number of events using uniform distribution

## BSC
This folder contains implementation of Single-Center case of Backstrom's algorithm plus a SyntheticDataGenerator which generates a dataset to show weakness of Backstrom's work

### Contributors
Hamed Mirzaei<br/>
Sarah Davis
