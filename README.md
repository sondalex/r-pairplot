# pairplot

Library for making pairplots in R, inspired by [seaborn.pairplot](https://seaborn.pydata.org/generated/seaborn.pairplot.html) and [GGally::ggpairs()].
and based on `ggplot2` and `patchwork`.
This package is small but designed to be general. 

The grid region is made of three regions:
1) Lower triangle
2) Diagonal
3) Upper diagonal

Each region can be mapped to a plot function or NULL.
IF NULL, the element_blank() is assigned. 
This way you can build triangle plots.


