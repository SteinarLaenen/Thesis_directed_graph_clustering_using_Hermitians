Python version used during this project 2.7.12

All package requirements are listed in the requirements.txt file. 

This folder contains all the code needed to reproduce the results of Steinar Laenens AI MSc thesis "Clustering Directed Graphs Using Hermitian Laplacians" - Can be found here: https://project-archive.inf.ed.ac.uk/msc/20193694/msc_proj.pdf

There are two folders containing all the files. One is the theory folder, which contains some simple scripts to run theory experiments on. 

The folder to reproduce the experimental section is ./UN_trade

The ./UN_trade/main.py file contains all the code written to perform experiments on the UN Comtrade database. For the interested reader, we recommend them to go through this code if they want to reproduce results and understand what the functions do. This file is copied in the plot_configurations with different parameters and plotting settings. Most functions in the code are to ensure easy plotting of the figures, and thus relate to tracking averages, maintaining results and maintaining objective function scores. The spectral clustering algorithm itself only is a few lines long.

The two different plots the files can output are the line plots, and these are outputted in the "plot_cut_imbalance_plots()" function.

Plots can also be visualised in the "plot_top_clusters()" function.

The adjacency matrices are created in the "create_adjacency_matrix" function, and can output Hermitian adjacency matrices or regular ones.

To plot the figures in the thesis, we have seperate files in the ./UN_trade/plot_configurations folder, which need to be moved up one directory (mv ./UN_trade/plot_configurations/__wantedplot__.py ..) and afterwards they can be run using (python __wantedplot__.py) to reproduce the corresponding figures.

The other files in the UN_trade folder are auxillary files to perform spectral clustering and spectral clustering using the Hermitian adjacency matrix. The code written in these files is taken from this GitHub repository:

https://github.com/Jonas1312/CommunityDetection

I recommend going through their documentation if you are looking for Python
code to do spectral clustering using Stochastic Block Models


Which files correspond to which figure to plot:

./code/UN_trade/cut_imbalance_plots.py plots Figure 4.3
./code/UN_trade/HS27_main_results_all_appendix  plot all figures from appendix B.2
./code/UN_trade/HS27_main_results_HERMRW.py plots Figures 4.4 and 4.5
./code/UN_trade/HS4_effect_of_k plots Figure 4.7
./code/UN_trade/k=2_wood_plot_for_k2_section_appendix.py plots Figure B.1
./code/UN_trade/k=2_wood_plot_for_k2_section_mainsection.py plots Figure 4.2
./code/UN_trade/randomsampling.py plots figure 4.8

For any questions, please contact steinar9@gmail.com 
