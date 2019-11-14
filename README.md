# Displacement arrays of rendered tractions

Here is the code to compute and analyse displacement arrays of rendered tractions for mechanophenotyping of 3D multicellular clusters.

## Requirements:
1. Matlab (The code was tested to run in Matlab 2017b version)
2. Python 3 (Prefer anaconda installation with all the packages)
3. Jupyter notebooks

### Mathwork packages
Matlab packages used from the Mathworks file exhanchange, are in the Supporting MFiles subfolder. 

## Instruction to analysis and compute DARTs

#### Step 1: Convert images files into separate .mat files

1. Save the multi-timepoing and multi-step images in folder with date as the name of the folder.  Example:
```
project directory
│   README.md
│   ...    
│
└───supporting_mfiles
│   │   ...
|
└───20190101
    │   example.nd2
```
 
2. Now, conver the .nd2 file in the folder '20190101' into the subfolders according to the name of the well position. And in each well position folder, save individual channel for each time point. You can then delete the raw .nd2 file. At the end, the file strucute should look like this:
```
project directory
│   README.md
│   ...    
│
└───supporting_mfiles
│   │   ...
|
└───20190101
    └───A01   (Multi-well position subfolder)
    |   20190101_A01_bead_t00.mat  (Each channel saved at each multi-well position and time point as .mat file)
    |   20190101_A01_bead_t01.mat
    |   20190101_A01_bead_t02.mat
    |   20190101_A01_bead_t03.mat
    |   20190101_A01_cell_t00.mat
    |   20190101_A01_cell_t01.mat
    |   20190101_A01_cell_t02.mat
    |   20190101_A01_cell_t03.mat
    |
    └───A02
    └───A03
    └───...
```

Note, the files at timepoint t=00, refer to the image captured after cell lyses. 

Here, I decribe the process to convert multi-timepoint and multi-step .nd2 files into .mat files. You will need to adapt this section depending on your imaging format. 

3. Run `run_nd2.mat` file, to convert all .nd2 file to .mat file in the same './20190101' folder

4. Run `structure_files.ipynb`, to order the .mat files into substructures are described in step 2. Then delete the raw .nd2 file to save space in the parent directory. 

#### Step 2: Use T-PT to track particle displacement from bead channel images

Run `run_TPT.m` to run the Topology-based particle tracking algorithm to track the particle displacement from bead position. The version of TPT used is in the following subfolder `./supporting_mfiles/TPT`

If used, please cite:
[Paper: Patel, M., Leggett, S. E., Landauer, A. K., Wong, I. Y., & Franck, C. (2018). Rapid, topology-based particle tracking for high-resolution measurements of large complex 3D motion fields. Scientific reports, 8(1), 5581.](https://www.nature.com/articles/s41598-018-23488-y)
```bibtex
@article{patel2018rapid,
  title={Rapid, topology-based particle tracking for high-resolution measurements of large complex 3D motion fields},
  author={Patel, Mohak and Leggett, Susan E and Landauer, Alexander K and Wong, Ian Y and Franck, Christian},
  journal={Scientific reports},
  volume={8},
  number={1},
  pages={5581},
  year={2018},
  publisher={Nature Publishing Group}
}
```

Github link: https://github.com/FranckLab/T-PT

#### Step 3: Segment cell clusters from cell channel images

1. Filter cell images with a median filter and save it back to the same image cell. Median filtering is a computational expensive process, so we do not want to repeat this step multiple time. Hence, we save the results. Perform this step by running `run_filter_cell_img.m`.

2. First perform substep 1 and 2 from the Step 4. Then, run `run_analysis.m` file to segment the cells. Refer to the comment in line 109 in `run_analysis.m`. After finding the appropriate thresvold value, fill that value in the `thres.xlsx` file. 

#### Step 4: Compute DART and other mechanical metrics. 

1. In the cell image folder (Ex. `./20190101`), add the `Conditions.xlsx` file. In this file, save the basic information for each well like the date of experiment, well name, experimental start time, induce type, drug treatment, and drug concentration. 

2. In the same cell image folder (Ex. `./20190101`), add the `thres.xlsx` file. In this file, initiate the cell threshold to 0. 

3. Run `run_analysis.m` file to segment the cells and compute DART and other cell matrics. Remeber to comment out the section in `analysis_dart.m` between lines 116 and 132. The code saves the results for each multi-well in `dart.csv` file within each multi-well folder. 

#### Step 5: Filter cell clusters for data processing. 

Visualize the displacements for individual cell clusters and the whole cell clusters. Remove the cell clusters from the analysis based on exclusion critiera in your study (Ex. segementation of fragments of a bigger cell clusters, cell clusters sitting at the bottom of the coverslip, etc.). For all the cell clusters that you wish to keep in the analysis, open the 'dart.csv' file and mark `1` in the `is_disp_checked` column for the corresponding cluter. Only these cluster will be used in the analysis later on. 

Some helpful scripts to help you in this process:
- `vis_radial_disp`: Visualizes the radial displacement around a given cell cluster. 
- `vis_xy_cell_projection`: Visualizes the xy projection of the cell in labelled 2D image to connect the cell from the stored `dart.csv` file to the actual image. 

#### Step 6: Post-processing. 

1. Run the `consolidate_data.ipynb` to combine and save all the data from multiple `dart.csv` into a single file. 

2. Run `publicaiton_figures.ipynb` to perform the post processing, figure making and ML modeling of the data. Refer to the comments in the file for details.

## Cite
If used please cite:
[Leggett, Susan E., Mohak Patel, Thomas M. Valentin, Lena Gamboa, Amanda S. Khoo, Evelyn Kendall Williams, Christian Franck, and Ian Y. Wong. "Mechanophenotyping of 3D Multicellular Clusters using Displacement Arrays of Rendered Tractions." bioRxiv (2019): 809871.](https://www.biorxiv.org/content/10.1101/809871v1)

```bibtex
@article{leggett2019mechanophenotyping,
  title={Mechanophenotyping of 3D Multicellular Clusters using Displacement Arrays of Rendered Tractions},
  author={Leggett, Susan E and Patel, Mohak and Valentin, Thomas M and Gamboa, Lena and Khoo, Amanda S and Williams, Evelyn Kendall and Franck, Christian and Wong, Ian Y},
  journal={bioRxiv},
  pages={809871},
  year={2019},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Questions
For questions, please first refer to and [Questions/Issues](https://github.com/FranckLab/DART/issues) (make sure to look through the closed Issues too!). Add a new question if similar issue hasn't been reported. We ask you to post questions on Github, so that future users with similar issues can benefit from the answers. We shall try our best to help you at the earliest. The author's contact information can be found at [Franck Lab](https://www.franck.engr.wisc.edu/).
