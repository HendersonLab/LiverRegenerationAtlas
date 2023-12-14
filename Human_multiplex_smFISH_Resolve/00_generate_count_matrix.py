'''
Interpolate a set of rois with a transcript coordinates file  to generate a count matrix per cell and transcript.

The script takes three arguments:
    - A *_results.txt file with transcript coordinates
    - A <roiSet>.zip file with rois to interpolate with transcripts
    - The filenames output prefix


It creates three outputs in the same folder as the roiSet.zip file:
    - A <roiSet>_countmatrix.csv - The countmatrix for the roi set
    - A <roiSet>_renamed.zip     - Renamed roifile that has usefull cell names in the format of "CellNUM_XCOORD_YCOORD"
    - A <roiSet>_cellpos_mat.csv - A position matrix with x,y coordinates and total transcript count  (potentially redundant)


Ricardo Guerreiro, Resolve Biosciences 02/09/2022
'''


from shapely.geometry import Polygon, Point, MultiPoint
import glob, roifile
import numpy as np
import pandas as pd
import argparse, os


# Parse input arguments
parser = argparse.ArgumentParser(description='''
    Takes Results.txt and roi.zip sets to create count matrix for those rois. Multiple Results.txt and roi.zip can be provided separated by commas. ''') 

parser.add_argument("coords", help="Input *_results.txt (coordinates) file")
parser.add_argument("rois",   help= "Corresponding cells or nucleii roiSet.zip files")
parser.add_argument("filenames_prefix",   help= "Filename prefix for output files")
args=parser.parse_args()

# Split arguments into lists of inputs
args.coords = args.coords.split(',')
args.rois = args.rois.split(',')
args.filenames_prefix = args.filenames_prefix.split(',')


# Main iteration
for i,slide in enumerate(args.coords):
    
    cellnum = 0 
    cell_rois = roifile.roiread(args.rois[i])

    df = pd.read_csv(slide, sep = '\t')
    df.columns = ["x","y","z","name","qual"]

    # Initialize count matrix and position matrix for cells
    count_ma = pd.DataFrame(index = df.name.unique())
    cell_pos = pd.DataFrame()
    rois_to_write = [] 
    gen_counts = []
    roi_names = []
    
    for roi in cell_rois: 
        cellnum += 1

         # Prepare roi coordinates 
        cell_pgon  = np.array(roi.coordinates()).T
        roiX,roiY = cell_pgon[0],cell_pgon[1]
        xcenter,ycenter = map(lambda coord: round(np.mean(coord)),cell_pgon)

        roi.name = f"{args.filenames_prefix[i]}_Cell{cellnum}_{xcenter}_{ycenter}"

        # Pre-filter DF for roi space only 
        df_subset =  df [((df.x < max(roiX)) * (df.x > min(roiX)))*
                         ((df.y < max(roiY)) * (df.y > min(roiY))) ]

        if len(df_subset.index)> 3:
            # Polygon intersection of rois with DF
            pgon = Polygon(roi.coordinates())
            points = MultiPoint(list(zip(df_subset.x,df_subset.y)))
            boolio = [pgon.intersects(point) for point in points.geoms]   

            df_subset = df_subset[boolio]
                
            # Prepared data that is used for count matrix
            gen_counts.append(df_subset.name.value_counts())
            roi_names.append(roi.name)

            # Cell position matrix
            cell_pos.at[roi.name,['x','y']] = xcenter,ycenter
            cell_pos.at[roi.name,'n_counts'] = int(sum(df_subset.name.value_counts()))

            rois_to_write.append(roi)
    
    # Create Countmatrix with cellnames as index,  Set NaNs to 0
    count_ma = pd.DataFrame(gen_counts)
    count_ma.index = roi_names
    count_ma [count_ma != count_ma] = 0
    count_ma = count_ma.astype(int)

    # Write count matrix and cell position matrix to same place as input roiSet.zip file 
    #outname = os.path.abspath(args.rois[i][:-4])
    outname = os.path.abspath(args.filenames_prefix[i] + "/")

    count_ma.to_csv(f"{outname}/count_mat.csv")
    cell_pos.astype(int).to_csv(f"{outname}/cellpos_mat.csv")

    #Rewrite rois with new names (QuPATH does not create unique roi.name instances)
    roifile.roiwrite(outname + '/roiset.zip', rois_to_write) 

