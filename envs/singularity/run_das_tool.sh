#!/bin/bash
#To generate the das_tool.sif, run the following command: singularity build das_tool.sif docker://shengwei/das_tool

read -p "
Enter the fullpath to a folder to mount to '/mnt' within the Das_Tool container
 (default:/lustre/groups/hpc/urban_lab/tools/Das_Tool/mount ): " mountfolder

# Check if the user provided a folder name
if [ -z "$mountfolder" ]; then
  # Assign a default folder name if no input provided
  mountfolder="./mount"
fi

echo Mounting DAS_Tool container to "$mountfolder"
echo If you have trouble mounting this folder, try another one near the default directory.
apptainer shell --mount type=bind,src=$mountfolder,dst=/mnt das_tool.sif
