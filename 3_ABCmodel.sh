# extract cell-type average Hi-C data
nohup python /data/whu/software/abc/ABC-Enhancer-Gene-Prediction-1.1.2/workflow/scripts/extract_avg_hic.py --avg_hic_bed_file /data/whu/abc_model_project/ref/average_hic/ENCFF134PUN.bed.gz --output_dir /data/whu/abc_model_project/ref/average_hic &

# running
mamba activate abc-env
nohup snakemake -j1 > cell_type.log 2>&1 &


