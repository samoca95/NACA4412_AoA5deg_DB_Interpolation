#rsync -vhra oner751450@glogin1.bsc.es:~/scratch/Nek5000_Interpolation /data2/wing --include='**.gitignore' --exclude='/.git' --filter=':- .gitignore' #--delete-after
rsync -vhra bsc:~/scratch/NACA4412_5deg_Interpolation /data2/wing --include='**.gitignore' --exclude='/.git' --filter=':- .gitignore' #--delete-after
