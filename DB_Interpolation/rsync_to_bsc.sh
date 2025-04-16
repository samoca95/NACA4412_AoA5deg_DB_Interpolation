# rsync -vhra /data2/wing/Nek5000_Interpolation oner751450@glogin1.bsc.es:~/scratch --include='**.gitignore' --exclude='/.git' --filter=':- .gitignore' #--delete-after
rsync -vhra /data2/wing/NACA4412_5deg_Interpolation bsc:~/scratch --include='**.gitignore' --exclude='/.git' --filter=':- .gitignore' #--delete-after
