clear
clc

mex ./utilities/convc.c -outdir ./utilities
mex ./utilities/matchTemplate.c -outdir ./utilities
mex ./utilities/mkCumDist.c -outdir ./utilities
mex ./utilities/mkCumProd.c -outdir ./utilities
mex ./utilities/xcorr21.c -outdir ./utilities
mex ./utilities/xcorr212.c -outdir ./utilities
mex ./utilities/xseqc.c -outdir ./utilities